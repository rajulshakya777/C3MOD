import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
import warnings

# Suppress pandas SettingWithCopyWarning
pd.options.mode.chained_assignment = None

# Suppress all FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Suppress all UserWarnings
warnings.simplefilter(action='ignore', category=UserWarning)

# Suppress any other warnings (e.g., all warnings)
warnings.filterwarnings('ignore')
YELLOW = "\033[1;33m"
CYAN = "\033[1;36m"
BLUE = "\033[1;34m"
BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[1;32m"
CHECK_EMOJI = f"{GREEN}✔{RESET}"

# Create output folder if it doesn't exist
output_folder = os.path.join(os.path.dirname(__file__), '..', 'output', 'stage_analysis')
os.makedirs(output_folder, exist_ok=True)

def stage_analysis(cancer_type, algorithm_choice):
    # Loading stage data and classification data to find which patient is in which cancer stage (Used common samples which)
    # For which we are having the stage information

    # Step 1: Load the mutation file and transpose it
    stage_file = os.path.join("data", "input_data", "TCGA_data", cancer_type, "stage_data", f"clinical_patient_{cancer_type}_info.tsv")
    stage_df = pd.read_csv(stage_file, sep='\t')

    # Filter the data to include only rows where the tumor stage is not empty
    stage_data = stage_df[[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode', f'clinical_patient_{cancer_type.lower()}.ajcc_pathologic_tumor_stage']]
    stage_data = stage_data[stage_data[f'clinical_patient_{cancer_type.lower()}.ajcc_pathologic_tumor_stage'].notna()]

    # Define a dictionary to map Roman numerals to Arabic numerals
    roman_to_arabic = {
        'I': 1,
        'II': 2,
        'III': 3,
        'IV': 4,
        'V': 5,
        'VI': 6,
        'VII': 7,
        'VIII': 8,
        'IX': 9,
        'X': 10
    }

    # Function to convert Roman numerals to Arabic numerals, ignoring the last alphabet if present
    def convert_stage(stage):
        stage_roman = stage.rstrip('ABCDEFG')
        return roman_to_arabic.get(stage_roman, stage_roman)

    # Apply the conversion to the stage column
    stage_data['ajcc_pathologic_tumor_stage_numeric'] = stage_data[f'clinical_patient_{cancer_type.lower()}.ajcc_pathologic_tumor_stage'].apply(convert_stage)

    # Drop the rows with headers or metadata if they exist
    stage_data = stage_data[stage_data[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode'] != 'bcr_patient_barcode']

    # Select the patient barcode and the numeric stage
    result_data = stage_data[[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode', 'ajcc_pathologic_tumor_stage_numeric']]
    
    # Define the mapping
    stage_mapping = {
        'Stage I': 1,
        'Stage II': 2,
        'Stage III': 3,
        'Stage IV': 4
    }
    result_data['ajcc_pathologic_tumor_stage_numeric'] = result_data['ajcc_pathologic_tumor_stage_numeric'].replace(stage_mapping)

    # Display the updated data
    result_data = result_data.reset_index(drop=True)
    result_data[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode'] = result_data[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode'].astype(str)

    result_data[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode'] = (
        result_data[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode']
        .str.replace('-', '.')
        .apply(lambda x: x if x.endswith('A') else f"{x}.01A")
    )
    stage_df = result_data
    stage_df = stage_df[stage_df['ajcc_pathologic_tumor_stage_numeric'].isin([1, 2, 3, 4])]
    stage_df.head()


    classification_df = pd.read_csv(os.path.join("output", "clustering_results", f"{cancer_type}_classification_{algorithm_choice}.txt"), sep='\t')


    classification_df['sample_id'] = classification_df['samples']
    stage_df['sample_id'] = stage_df[f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode']
    common_samples_df = pd.merge(classification_df, stage_df, on='sample_id', how='inner')
    common_samples_df = common_samples_df.drop(columns=['sample_id', f'clinical_patient_{cancer_type.lower()}.bcr_patient_barcode'])

    # Determine the number of clusters (value of K)
    K = common_samples_df['K'].iloc[0]

    # Ensure the 'ajcc_pathologic_tumor_stage_numeric' column is numeric
    common_samples_df['ajcc_pathologic_tumor_stage_numeric'] = pd.to_numeric(common_samples_df['ajcc_pathologic_tumor_stage_numeric'], errors='coerce')

    # Drop rows where the conversion failed
    common_samples_df = common_samples_df.dropna(subset=['ajcc_pathologic_tumor_stage_numeric'])
    common_samples_df.head()

    # Create a figure with subplots for each cluster
    fig, axes = plt.subplots(1, K, figsize=(15, 5), sharey=True)

    fig.suptitle(f'{algorithm_choice} : Stage Distribution Across Different Clusters', fontsize=14)

    # Create a mapping from numeric stages to Roman numerals
    stage_mapping = {1: 'I', 2: 'II', 3: 'III', 4: 'IV'}

    # Initialize list to collect wedges for the global legend
    global_wedges = None  # We will set this only once
    stage_colors = {
        'I': '#FF9999',    # Light Red (Stage I)
        'II': '#FFCC99',   # Light Orange (Stage II)
        'III': '#99FF99',  # Light Green (Stage III)
        'IV': '#9999FF'    # Light Blue (Stage IV)
    }

    for cluster in range(1, K + 1):
        # Filter data for the current cluster
        cluster_data = common_samples_df[common_samples_df['cluster'] == cluster]
        
        # Count the number of patients in each stage
        stage_counts = cluster_data['ajcc_pathologic_tumor_stage_numeric'].value_counts().sort_index()
        
        # Map numeric stages to Roman numerals
        stage_counts.index = stage_counts.index.map(stage_mapping)
        
        # Get the colors for the current stage counts
        colors = [stage_colors[stage] for stage in stage_counts.index]
        
        # Plot pie chart
        wedges, texts, autotexts = axes[cluster - 1].pie(
            stage_counts,
            labels=stage_counts.index,
            autopct='%1.1f%%',
            startangle=140,
            colors=colors  # Use the mapped colors from stage_colors
        )
        
        # Add title with total samples
        total_samples = cluster_data.shape[0]
        axes[cluster - 1].set_title(f'Cluster {cluster} (Total samples: {total_samples})')

        # Store wedges for global legend only once
        if global_wedges is None:
            global_wedges = wedges

    # Create a global legend with fixed labels after all pie charts are plotted
    axes[0].legend(
        handles=global_wedges,
        labels=['Stage I', 'Stage II', 'Stage III', 'Stage IV'],
        title="Stages",
        loc="best",
        bbox_to_anchor=(1, 1),
        bbox_transform=plt.gcf().transFigure
    )

    # Adjust layout to prevent overlap
    # plt.tight_layout()
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust to make space for the global titl
    plot_file = os.path.join("output", "stage_analysis", f"{cancer_type}_{algorithm_choice}_stage_pi_chart_clusters.png")
    plt.savefig(plot_file)

def main():
    with open(os.path.join('.', 'modules', 'data.pkl'), 'rb') as f:
        X_scaled, pca_df, cancer_type, k, algorithm_choice = pickle.load(f)  # Load the entire object

    algorithms = ["SNF", "KMeans", "Hierarchical", "SpectralClustering", "FuzzyCMeans", "All"]
    algorithm_choice = int(algorithm_choice)

    if algorithm_choice == 6:
        print(f"{BOLD}🧬 Running stage analysis for Clusters: SNF, KMeans, Hierarchical, SpectralClustering, and FuzzyCMeans{RESET}")
        for i, algorithm in enumerate(algorithms[:-1], 1):  # Iterate over all algorithms except "All"
            stage_analysis(cancer_type, algorithm)
    else:
        algorithm = algorithms[algorithm_choice-1]
        print(f"{BOLD}🧬 Running stage analysis for Cluster: {RESET}{algorithm}")
        stage_analysis(cancer_type, algorithm)

    print("Analysis Summary:")
    print("Stage analysis resutls saved to stage_analysis folder")
# Call the main function
if __name__ == "__main__":
    main()
