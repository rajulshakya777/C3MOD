import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import os
import pickle

YELLOW = "\033[1;33m"
CYAN = "\033[1;36m"
BLUE = "\033[1;34m"
BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[1;32m"
CHECK_EMOJI = f"{GREEN}✔{RESET}"

# Create output folder if it doesn't exist
output_folder = os.path.join(os.path.dirname(__file__), '..', 'output', 'immune_analysis')
os.makedirs(output_folder, exist_ok=True)

def immune_characterization(cancer_type, algorithm_choice):

    # Step 1: Load the mutation file and transpose it
    stage_file = os.path.join("data", "input_data", "TCGA_data", cancer_type, "immune_data", "Scores_160_Signatures.tsv")
    stage_df = pd.read_csv(stage_file, sep='\t')

    # Load the data and transpose it
    stage_df = pd.read_csv(stage_file, sep='\t', header=0).transpose()

    # Drop the first row (which might have been incorrectly treated as data)
    stage_df = stage_df.drop(stage_df.index[0])

    # Set the first row as the new header after transposing
    stage_df.columns = stage_df.iloc[0]  # The first row becomes the header
    stage_df = stage_df.drop(stage_df.index[0])  # Drop the now-redundant first row
    # Display the DataFrame's head to verify the result

    # Function to modify index
    def modify_index(id):
        # Remove the last part after the second-to-last hyphen
        base_id = '-'.join(id.split('-')[:4])
        # Replace hyphens with dots
        modified_id = base_id.replace('-', '.')
        return modified_id

    # Apply the function to the DataFrame's index
    stage_df.index = stage_df.index.map(modify_index)

    classification_file = os.path.join("output", "clustering_results", f"{cancer_type}_classification_{algorithm_choice}.txt")
    classification_df = pd.read_csv(classification_file, sep='\t').drop(columns=['Unnamed: 0'])
    classification_df = classification_df.reset_index(drop=True)
    merged_df = classification_df.merge(stage_df, left_on='samples', right_index=True, how='inner')
    merged_df.head()


    # Assuming merged_df is your DataFrame

    # Get the value of K from the column
    K = merged_df['K'].iloc[0]  # Assuming K is consistent across the DataFrame

    # Initialize DataFrames to store the R-values and p-values
    r_value_matrix = pd.DataFrame(index=range(1, K + 1), columns=range(1, K + 1))
    p_value_matrix = pd.DataFrame(index=range(1, K + 1), columns=range(1, K + 1))

    # Loop over each pair of clusters
    for i in range(1, K):
        for j in range(i + 1, K + 1):
            # Filter the data for the two clusters
            cluster_i = merged_df[merged_df['cluster'] == i].drop(columns=['samples', 'K', 'cluster'])
            cluster_j = merged_df[merged_df['cluster'] == j].drop(columns=['samples', 'K', 'cluster'])

            # Aggregate data by calculating the mean for each feature within the cluster
            # cluster_i_mean = cluster_i.mean()
            # cluster_j_mean = cluster_j.mean()

            # Convert to numeric and drop any non-numeric data
            cluster_i_mean = pd.to_numeric(cluster_i.mean(), errors='coerce').dropna()
            cluster_j_mean = pd.to_numeric(cluster_j.mean(), errors='coerce').dropna()

            # Check if both series have the same length after cleanup
            if len(cluster_i_mean) != len(cluster_j_mean):
                raise ValueError("Cluster means have different lengths after removing non-numeric values.")

            # Calculate Pearson correlation between the mean feature values of the two clusters
            r_value, p_value = pearsonr(cluster_i_mean, cluster_j_mean)

            # Store the R-value and p-value in the matrices
            r_value_matrix.loc[i, j] = r_value
            r_value_matrix.loc[j, i] = r_value  # Matrix is symmetric
            p_value_matrix.loc[i, j] = p_value
            p_value_matrix.loc[j, i] = p_value  # Matrix is symmetric

    # Convert R-value matrix to float type for heatmap
    r_value_matrix = r_value_matrix.astype(float)

    # Create mask for upper triangle (excluding diagonal)
    mask = np.triu(np.ones_like(r_value_matrix, dtype=bool), k=1)

    # Calculate the min and max values for dynamic scaling
    r_min = r_value_matrix.min().min()
    r_max = r_value_matrix.max().max()

    # Functions to calculate nearest lower and upper limits
    def round_down(value, decimal_places):
        factor = 10 ** decimal_places
        return np.floor(value * factor) / factor

    def round_up(value, decimal_places):
        factor = 10 ** decimal_places
        return np.ceil(value * factor) / factor

    # Define thresholds
    lower_decimal_places = 2
    upper_decimal_places = 2

    # Calculate lower and upper limits based on the thresholds
    lower_limit = round_down(r_min, lower_decimal_places)
    upper_limit = round_up(r_max, upper_decimal_places)

    # Create a masked array for the heatmap
    masked_array = np.ma.masked_where(mask, r_value_matrix)

    # Create a heatmap using matplotlib
    plt.figure(figsize=(10, 8))
    heatmap = plt.imshow(masked_array, cmap='viridis', interpolation='nearest', vmin=0, vmax=1)

    # Add colorbar with ticks formatted to three decimal places
    cbar = plt.colorbar(heatmap, format='%.3f')
    cbar.set_label('R-value')

    # Format colorbar ticks
    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

    # Add labels
    plt.xticks(ticks=np.arange(len(r_value_matrix.columns)), labels=r_value_matrix.columns)
    plt.yticks(ticks=np.arange(len(r_value_matrix.index)), labels=r_value_matrix.index)

    # Add annotations with R-values and exact p-values
    for i in range(len(r_value_matrix.index)):
        for j in range(len(r_value_matrix.columns)):
            if not mask[i, j]:  # Only annotate non-masked values
                r_val = r_value_matrix.iloc[i, j]
                p_val = p_value_matrix.iloc[i, j]
                # Format p-value with scientific notation if very small
                if p_val < 0.001:
                    p_val_text = f'{p_val:.3e}'  # Scientific notation
                else:
                    p_val_text = f'{p_val:.3f}'  # Decimal notation
                plt.text(j, i, f'R-val: {r_val:.3f}\n(p-val: {p_val_text})', ha='center', va='center', color='black')

    # Overlay light grey rectangles on the diagonal cells
    for i in range(len(r_value_matrix.index)):
        plt.gca().add_patch(patches.Rectangle((i - 0.5, i - 0.5), 1, 1, color='lightgrey', zorder=10))

    plt.title('Similarity between clusters (ImmuneScores-based)')
    plt.xlabel('Cluster')
    plt.ylabel('Cluster')

    plot_file = os.path.join("output", "immune_analysis", f"{cancer_type}_{algorithm_choice}_immune_scores_based_similarity_between_clusters_heatmap.png")
    plt.savefig(plot_file)


def main():
    with open(os.path.join('.', 'modules', 'data.pkl'), 'rb') as f:
        X_scaled, pca_df, cancer_type, k, algorithm_choice = pickle.load(f)  # Load the entire object

    # Create output folder if it doesn't exist
    # print("Performing mutation analysis on cancer Type:", cancer_type)

    algorithms = ["SNF", "KMeans", "Hierarchical", "SpectralClustering", "FuzzyCMeans", "All"]
    algorithm_choice = int(algorithm_choice)

    if algorithm_choice == 6:
        print(f"{BOLD}🧬 Running Immune Analysis for Clusters: {RESET} SNF, KMeans, Hierarchical, SpectralClustering, and FuzzyCMeans")
        for i, algorithm in enumerate(algorithms[:-1], 1):  # Iterate over all algorithms except "All"
            immune_characterization(cancer_type, algorithm)
    else:
        algorithm = algorithms[algorithm_choice-1]
        print(f"{BOLD}🧬 Running Immune Analysis for Cluster: {RESET}{algorithm}")
        immune_characterization(cancer_type, algorithm)

    print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW} All Immune Analysis results saved to immune_analysis folder.{RESET}")

# Call the main function
if __name__ == "__main__":
    main()
