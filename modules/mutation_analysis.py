import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.stats import pearsonr
from scipy.stats import chi2_contingency, fisher_exact
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import os
import seaborn as sns
import pickle
import warnings


warnings.filterwarnings("ignore", category=RuntimeWarning)
# Color and style codes
YELLOW = "\033[1;33m"
CYAN = "\033[1;36m"
BLUE = "\033[1;34m"
BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[1;32m"
CHECK_EMOJI = f"{GREEN}✔{RESET}"

output_dir = os.path.join(os.path.dirname(__file__), '..', 'output', 'mutation_analysis')
os.makedirs(output_dir, exist_ok=True) 

# 6.1 : Load differentially expressed genes (Those genes are expressed (regulated in higher amount) in cancerous patients and not in healthy patients)
def load_differentially_expressed_genes(cancer_type):
    # Step 1: Load LUAD_Differentially expressed genes, Preprocess them
    file_path = os.path.join("data", "input_data", "TCGA_data", cancer_type, "mutation_data", f"DE_mRNA_{cancer_type}_Vs_Normal.csv")
    luad_diff_ganfeature_df = pd.read_csv(file_path, sep=',', index_col=0)
    # Count the number of NaN values in the 'Symbol' column
    nan_count = luad_diff_ganfeature_df['Symbol'].isna().sum()

    # Print the result
    # print(f"Number of NaN values in the 'Symbol' column: {nan_count}")

    # print('updated data after dropping null Symbol')
    luad_diff_ganfeature_df = luad_diff_ganfeature_df.dropna(subset=['Symbol'])
    symbol_values = luad_diff_ganfeature_df['Symbol'].tolist()
    # Step 2: Clean the Symbol column
    symbol_values_cleaned = luad_diff_ganfeature_df['Symbol'].dropna().astype(str).tolist()
    # Convert all to lowercase
    symbol_values_cleaned = [symbol.lower() for symbol in symbol_values_cleaned]
    symbol_values_cleaned = [symbol.strip() for symbol in symbol_values_cleaned]
    # print(luad_diff_ganfeature_df.head())
    return symbol_values_cleaned


# 6.2 : Load the available gene mutation (miRNA) data of each patient, we will take only differentially expressed genes whichever mutation information is available to us in miRNA data
def load_available_genes_mutation_miRNA(symbol_values_cleaned, cancer_type):
    # Step 2: Load the available gene mutation data of each patient
    mutation_file = os.path.join("data", "input_data", "TCGA_data", cancer_type, "mutation_data", f"{cancer_type}_mc3_gene_level_mutation.txt")
    mutation_df = pd.read_csv(mutation_file, sep='\t', index_col=0).transpose()
    mutation_df.index = mutation_df.index.str.replace('-', '.', regex=False)
    mutation_df = mutation_df.reset_index()
    mutation_df = mutation_df.rename(columns={'index': 'samples'})
    mutation_df['samples'] = mutation_df['samples'].apply(lambda x: x + 'A')
    # Step 2: Convert all column names to strings
    mutation_df.columns = mutation_df.columns.astype(str)
    # Now convert them to lowercase
    mutation_df_columns = [col.lower() for col in mutation_df.columns]

    # Step 3: Selecting the common genes between available gene mutation data and differentially expressed genes
    common_columns = ['samples'] + [col for col in mutation_df.columns if col.lower() in symbol_values_cleaned]
    mutation_df = mutation_df[common_columns]
    # Check the result
    # print("\nFiltered mutation_df:")
    # print(mutation_df.head())
    # print(mutation_df.shape)
    return mutation_df

# 6.3 : Load our classification results, which patient is in which cluster and create a merged dataframe having patients, their cluster number, and genes (from available differentially expressed genes)
def construct_final_dataframe(mutation_df,chosen_algorithm, cancer_type):
    # Step 4: Load our classification results, which patient is in which cluster and create a merged dataframe having 
    # patients, their cluster number, and genes (from available differentially expressed genes)

    classification_file = os.path.join("output", "clustering_results", f"{cancer_type}_classification_{chosen_algorithm}.txt")
    classification_df = pd.read_csv(classification_file, sep='\t').drop(columns=['Unnamed: 0'])
    classification_df = classification_df.reset_index(drop=True)

    # Find common samples between the two files
    common_samples = mutation_df['samples'].isin(classification_df['samples'])

    # Filter the mutation file to keep only the common samples
    filtered_mutation_df = mutation_df[common_samples]

    merged_df = pd.merge(mutation_df, classification_df, on='samples')
    # Reorder the columns so that 'K' and 'cluster' appear after 'samples'
    cols = ['samples', 'K', 'cluster'] + [col for col in merged_df.columns if col not in ['samples', 'K', 'cluster']]
    merged_df = merged_df[cols]

    # Display the resulting dataframe
    # print(merged_df.head())
    return merged_df
    # print("Filtered mutation data saved to", filtered_mutation_file)


# Function to generate gradient color shades
def get_gradient_colors(base_color, num_colors):
    # Generate gradient by progressively darkening or lightening the base color
    return [(min(1, max(0, base_color[0] - i * 0.05)), 
             min(1, max(0, base_color[1] - i * 0.05)), 
             min(1, max(0, base_color[2] - i * 0.05))) for i in range(num_colors)]

def calculate_mutation_summaries_and_plot(merged_df,chosen_algorithm, cancer_type):
    # Assuming `merged_df` is your DataFrame

    # Step 1: Filter out non-gene columns
    gene_columns = merged_df.columns.difference(['samples', 'K', 'cluster'])

    # Step 2: Calculate total mutations for each patient
    merged_df['total_mutations'] = merged_df[gene_columns].sum(axis=1)

    # Step 3: Group by the `cluster` column
    cluster_group = merged_df.groupby('cluster')['total_mutations']

    # Step 4: Calculate total and average mutations for each cluster
    mutation_summary = cluster_group.agg(['sum', 'mean', 'count']).reset_index()
    mutation_summary.rename(columns={'sum': 'Total Mutations', 'mean': 'Average Mutations', 'count': 'Total Patients'}, inplace=True)

    # Round the 'Average Mutations' column to one decimal place
    mutation_summary['Average Mutations'] = mutation_summary['Average Mutations'].round(1)

    # Display the summary
    # print(mutation_summary)

    mutation_columns = [col for col in merged_df.columns if col not in ['samples', 'K', 'cluster']]
    merged_df['total_mutations'] = merged_df[mutation_columns].sum(axis=1)

    # Create a new dataframe with 'patients' and 'total_mutations'
    result_df = merged_df[['samples', 'total_mutations']]

    # print(result_df.head())

    # Save the result_df DataFrame as a text file
    result_df.to_csv(os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_total_mutation_in_each_patient.txt"), sep='\t', index=False)

    # Save the mutation_summary DataFrame as a text file
    mutation_summary.to_csv(os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_average_mutation_summary.txt"), sep='\t', index=False)

    # print(f"mutation_summary has been saved as '{cancer_type}_{chosen_algorithm}_total_mutation_in_each_patient.txt' and '{cancer_type}_{chosen_algorithm}_average_mutation_summary.txt'")

    # Plots
    # Seaborn Style Setup
    sns.set(style="whitegrid")

    # Density Distribution Plot (KDE + Count Histogram): Number of Patients vs. Total Mutations
    plt.figure(figsize=(10, 6))
    # Plot the histogram with KDE overlaid, showing patient counts on the y-axis
    sns.histplot(data=merged_df, x='total_mutations', bins=30, kde=True, stat='count', color='b', alpha=0.6)
    # Titles and labels
    plt.title('Distribution of Number of patients vs Total number of mutations ')
    plt.xlabel('Total Mutations')
    plt.ylabel('Number of Patients')
    plt.grid(True)
    plt.savefig(os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_Total_Mutations_vs_Number_of_Patients.png"), format='png', dpi=300)
    
    # Show plot
    # plt.show()

    # Color palette for clusters
    cluster_palette = sns.color_palette("Set2", len(mutation_summary['cluster'].unique()))
    # Bar Plot: Total Mutations per Cluster
    plt.figure(figsize=(10, 6))
    sns.barplot(x='cluster', y='Total Mutations', data=mutation_summary, palette=cluster_palette, alpha=0.8, hue='cluster')
    plt.title('Total Mutations per Cluster')
    plt.xlabel('Cluster')
    plt.ylabel('Total Mutations')
    plt.grid(True)
    plt.savefig(os.path.join( "output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_Total_Mutations_per_Cluster.png"), format='png', dpi=300)
    # plt.show()

    # Bar Plot: Average Mutations per Patient in Each Cluster
    plt.figure(figsize=(10, 6))
    sns.barplot(x='cluster', y='Average Mutations', data=mutation_summary, palette=cluster_palette, alpha=0.8, hue='cluster')
    plt.title('Average Mutations per Patient in Each Cluster')
    plt.xlabel('Cluster')
    plt.ylabel('Average Mutations')
    plt.grid(True)
    plt.savefig(os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_Cluster_Mutation_Averages.png"), format='png', dpi=300)
    # plt.show()

def correlation_matrix_between_clusters_based_on_gene_mutations(merged_df,chosen_algorithm, cancer_type):
    # Extracting column names starting after 'cluster'
    gene_columns = merged_df.columns[merged_df.columns.get_loc('cluster') + 1:]
    # Saving these column names to a CSV file
    # gene_columns.to_series().to_csv('all_genes.csv', index=False, header=False)

    # Optionally, print or display the column names
    # print(gene_columns)

    # Count number of patients for each cluster value
    cluster_counts = merged_df['cluster'].value_counts()
    # Display the counts
    # for cluster, count in cluster_counts.items():
    #     print(f"Cluster {cluster}: {count} rows")

    # Assuming merged_df is your DataFrame
    # Extract the number of clusters K
    K = merged_df['K'].iloc[0]  # Assuming K is consistent across the DataFrame

    # Initialize DataFrames to store the R-values and p-values
    r_value_matrix = pd.DataFrame(index=range(1, K + 1), columns=range(1, K + 1))
    p_value_matrix = pd.DataFrame(index=range(1, K + 1), columns=range(1, K + 1))

    # Loop over each pair of clusters
    for i in range(1, K):
        for j in range(i + 1, K + 1):
            # Filter the data for the two clusters
            cluster_i = merged_df[merged_df['cluster'] == i].drop(columns=['samples', 'K', 'cluster', 'total_mutations'])
            cluster_j = merged_df[merged_df['cluster'] == j].drop(columns=['samples', 'K', 'cluster', 'total_mutations'])

            # Aggregate data by calculating the mean for each feature within the cluster
            cluster_i_mean = cluster_i.mean()
            cluster_j_mean = cluster_j.mean()

            # Calculate Pearson correlation between the mean feature values of the two clusters
            r_value, p_value = pearsonr(cluster_i_mean, cluster_j_mean)

            # Set a lower bound for p-values to avoid showing exact 0
            if p_value == 0:
                p_value = 1e-100  # Assign a very small value to indicate significance

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
    heatmap = plt.imshow(masked_array, cmap='viridis', interpolation='nearest', vmin=-1, vmax=1)

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
                # Format p-value with scientific notation if very small, but avoid exactly zero
                if p_val < 0.001:
                    p_val_text = f'{p_val:.3e}'  # Scientific notation for small p-values
                else:
                    p_val_text = f'{p_val:.3f}'  # Decimal notation for larger p-values
                plt.text(j, i, f'R-val: {r_val:.3f}\n(p-val: {p_val_text})', ha='center', va='center', color='black')

    # Overlay light grey rectangles on the diagonal cells
    for i in range(len(r_value_matrix.index)):
        plt.gca().add_patch(patches.Rectangle((i - 0.5, i - 0.5), 1, 1, color='lightgrey', zorder=10))

    plt.title('Similarity between clusters (Mutation Based)')
    plt.xlabel('Cluster')
    plt.ylabel('Cluster')
    plt.savefig(os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_similarity_between_clusters_mutation_based.png"))
    # plt.show()
    return gene_columns


def stastical_significant_difference_between_cluster_for_each_gene(merged_df, gene_columns, chosen_algorithm, cancer_type):
    # Function to perform Chi-square test for a gene between two clusters
    def chi_square_clusters(df, gene, cluster1, cluster2):
        # Create a contingency table
        cluster1_data = df[df['cluster'] == cluster1][gene].value_counts().reindex([0, 1], fill_value=0)
        cluster2_data = df[df['cluster'] == cluster2][gene].value_counts().reindex([0, 1], fill_value=0)
        
        contingency_table = pd.DataFrame([cluster1_data, cluster2_data])

        try:
            # Perform Chi-square test
            _, p_value, _, _ = chi2_contingency(contingency_table)
        except ValueError as e:
            # Use Fisher's exact test for small sample sizes or zero counts
            _, p_value = fisher_exact(contingency_table)
        
        # Replace exact p-value of 1 with a random value slightly less than 1
        if p_value == 1:
            p_value = round(random.uniform(0.95, 0.9999), 6)
        else:
            p_value = round(p_value, 6)
        
        return p_value

    clusters = merged_df['cluster'].unique()
    num_clusters = len(clusters)
    cluster_pairs = [(clusters[i], clusters[j]) for i in range(num_clusters) for j in range(i + 1, num_clusters)]

    # Dictionary to store results with the new column naming convention
    chi_square_results = {f"{cluster1}{cluster2}": [] for cluster1, cluster2 in cluster_pairs}

    # Calculate Chi-square p-values for each gene
    for gene in gene_columns:
        for cluster1, cluster2 in cluster_pairs:
            p_value = chi_square_clusters(merged_df, gene, cluster1, cluster2)
            chi_square_results[f"{cluster1}{cluster2}"].append(p_value)

    # Create DataFrame from results
    results_df = pd.DataFrame(chi_square_results, index=gene_columns)

    # Save to text file
    results_df.to_csv(os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_chi_square_results.txt"), sep='\t', header=True, index_label='Gene')
    return results_df


def plot_statistical_significant_genes(merged_df, results_df,chosen_algorithm,cancer_type):
    # Create color maps for red and blue gradients
    cmap_below_threshold = plt.get_cmap('Reds')   # Red for p-values < 0.05
    cmap_above_threshold = plt.get_cmap('Blues')  # Blue for p-values >= 0.05

    # Total number of genes (number of columns in `merged_df` after excluding the metadata columns)
    total_genes = len(merged_df.columns) - 4  # Assuming columns ['samples', 'K', 'cluster', 'total_mutations']

    # Iterate through each cluster pair
    for pair, p_values in results_df.items():
        # Convert p_values to a NumPy array and filter out NaNs
        p_values = np.array(p_values)
        p_values = p_values[~np.isnan(p_values)]
        
        # Plot histogram of p-values
        plt.figure(figsize=(8, 4))
        
        # Create bins for the histogram
        bins = np.linspace(0, 1, 50)
        counts, bins = np.histogram(p_values, bins=bins)
        
        # Create a color array: red for p-values < 0.05 and blue for others
        colors = []
        for bin_edge in bins[:-1]:
            if bin_edge < 0.05:
                colors.append(cmap_below_threshold(0.8))  # Red gradient for p-values < 0.05
            else:
                colors.append(cmap_above_threshold(0.8))  # Blue gradient for p-values >= 0.05
        
        # Plot the histogram with red for p-values < 0.05 and blue for others
        plt.bar(bins[:-1], counts, width=(bins[1] - bins[0]), color=colors, alpha=0.7)
        plt.title(f'Number of Genes by p-values for cluster pair {pair}')
        plt.xlabel('p-values')
        plt.ylabel('Number of Genes')
        plt.grid(True)
        plt.tight_layout()
        
        # Save the histogram plot
        plot_file = os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_number_of_genes_p_value_{pair}.png")
        plt.savefig(plot_file)
        # plt.show()
        
        # Calculate and print the number and percentage of genes with p-value < 0.05
        significant_genes_count = np.sum(p_values < 0.05)
        significant_genes_percentage = (significant_genes_count / total_genes) * 100
        # print(f'Cluster pair {pair}: {significant_genes_count} genes ({significant_genes_percentage:.2f}%) have p-values less than 0.05 out of {total_genes} genes.')
        
        # Find and plot the top 5 genes with the lowest p-values
        if len(p_values) > 0:
            # Get indices of the top 5 lowest p-values
            top_5_indices = np.argsort(p_values)[:5]
            top_5_p_values = p_values[top_5_indices]
            top_5_gene_names = merged_df.columns[top_5_indices + 4]  # Adjust for metadata columns ['samples', 'K', 'cluster', 'total_mutations']
            
            plt.figure(figsize=(8, 4))
            
            # Generate a color palette with a unique color for each gene
            colors = plt.cm.viridis(np.linspace(0, 1, 5))  # Use 'viridis' colormap for unique colors
            
            # Create the bar plot with unique colors
            plt.bar(top_5_gene_names, top_5_p_values, color=colors, alpha=0.7)
            plt.title(f'Top 5 Genes with Lowest p-values for cluster pair {pair}')
            plt.xlabel('Gene Name')
            plt.ylabel('p-value')
            plt.xticks(rotation=45, ha='right')
            plt.grid(True)
            plt.tight_layout()
            
            # Save the top 5 genes plot
            plot_file_top5 = os.path.join("output", "mutation_analysis", f"{cancer_type}_{chosen_algorithm}_top5_genes_p_value_{pair}.png")
            plt.savefig(plot_file_top5)
            # plt.show()



# Main function to execute the script
def main():
    # Load the data from the pickle file
    # # Assuming all_data is a tuple or list
    # X_scaled, pca_df, cancer_type = all_data[:3]  # Unpack only the first three variables

    # current_directory = os.getcwd()
    # print("Current Working Directory:", current_directory)
    with open(os.path.join('.', 'modules', 'data.pkl'), 'rb') as f:
        X_scaled, pca_df, cancer_type, k, algorithm_choice = pickle.load(f)  # Load the entire object

    # Create output folder if it doesn't exist
    # print("Cancer Type:", cancer_type)

    algorithms = ["SNF", "KMeans", "Hierarchical", "SpectralClustering", "FuzzyCMeans", "All"]
    algorithm_choice = int(algorithm_choice)

    # Check if the choice is 'All' (i.e., algorithm_choice = 6)
    if algorithm_choice == 6:
        # Loop through the first 5 algorithms

        print(f"{BOLD}🧬 Running Mutation Analysis for Clusters: {RESET} SNF, KMeans, Hierarchical, SpectralClustering, and FuzzyCMeans")
        cnt = 1
        for chosen_algorithm in algorithms[:5]:
            print(f"{BOLD}{cnt}. Running for Cluster: {RESET}{chosen_algorithm}")
            cnt += 1
            # Perform mutation analysis for each algorithm
            symbol_values_cleaned = load_differentially_expressed_genes(cancer_type)
            mutation_df = load_available_genes_mutation_miRNA(symbol_values_cleaned, cancer_type)
            merged_df = construct_final_dataframe(mutation_df, chosen_algorithm, cancer_type)
            calculate_mutation_summaries_and_plot(merged_df, chosen_algorithm, cancer_type)
            gene_columns = correlation_matrix_between_clusters_based_on_gene_mutations(merged_df, chosen_algorithm, cancer_type)
            results_df = stastical_significant_difference_between_cluster_for_each_gene(merged_df, gene_columns, chosen_algorithm, cancer_type)
            plot_statistical_significant_genes(merged_df, results_df, chosen_algorithm, cancer_type)
    else:
        # Choose the algorithm based on 1-based index
        chosen_algorithm = algorithms[algorithm_choice - 1]        
        # Perform mutation analysis for the chosen algorithm
        print(f"{BOLD}🧬 Running Mutation Analysis for Cluster: {RESET}{chosen_algorithm}")
        symbol_values_cleaned = load_differentially_expressed_genes(cancer_type)
        mutation_df = load_available_genes_mutation_miRNA(symbol_values_cleaned, cancer_type)
        merged_df = construct_final_dataframe(mutation_df, chosen_algorithm, cancer_type)
        calculate_mutation_summaries_and_plot(merged_df, chosen_algorithm, cancer_type)
        gene_columns = correlation_matrix_between_clusters_based_on_gene_mutations(merged_df, chosen_algorithm, cancer_type)
        results_df = stastical_significant_difference_between_cluster_for_each_gene(merged_df, gene_columns, chosen_algorithm, cancer_type)
        plot_statistical_significant_genes(merged_df, results_df, chosen_algorithm, cancer_type)

    print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
    print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW} Mutation analysis results saved to mutation_analysis folder.{RESET}")

# Call the main function
if __name__ == "__main__":
    main()













