import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import skfuzzy as fuzz
import pickle
import warnings
from rpy2 import robjects
from rpy2.robjects import pandas2ri
import os

warnings.filterwarnings("ignore")
# Activate the automatic conversion of pandas DataFrames to R DataFrames
pandas2ri.activate()

# Color and style codes
YELLOW = "\033[1;33m"
CYAN = "\033[1;36m"
BLUE = "\033[1;34m"
BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[1;32m"
CHECK_EMOJI = f"{GREEN}✔{RESET}"

# Ensure the output directory exists
# output_directory = './../output/recommended_K'

# Load other libraries (like SNFtool) in R
# robjects.r('''
#     library(cluster)
#     standardNormalization = function(x) {
#         x = as.matrix(x)
#         mean = apply(x, 2, mean)
#         sd = apply(x, 2, sd)
#         sd[sd == 0] = 1
#         xNorm = t((t(x) - mean) / sd)
#         return(xNorm)
#     }

#     plot_silhouette <- function(dist_matrix, maxK = 10) {
#         sil_scores = numeric(maxK - 1)

#         for (k in 2:maxK) {
#             clus = kmeans(dist_matrix, centers = k, nstart = 25)
#             sil = silhouette(clus$cluster, dist_matrix)
#             sil_scores[k - 1] = mean(sil[, 3])
#         }

#         best_K = which.max(sil_scores) + 1
#         print(paste("Recommended K for SNF algorithm based on Silhouette Score:", best_K))

#         plot(2:maxK, sil_scores, type = "b", pch = 19, xlab = "Number of Clusters (K)", ylab = "Average Silhouette Score",
#              main = "Silhouette Score for Different K")

#         return(best_K)
#     }
# ''')

# Define a range of K values
k_values = range(2, 7)

# Initialize lists to store results
results = {
    'KMeans': {'K': [], 'Silhouette Score': [], 'WCSS': []},
    'Hierarchical': {'K': [], 'Silhouette Score': [], 'Total Within-Cluster Sum of Squares': []},
    'SpectralClustering': {'K': [], 'Silhouette Score': [], 'Gap Statistic': []},
    'FuzzyCMeans': {'K': [], 'Fuzzy Partition Coefficient': [], 'Fuzzy Silhouette Index': []}
}

def compute_gap_statistic(X, k, n_refs=10):
    """
    Compute the Gap Statistic for a given number of clusters k.
    
    Parameters:
    - X: The data array
    - k: Number of clusters
    - n_refs: Number of reference datasets to generate
    
    Returns:
    - gap_stat: Gap Statistic value
    """
    spectral = SpectralClustering(n_clusters=k, assign_labels="discretize", random_state=42)
    labels = spectral.fit_predict(X)

    def compute_wcss(X, labels):
        return np.sum([np.sum((X[labels == i] - np.mean(X[labels == i], axis=0)) ** 2) for i in np.unique(labels)])

    wcss_actual = compute_wcss(X, labels)

    wcss_refs = []
    for _ in range(n_refs):
        X_ref = np.random.uniform(low=np.min(X, axis=0), high=np.max(X, axis=0), size=X.shape)
        labels_ref = spectral.fit_predict(X_ref)
        wcss_ref = compute_wcss(X_ref, labels_ref)
        wcss_refs.append(wcss_ref)

    wcss_refs = np.array(wcss_refs)
    gap_stat = np.mean(np.log(wcss_refs)) - np.log(wcss_actual)

    return gap_stat

def run_kmeans(X_scaled):
    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=42)
        labels = kmeans.fit_predict(X_scaled)
        results['KMeans']['K'].append(k)
        results['KMeans']['Silhouette Score'].append(silhouette_score(X_scaled, labels))
        results['KMeans']['WCSS'].append(kmeans.inertia_)

def run_hierarchical(X_scaled):
    for k in k_values:
        hierarchical = AgglomerativeClustering(n_clusters=k)
        labels = hierarchical.fit_predict(X_scaled)
        results['Hierarchical']['K'].append(k)
        results['Hierarchical']['Silhouette Score'].append(silhouette_score(X_scaled, labels))
        within_cluster_ss = np.sum([np.sum((X_scaled[labels == i] - np.mean(X_scaled[labels == i], axis=0)) ** 2) for i in range(k)])
        results['Hierarchical']['Total Within-Cluster Sum of Squares'].append(within_cluster_ss)

def run_spectral(X_scaled):
    for k in k_values:
        gap_stat = compute_gap_statistic(X_scaled, k)
        spectral = SpectralClustering(n_clusters=k, assign_labels="discretize", random_state=42)
        labels = spectral.fit_predict(X_scaled)
        results['SpectralClustering']['K'].append(k)
        results['SpectralClustering']['Silhouette Score'].append(silhouette_score(X_scaled, labels))
        results['SpectralClustering']['Gap Statistic'].append(gap_stat)

def run_fuzzycmeans(X_scaled):
    for k in k_values:
        cmeans = fuzz.cluster.cmeans(X_scaled.T, c=k, m=2, error=0.005, maxiter=1000, init=None)
        labels = np.argmax(cmeans[1], axis=0)
        results['FuzzyCMeans']['K'].append(k)
        fpc = np.mean(np.sum(cmeans[1]**2, axis=0))
        results['FuzzyCMeans']['Fuzzy Partition Coefficient'].append(fpc)
        fsi = np.mean([silhouette_score(X_scaled, labels)])
        results['FuzzyCMeans']['Fuzzy Silhouette Index'].append(fsi)


def run_snf(cancer_type):

    # Load other libraries (like SNFtool) in R
    robjects.r('''
        library(cluster)
        standardNormalization = function(x) {
            x = as.matrix(x)
            mean = apply(x, 2, mean)
            sd = apply(x, 2, sd)
            sd[sd == 0] = 1
            xNorm = t((t(x) - mean) / sd)
            return(xNorm)
        }

        plot_silhouette <- function(dist_matrix, maxK = 10) {
            sil_scores = numeric(maxK - 1)

            for (k in 2:maxK) {
                clus = kmeans(dist_matrix, centers = k, nstart = 25)
                sil = silhouette(clus$cluster, dist_matrix)
                sil_scores[k - 1] = mean(sil[, 3])
            }

            best_K = which.max(sil_scores) + 1
            print(paste("Recommended K for SNF algorithm based on Silhouette Score:", best_K))

            plot(2:maxK, sil_scores, type = "b", pch = 19, xlab = "Number of Clusters (K)", ylab = "Average Silhouette Score",
                 main = "Silhouette Score for Different K")

            return(best_K)
        }
    ''')
    # Load and transpose your data files
    base_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'input_data', 'TCGA_data')
    try:
        # Load and transpose your data files
        miRNA_data = pd.read_csv(os.path.join(base_path, cancer_type, f"{cancer_type}_miRNA_GANfeature.txt"), sep="\t", header=0, index_col=0).T
        mRNA_data = pd.read_csv(os.path.join(base_path, cancer_type, f"{cancer_type}_mRNA_GANfeature.txt"), sep="\t", header=0, index_col=0).T
        methyl_data = pd.read_csv(os.path.join(base_path, cancer_type, f"{cancer_type}_methyl_GANfeature.txt"), sep="\t", header=0, index_col=0).T

    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # Normalize the data
    try:
        # print("Normalizing miRNA data...")
        miRNA_norm = pandas2ri.py2rpy(miRNA_data)
        robjects.globalenv['data_miRNA_Norm'] = robjects.r['standardNormalization'](miRNA_norm)
        # print("miRNA normalization complete.")

        # print("Normalizing mRNA data...")
        mRNA_norm = pandas2ri.py2rpy(mRNA_data)
        robjects.globalenv['data_mRNA_Norm'] = robjects.r['standardNormalization'](mRNA_norm)
        # print("mRNA normalization complete.")

        # print("Normalizing methylation data...")
        methyl_norm = pandas2ri.py2rpy(methyl_data)
        robjects.globalenv['data_methyl_Norm'] = robjects.r['standardNormalization'](methyl_norm)
        # print("Methylation normalization complete.")

    except Exception as e:
        print(f"Error during normalization: {e}")
        return

    # Check if all datasets have the same number of rows
    try:
        miRNA_rows = robjects.r['dim'](robjects.globalenv['data_miRNA_Norm'])[0]
        mRNA_rows = robjects.r['dim'](robjects.globalenv['data_mRNA_Norm'])[0]
        methyl_rows = robjects.r['dim'](robjects.globalenv['data_methyl_Norm'])[0]

        # print(f"Row counts - miRNA: {miRNA_rows}, mRNA: {mRNA_rows}, methylation: {methyl_rows}")

        if miRNA_rows == mRNA_rows == methyl_rows:
            # Combine normalized datasets
            # print("Combining normalized datasets...")
            robjects.r('''
                library(cluster)
                standardNormalization = function(x) {
                    x = as.matrix(x)
                    mean = apply(x, 2, mean)
                    sd = apply(x, 2, sd)
                    sd[sd == 0] = 1
                    xNorm = t((t(x) - mean) / sd)
                    return(xNorm)
                }

                plot_silhouette <- function(dist_matrix, maxK = 10) {
                    sil_scores = numeric(maxK - 1)

                    for (k in 2:maxK) {
                        clus = kmeans(dist_matrix, centers = k, nstart = 25)
                        sil = silhouette(clus$cluster, dist_matrix)
                        sil_scores[k - 1] = mean(sil[, 3])
                    }

                    best_K = which.max(sil_scores) + 1
                    print(paste("Recommended K for SNF algorithm based on Silhouette Score:", best_K))

                    plot(2:maxK, sil_scores, type = "b", pch = 19, xlab = "Number of Clusters (K)", ylab = "Average Silhouette Score",
                         main = "Silhouette Score for Different K")

                    return(best_K)
                }

                data_combined = cbind(data_miRNA_Norm, data_mRNA_Norm, data_methyl_Norm)
                dist_matrix = dist(data_combined)

                # Define function to compute and plot Silhouette Scores
                plot_silhouette <- function(dist_matrix, maxK = 10) {
                    sil_scores = numeric(maxK - 1)
                    
                    for (k in 2:maxK) {
                        clus = kmeans(dist_matrix, centers = k, nstart = 25)
                        sil = silhouette(clus$cluster, dist_matrix)
                        sil_scores[k - 1] = mean(sil[, 3])
                    }
                    
                    # Find the best K based on maximum silhouette score
                    best_K = which.max(sil_scores) + 1
                    
                    # Print the best K
                    # print(paste("Recommended K for SNF algorithm based on Silhouette Score:", best_K))
                    
                    # Set plot background to white
                    par(bg = "white")
                    
                    # Plot the silhouette scores
                    plot(2:maxK, sil_scores, type = "b", pch = 19, xlab = "Number of Clusters (K)", ylab = "Average Silhouette Score",
                         main = "SNF - Silhouette Score for Different K", col = "blue", bg = "white")
                    
                    # print("---getwd----")
                    # print(getwd())

                    if (!dir.exists("./output/recommended_K")) {
                      dir.create("./output/recommended_K", recursive = TRUE)
                    }
                    # dev.copy(png, filename = './../output/recommended_K/SNF_recommended_k.png', bg = "white")
                    dev.copy(png, filename = file.path("output", "recommended_K", "SNF_recommended_k.png"), bg = "white")
                    dev.off()
                    # print("SNF best recommended K plot saved as SNF_recommended_k_plots.png.")
                    print(paste("Recommended K for SNF algorithm based on Silhouette Score:", best_K))
                    return(best_K)
                }

                # Plot Silhouette Method and print best K
                best_K <- plot_silhouette(dist_matrix, maxK = 10)
            ''')

        else:
            print(f"Row mismatch: miRNA ({miRNA_rows}), mRNA ({mRNA_rows}), methylation ({methyl_rows}). Ensure all datasets have the same number of samples.")

    except Exception as e:
        print(f"Error during combination or silhouette plotting: {e}")


def plot_results(df, method_name, metrics, titles):
    plt.figure(figsize=(12, 8))
    colors = ['b', 'g', 'r', 'c']  # Define different colors for each metric
    for idx, metric in enumerate(metrics):
        plt.subplot(2, 1, idx + 1)
        plt.plot(df['K'], df[metric], marker='o', color=colors[idx], label=metric)
        plt.xlabel('Number of Clusters (k)')
        plt.ylabel(metric)
        plt.title(f'{method_name} - {titles[idx]}')
        plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, f"{method_name}_recommended_k_plots.png")) # Save the plot as a PNG file
    # plt.savefig(f'{method_name}_results.png')  # Save the plot as a PNG file
    # plt.show()

def find_best_k(df, metric):
    if metric == 'WCSS':
        return df.loc[df['WCSS'].diff().abs().idxmax(), 'K']
    elif metric == 'Silhouette Score':
        return df.loc[df['Silhouette Score'].idxmax(), 'K']
    elif metric == 'Total Within-Cluster Sum of Squares':
        return df.loc[df['Total Within-Cluster Sum of Squares'].diff().abs().idxmax(), 'K']
    elif metric == 'Gap Statistic':
        return df.loc[df['Gap Statistic'].idxmax(), 'K']
    elif metric == 'Fuzzy Partition Coefficient':
        return df.loc[df['Fuzzy Partition Coefficient'].idxmin(), 'K']
    elif metric == 'Fuzzy Silhouette Index':
        return df.loc[df['Fuzzy Silhouette Index'].idxmax(), 'K']

# Existing code...

def select_algorithm():
    with open(os.path.join('.', 'modules', 'data.pkl'), 'rb') as f:
        all_data = pickle.load(f)  # Load the entire object
    # Assuming all_data is a tuple or list
    X_scaled, pca_df, cancer_type = all_data[:3]

    global output_directory
    output_directory = os.path.join(os.path.dirname(__file__), '..', 'output', 'recommended_K')
    os.makedirs(output_directory, exist_ok=True)

    print("Select the algorithm for which you want to find the recommended K:")
    print("1. SNF")
    print("2. KMeans")
    print("3. Hierarchical Clustering")
    print("4. Spectral Clustering")
    print("5. Fuzzy C-Means")
    print("6. Run All")

    # choice = input("Enter your choice (1-6): ")
    while True:
        try:
            choice = int(input("Enter the number corresponding to your choice (1-6): "))

            if choice == 1:
                print(f"{BOLD}🔍 Selected algorithm: {RESET}SNF")
                print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
                run_snf(cancer_type)
                # Assuming SNF doesn't involve multiple K metrics, you'd print the result as required.

            elif choice == 2:
                print(f"{BOLD}🔍 Selected algorithm: {RESET}KMeans")
                print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
                run_kmeans(X_scaled)
                kmeans_df = pd.DataFrame(results['KMeans'])
                plot_results(kmeans_df, 'KMeans', ['Silhouette Score', 'WCSS'], ['Silhouette Score', 'WCSS'])

                best_k_silhouette = find_best_k(kmeans_df, 'Silhouette Score')
                best_k_wcss = find_best_k(kmeans_df, 'WCSS')
                print(f"  {YELLOW} • Based on Silhouette Score, recommended K is:{RESET} {best_k_silhouette}")
                print(f"  {YELLOW} • Based on WCSS, recommended K is:{RESET} {best_k_wcss}")

            elif choice == 3:
                print(f"{BOLD}🔍 Selected algorithm: {RESET}Hierarchical Clustering")
                print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
                run_hierarchical(X_scaled)
                hierarchical_df = pd.DataFrame(results['Hierarchical'])
                plot_results(hierarchical_df, 'Hierarchical Clustering', ['Silhouette Score', 'Total Within-Cluster Sum of Squares'], ['Silhouette Score', 'Total Within-Cluster Sum of Squares'])

                best_k_silhouette = find_best_k(hierarchical_df, 'Silhouette Score')
                best_k_wcss = find_best_k(hierarchical_df, 'Total Within-Cluster Sum of Squares')

                print(f"  {YELLOW} • Based on Silhouette Score, recommended K is:{RESET} {best_k_silhouette}")
                print(f"  {YELLOW} • Based on WCSS, recommended K is:{RESET} {best_k_wcss}")

            elif choice == 4:
                print(f"{BOLD}🔍 Selected algorithm: {RESET}Spectral Clustering")
                print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
                run_spectral(X_scaled)
                spectral_df = pd.DataFrame(results['SpectralClustering'])
                plot_results(spectral_df, 'Spectral Clustering', ['Silhouette Score', 'Gap Statistic'], ['Silhouette Score', 'Gap Statistic'])

                best_k_silhouette = find_best_k(spectral_df, 'Silhouette Score')
                best_k_gap = find_best_k(spectral_df, 'Gap Statistic')

                print(f"  {YELLOW} • Based on Silhouette Score, recommended K is:{RESET} {best_k_silhouette}")
                print(f"  {YELLOW} • Based on Gap Statistic, recommended K is:{RESET} {best_k_gap}")

            elif choice == 5:
                print(f"{BOLD}🔍 Selected algorithm: {RESET}Fuzzy C-Means")
                print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
                run_fuzzycmeans(X_scaled)
                fuzzy_df = pd.DataFrame(results['FuzzyCMeans'])
                plot_results(fuzzy_df, 'Fuzzy C-Means', ['Fuzzy Partition Coefficient', 'Fuzzy Silhouette Index'], ['Fuzzy Partition Coefficient', 'Fuzzy Silhouette Index'])

                best_k_fpc = find_best_k(fuzzy_df, 'Fuzzy Partition Coefficient')
                best_k_silhouette = find_best_k(fuzzy_df, 'Fuzzy Silhouette Index')

                print(f"  {YELLOW} • Based on Fuzzy Partition Coefficient, recommended K is:{RESET} {best_k_fpc}")
                print(f"  {YELLOW} • Based on Fuzzy Silhouette Index, recommended K is:{RESET} {best_k_silhouette}")

            elif choice == 6:
                # Run all algorithms
                run_kmeans(X_scaled)
                run_hierarchical(X_scaled)
                run_spectral(X_scaled)
                run_fuzzycmeans(X_scaled)

                # Gather results and plot for all algorithms
                kmeans_df = pd.DataFrame(results['KMeans'])
                hierarchical_df = pd.DataFrame(results['Hierarchical'])
                spectral_df = pd.DataFrame(results['SpectralClustering'])
                fuzzy_df = pd.DataFrame(results['FuzzyCMeans'])

                print(f"{BOLD}🔍 Selected algorithm: {RESET}Run All")
                print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
                print("1. KMeans")
                plot_results(kmeans_df, 'KMeans', ['Silhouette Score', 'WCSS'], ['Silhouette Score', 'WCSS'])
                best_k_kmeans_silhouette = find_best_k(kmeans_df, 'Silhouette Score')
                best_k_kmeans_wcss = find_best_k(kmeans_df, 'WCSS')
                print(f"  {YELLOW} • Based on Silhouette Score, recommended K is:{RESET} {best_k_kmeans_silhouette}")
                print(f"  {YELLOW} • Based on WCSS, recommended K is:{RESET} {best_k_kmeans_wcss}")

                print("2. Hierarchical Clustering")
                plot_results(hierarchical_df, 'Hierarchical Clustering', ['Silhouette Score', 'Total Within-Cluster Sum of Squares'], ['Silhouette Score', 'Total Within-Cluster Sum of Squares'])
                best_k_hierarchical_silhouette = find_best_k(hierarchical_df, 'Silhouette Score')
                best_k_hierarchical_wcss = find_best_k(hierarchical_df, 'Total Within-Cluster Sum of Squares')
                print(f"  {YELLOW} • Based on Silhouette Score, recommended K is:{RESET} {best_k_hierarchical_silhouette}")
                print(f"  {YELLOW} • Based on Total Within-Cluster Sum of Squares, recommended K is:{RESET} {best_k_hierarchical_wcss}")        

                print("3. Spectral Clustering")
                plot_results(spectral_df, 'Spectral Clustering', ['Silhouette Score', 'Gap Statistic'], ['Silhouette Score', 'Gap Statistic'])
                best_k_spectral_silhouette = find_best_k(spectral_df, 'Silhouette Score')
                best_k_spectral_gap = find_best_k(spectral_df, 'Gap Statistic')   
                print(f"  {YELLOW} • Based on Silhouette Score, recommended K is:{RESET} {best_k_spectral_silhouette}")
                print(f"  {YELLOW} • Based on Gap Statistic, recommended K is:{RESET} {best_k_spectral_gap}")

                print("4. Fuzzy C-Means")
                plot_results(fuzzy_df, 'Fuzzy C-Means', ['Fuzzy Partition Coefficient', 'Fuzzy Silhouette Index'], ['Fuzzy Partition Coefficient', 'Fuzzy Silhouette Index'])
                best_k_fuzzy_fpc = find_best_k(fuzzy_df, 'Fuzzy Partition Coefficient')
                best_k_fuzzy_silhouette = find_best_k(fuzzy_df, 'Fuzzy Silhouette Index')
                print(f"  {YELLOW} • Based on Fuzzy Partition Coefficient, recommended K is:{RESET} {best_k_fuzzy_fpc}")
                print(f"  {YELLOW} • Based on Fuzzy Silhouette Index, recommended K is:{RESET} {best_k_fuzzy_silhouette}")

                print("5. SNF")
                run_snf(cancer_type)
            else:
                print(f"{YELLOW}⚠️ Invalid choice. Please select a number between 1 and 6.{RESET}")
                continue  # Continues the loop to ask for valid input
            break  # Exit the loop after a valid choice is made

        except ValueError:
            print(f"{YELLOW}⚠️ Invalid input. Please enter a valid number.{RESET}")

    print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW}All WCSS and elbow method plots saved to recommended_K folder.{RESET}")

if __name__ == "__main__":
    select_algorithm()
