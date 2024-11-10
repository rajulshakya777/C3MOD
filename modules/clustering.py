import pickle
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, ClusterMixin
import skfuzzy as fuzz
import warnings
import rpy2.robjects as robjects

# Color and style codes
YELLOW = "\033[1;33m"
CYAN = "\033[1;36m"
BLUE = "\033[1;34m"
BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[1;32m"
CHECK_EMOJI = f"{GREEN}✔{RESET}"

warnings.filterwarnings("ignore")
output_folder = os.path.join(os.path.dirname(__file__), '..', 'output', 'clustering_results')
os.makedirs(output_folder, exist_ok=True)

# flag = False
# Custom wrapper for Fuzzy C-means for compatibility
class FuzzyCMeansWrapper(BaseEstimator, ClusterMixin):
    def __init__(self, n_clusters):
        self.n_clusters = n_clusters

    def fit(self, X, y=None):
        self.cntr, u, _, _, _, _, _ = fuzz.cluster.cmeans(X.T, self.n_clusters, 2, error=0.005, maxiter=1000)
        self.labels_ = np.argmax(u, axis=0)
        return self

    def predict(self, X):
        return self.labels_


# Function to get valid integer input from user
def get_integer_input(prompt, min_value, max_value):
    while True:
        try:
            value = int(input(prompt))
            if min_value <= value <= max_value:
                return value
            else:
                print(f"{YELLOW}⚠️ Please enter a number between {min_value} and {max_value}.{RESET}")
        except ValueError:
            print(f"{YELLOW}⚠️ Invalid input. Please enter a valid integer.{RESET}")

# Function to get float input from user
def get_float_input(prompt, min_value, max_value):
    while True:
        try:
            value = float(input(prompt))
            if min_value <= value <= max_value:
                return value
            else:
                print(f"{YELLOW}⚠️ Please enter a number between {min_value} and {max_value}.{RESET}")
        except ValueError:
            print(f"{YELLOW}⚠️ Invalid input. Please enter a valid number.{RESET}")



def run_snf(n_clusters, cancer_type):

    selected_cancer = cancer_type
    K = n_clusters
    n_neighbour = 20
    alpha = 0.5
    T = 15

    # Step 3: Change parameters
    print("Do you want to change any of the following parameters for SNF clustering?")
    print("1. Number of neighbors for constructing affinity matrix (Default is 20)")
    print("2. Hyperparameter alpha for constructing affinity matrix (Default is 0.5)")
    print("3. Number of iterations for SNF (Default is 15)")
    print("4. No (Keep the default values)")

    param_choice = get_integer_input("Select an option (1-4): ", 1, 4)

    if param_choice == 1:
        n_neighbour = get_integer_input("Enter the value of Number of neighbors for constructing affinity matrix (5-25): ", 5, 25)
    elif param_choice == 2:
        alpha = get_float_input("Enter the value of Hyperparameter alpha for constructing affinity matrix (0.1-0.9): ", 0.1, 0.9)
    elif param_choice == 3:
        T = get_integer_input("Enter the value of Number of iterations for SNF (5-20): ", 5, 20)

    # Print the final parameters
    # if(flag==True):
    print(f"{BOLD}🔍 Selected algorithm: {RESET} SNF")
    # else:

    print(" SNF Runnning with parameters:")
    print(f" {YELLOW}Cancer Dataset: {selected_cancer}{RESET}")
    print(f" {YELLOW}Number of clusters (K): {K}{RESET}")
    print(f" {YELLOW}Number of neighbors: {n_neighbour}{RESET}")
    print(f" {YELLOW}Hyperparameter alpha: {alpha}{RESET}")
    print(f" {YELLOW}Number of iterations for SNF: {T}{RESET}")

    # print(f"Cancer Dataset: {selected_cancer}")
    # print(f"Number of clusters (K): {K}")
    # print(f"Number of neighbors: {n_neighbour}")
    # print(f"Hyperparameter alpha: {alpha}")
    # print(f"Number of iterations for SNF: {T}")

    # Step 4: Running the R script with parameters
    r_code = f"""
    # library(SNFtool)
    # library(dplyr)
    # library(stringr)
    # library(ggplot2)

    invisible(capture.output({{
      library(SNFtool)
      library(dplyr)
      library(stringr)
      library(ggplot2)
    }}, type = "message"))
    options(warn = -1)

    # suppressMessages(library(dplyr))
    # suppressMessages(library(stringr))
    # suppressMessages(library(ggplot2))

    # Set memory limit based on OS
    if (.Platform$OS.type == "windows") {{
      memory.limit(size = 50 * 1024)  # Set the memory limit to 50GB for Windows
    }} else {{
      Sys.setenv(R_MAX_VSIZE = "50GB")  # Set a memory limit 50GB for macOS/Linux
    }}

    # Cancer Input
    cancers <- list("{selected_cancer}")

    # Maximum number of clusters to predict
    maxK = {K}

    # Number of neighbors for constructing affinity matrix
    n_neighbour = {n_neighbour}

    # Hyperparameter for constructing affinity matrix
    alpha = {alpha}

    # Number of iterations for SNF
    T = {T}

    # Define standard normalization function
    standardNormalization = function(x) {{
      x = as.matrix(x)
      mean = apply(x, 2, mean)
      sd = apply(x, 2, sd)
      sd[sd == 0] = 1
      xNorm = t((t(x) - mean) / sd)
      return(xNorm)
    }}

    # Loop through each cancer type
    for (project in cancers) {{
      tumor_type = project

      ##### Reading Omics Data ####
      ## miRNA data ####
      file_miRNA <- file.path("data", "input_data", "TCGA_data", tumor_type, paste(tumor_type, "miRNA_GANfeature1.txt", sep = "_"))
      miRNA <- read.csv(file_miRNA, header = TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE, check.names = FALSE)
      data_miRNA <- data.frame(t(as.matrix(miRNA)))

      ## mRNA data #####
      file_mRNA <- file.path("data", "input_data", "TCGA_data", tumor_type, paste(tumor_type, "mRNA_GANfeature1.txt", sep = "_"))
      mRNA <- read.csv(file_mRNA, header = TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE, check.names = FALSE)
      data_mRNA <- data.frame(t(as.matrix(mRNA)))

      ## methyl data ####
      file_methyl <- file.path("data", "input_data", "TCGA_data", tumor_type, paste(tumor_type, "methyl_GANfeature1.txt", sep = "_"))
      methyl <- read.csv(file_methyl, header = TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE, check.names = FALSE)
      data_methyl <- data.frame(t(as.matrix(methyl)))

      print("--Omics data reading completed--")

      ##### Run SNF ####
      # Normalize the data
      data_miRNA_Norm = lapply(data_miRNA, standardNormalization)
      data_mRNA_Norm = lapply(data_mRNA, standardNormalization)
      data_methyl_Norm = lapply(data_methyl, standardNormalization)
      print("--Data normalization completed--")

      ## Calculate the distances for each dataset
      dist_miRNA = lapply(as.matrix(data_miRNA_Norm), function(x) dist2(x, x))
      dist_mRNA = lapply(as.matrix(data_mRNA_Norm), function(x) dist2(x, x))
      dist_methyl = lapply(as.matrix(data_methyl_Norm), function(x) dist2(x, x))
      print("--Distance matrix calculation completed--")

      ## Construct the similarity graphs
      affinityL_miRNA = lapply(dist_miRNA, function(x) affinityMatrix(x, n_neighbour, alpha))
      affinityL_mRNA = lapply(dist_mRNA, function(x) affinityMatrix(x, n_neighbour, alpha))
      affinityL_methyl = lapply(dist_methyl, function(x) affinityMatrix(x, n_neighbour, alpha))
      print("--Similarity graph network constructed--")

      ## Construct the fused network
      W_miRNA = SNF(affinityL_miRNA, n_neighbour, T)
      W_mRNA = SNF(affinityL_mRNA, n_neighbour, T)
      W_methyl = SNF(affinityL_methyl, n_neighbour, T)
      W_integrated = SNF(list(W_miRNA, W_mRNA, W_methyl), n_neighbour, T)
      # print("--Network fusion completed--")

      ## Perform clustering on the fused network
      clus = spectralClustering(W_integrated, maxK)
      df = data.frame(samples = rownames(data_miRNA), K = rep(maxK, length(rownames(data_miRNA))), cluster = clus)

      # Save clustering results
      outFile_clus <- file.path("output", "clustering_results", paste(tumor_type, "classification_SNF.txt", sep = "_"))
      write.table(df, file = outFile_clus, sep = "\t", quote = FALSE, col.names = NA)
      print("---- SNF clusters saved ----")

      ## PCA and Plot
      clus <- spectralClustering(W_integrated, maxK)
      df <- data.frame(samples = rownames(data_miRNA), cluster = as.factor(clus))
      pca_result <- prcomp(W_integrated, center = TRUE, scale. = TRUE)
      pca_df <- as.data.frame(pca_result$x)
      pca_df$cluster <- df$cluster

      ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
          geom_point(size = 3) +
          labs(title = paste("PCA of Patients with SNF Clustering"),
              x = "Principal Component 1",
              y = "Principal Component 2") +
          theme_minimal() +
          theme(legend.title = element_blank(), 
                legend.position = "right",
                panel.background = element_rect(fill = "white", color = "white"),
                plot.background = element_rect(fill = "white", color = "white"),
                panel.grid.major = element_line(color = "grey80"),
                panel.grid.minor = element_line(color = "grey90"),
                plot.title = element_text(color = "black", hjust = 0.5),
                axis.title = element_text(color = "black"),
                axis.text = element_text(color = "black"),
                panel.border = element_rect(color = "black", fill = NA, size = 1))

      ggsave(filename = file.path("output", "clustering_results", paste(tumor_type, "SNF_pca.png", sep = "_")), plot = last_plot(), width = 10, height = 8)
      # print("---- PCA plot generated ----")
      # print("Analysis summary:")
      # print("SNF Clustering results saved in clustering_results folder")
    }}
    """
    # Run the R script using rpy2
    robjects.r(r_code)
    print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
    print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW}All Clustering results saved to clustering_results folder.{RESET}")

# Function to get user input for algorithm and number of clusters
def get_user_input():
    print("Select the algorithm to run:")
    print("1. SNF")
    print("2. KMeans")
    print("3. Hierarchical Clustering")
    print("4. Spectral Clustering")
    print("5. Fuzzy C-Means")
    print("6. Run All")

    while True:
        try:
            algorithm_choice = int(input("Enter the number corresponding to the algorithm (1-6): "))
            if 1 <= algorithm_choice <= 6:
                break
            else:
                print(f"{YELLOW}⚠️ Invalid choice. Please select a number between 1 and 6.{RESET}")
        except ValueError:
            print(f"{YELLOW}⚠️ Invalid input. Please enter a valid number.{RESET}")


    while True:
        try:
            k = int(input("Enter the number of clusters (k) between 2 and 7: "))
            if 2 <= k <= 7:
                break
            else:
                print(f"{YELLOW}⚠️ Invalid choice. Please select a number between 2 and 7.{RESET}")
        except ValueError:
            print(f"{YELLOW}⚠️ Invalid input. Please enter a valid number.{RESET}")
    

    return algorithm_choice, k


# Function to define the clustering methods based on user choice
def define_clustering_methods(algorithm_choice, k, cancer_type):

    if algorithm_choice == 1:
        run_snf(k,cancer_type)
        return "SNF"
    if algorithm_choice == 2:
        return {'KMeans': KMeans(n_clusters=k, random_state=42)}
    elif algorithm_choice == 3:
        return {'Hierarchical': AgglomerativeClustering(n_clusters=k)}
    elif algorithm_choice == 4:
        return {'SpectralClustering': SpectralClustering(n_clusters=k, assign_labels="discretize", random_state=42)}
    elif algorithm_choice == 5:
        return {'FuzzyCMeans': FuzzyCMeansWrapper(n_clusters=k)}
    elif algorithm_choice == 6 :
        # print(f"{BOLD}🔍 Selected algorithm: {RESET}{method_name}")
        # global flag
        # flag = True
        print(f"{BOLD}🔍 Selected Choice: {RESET}Run All")
        run_snf(k,cancer_type)
        return {
            'KMeans': KMeans(n_clusters=k, random_state=42),
            'Hierarchical': AgglomerativeClustering(n_clusters=k),
            'SpectralClustering': SpectralClustering(n_clusters=k, assign_labels="discretize", random_state=42),
            'FuzzyCMeans': FuzzyCMeansWrapper(n_clusters=k)
        }
    else:
        print("Invalid choice")
        exit()


# Function to perform clustering and save the results
def perform_clustering_and_save_results(clustering_methods, X_scaled, pca_df, output_folder):
    # Ensure reproducibility of colors
    np.random.seed(42)

    # Define a list of colors
    colors = ['r', 'g', 'b', 'y', 'c', 'm', 'k']

    # cnt = 2;
    for method_name, model in clustering_methods.items():
        # print(method_name)
        # print(f"Algorithm: {method_name}")
        # Fit the clustering algorithm to the original data
        model.fit(X_scaled)
        
        # Predict cluster labels, handle cases where model has no predict method
        if hasattr(model, 'predict'):
            cluster_labels = model.predict(X_scaled)
        else:
            cluster_labels = model.labels_

        # Add cluster labels to the PCA DataFrame
        pca_df['cluster'] = cluster_labels + 1

        # Rename 'PatId' column to 'samples'
        pca_df.rename(columns={'PatId': 'samples'}, inplace=True)

        # Reset index starting from 1
        pca_df.reset_index(drop=True, inplace=True)
        pca_df.index += 1  # Start index from 1

        # Add a column for the number of clusters used
        pca_df['K'] = model.n_clusters

        # Save classification file for each algorithm as a text file
        classification_file = os.path.join(output_folder, f'LUAD_classification_{method_name}.txt')
        pca_df[['samples', 'K', 'cluster']].to_csv(classification_file, sep='\t', index=True)

        # Sort the unique cluster labels to ensure consistent labeling
        unique_clusters = sorted(pca_df['cluster'].unique())

        # Create a new figure for each clustering method
        plt.figure(figsize=(10, 8))

        # print(f"Algorithm: {method_name}")
        print(f"{BOLD}🔍 Selected algorithm: {RESET}{method_name}")

        # Plot the PCA result
        for cluster_label in unique_clusters:
            cluster_data = pca_df[pca_df['cluster'] == cluster_label]
            plt.scatter(cluster_data['PC1'], cluster_data['PC2'], 
                        c=colors[cluster_label % len(colors)], 
                        label=f'cluster {cluster_label}')

        plt.title(f'PCA of Patients with {method_name}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend()
        plt.grid(True)

        # Save the plot to a file
        plot_file = os.path.join(output_folder, f'pca_{method_name}.png')
        plt.savefig(plot_file)
        print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
        print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW}All Clustering results saved to clustering_results folder.{RESET}")


# Main function to execute the script
def main():
    # Load the data from the pickle file
    with open(os.path.join('.', 'modules', 'data.pkl'), 'rb') as f:
        all_data = pickle.load(f)  # Load the entire object

    # Assuming all_data is a tuple or list
    X_scaled, pca_df, cancer_type = all_data[:3]  # Unpack only the first three variables

    # Create output folder if it doesn't exist
    # output_folder = './../output/clustering_results'
    # if not os.path.exists(output_folder):
    #     os.makedirs(output_folder)

    # Get user input
    algorithm_choice, k = get_user_input()
    print(f"{BOLD}🔍 Selected Number of Clusters(K): {RESET}{k}")
    # print(f"{BOLD}🔍 Selected algorithm: {RESET}SNF")
    # print(f"{BOLD}🔍 Selected algorithm: {RESET}SNF")
    
    
    # Add the new variables to the existing data
    data_to_save = (X_scaled, pca_df, cancer_type, k, algorithm_choice)

    # Save the updated data back to the pickle file
    with open(os.path.join('.', 'modules', 'data.pkl'), 'wb') as f:
        pickle.dump(data_to_save, f)

    # Define clustering methods
    clustering_methods = define_clustering_methods(algorithm_choice, k, cancer_type)

    # If SNF was chosen, exit after running SNF
    if clustering_methods == "SNF":
        return

    # Perform clustering and save results
    perform_clustering_and_save_results(clustering_methods, X_scaled, pca_df, output_folder)


# Call the main function
if __name__ == "__main__":
    main()
