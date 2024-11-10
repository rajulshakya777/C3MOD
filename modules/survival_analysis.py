import os
import pickle
import rpy2.robjects as robjects
from rpy2.robjects import r


# Color and style codes
YELLOW = "\033[1;33m"
CYAN = "\033[1;36m"
BLUE = "\033[1;34m"
BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[1;32m"
CHECK_EMOJI = f"{GREEN}✔{RESET}"

# Define the survival analysis function

output_dir = os.path.join(os.path.dirname(__file__), '..', 'output', 'survival_analysis')
os.makedirs(output_dir, exist_ok=True)

def survival_analysis(selected_cancer, K, algorithm_choice):
    # Set the selected cancer, K, and algorithm choice in the R environment
    r(f'cancers <- list("{selected_cancer}")')
    r(f'K <- {K}')
    # r(f'df_methods <- data.frame(methods = c("{algorithm_choice}"))')
    
    # R script for survival analysis
    r(f'''
      invisible(capture.output({{
      library(stringr)
      library(survival)
      library(survminer)
    }}, type = "message"))
    options(warn = -1)

    # Initialize a data frame to store p-values for the specified method
    # methods <- c("{algorithm_choice}")

    algorithms <- c("SNF", "KMeans", "Hierarchical", "SpectralClustering", "FuzzyCMeans")
    algorithm_choice <- {algorithm_choice}  # Change this to test different choices

    # Select the method(s) based on the algorithm_choice
    if (algorithm_choice == 6) {{
      # If choice is 6, select all methods
      selected_methods <- algorithms
    }} else {{
      # Otherwise, select only the chosen method (1-based indexing)
      selected_methods <- algorithms[algorithm_choice]
    }}

    # Define df_methods and methods based on the selected methods
    df_methods <- data.frame(methods = selected_methods)
    methods <- selected_methods

    for (project in cancers) {{
      tumor_type <- substring(project, 6)
      # print(project)
      pval_list <- c()
      selected_cancer <- "{selected_cancer}"
      # Survival data
      str1 <- paste(tumor_type, "survival_UCal.tsv", sep = ".")
      file_survival <- file.path("data", "input_data", "TCGA_data", selected_cancer, paste0(selected_cancer, str1))
      survival <- read.csv(file_survival, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
      dim(survival)
      head(survival, 3)
      replacements <- str_replace_all(rownames(survival), "-", ".")
      rownames(survival) <- replacements
      dim(survival)
      head(survival, 3)
      
      # Process the specified clustering method only
      for (method in methods) {{
        str1 <- paste0(project, "_classification_", method, ".txt")
        file_cl <- file.path("output", "clustering_results", str1)
        if (!file.exists(file_cl)) {{
          print(paste("File not found:", file_cl))
          pval_list <- append(pval_list, NA)
          next
        }}
        
        cluster_data <- read.csv(file_cl, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
        cluster_data <- cluster_data[cluster_data$K == K, ]
        rownames(cluster_data) <- cluster_data$samples
        cluster_data <- cluster_data[, c(which(colnames(cluster_data) == "cluster"), which(colnames(cluster_data) != "cluster"))]
        cluster_data <- cluster_data[1:(length(cluster_data) - 1)]
        dim(cluster_data)
        head(cluster_data, 3)
        
        common_samples_list <- Reduce(intersect, list(as.list(rownames(survival)), as.list(rownames(cluster_data))))
        common_samples_list <- unlist(common_samples_list)
        # print(paste("Length of common samples between survival data and", method, "clustering data:", length(common_samples_list)))
        head(common_samples_list)
        
        survival_reduced <- survival[common_samples_list, c("OS", "OS.time")]
        cluster_data_reduced <- cluster_data[common_samples_list, ]
        
        df_merge <- merge(survival_reduced, cluster_data_reduced, by = "row.names", all = TRUE)
        df_merge <- df_merge[order(row.names(df_merge)), ]
        head(df_merge)
        dim(df_merge)
        
        res.cox <- coxph(Surv(OS.time, OS) ~ cluster, data = df_merge)
        res.cox.sum <- summary(res.cox)$coefficients
        pval <- res.cox.sum[, 5]
        # print(paste(method, "p-value:", pval))
        pval_list <- append(pval_list, pval)
        
        fit <- survfit(Surv(OS.time, OS) ~ cluster, data = df_merge)
        
        # Additional survival plotting and p-value matrix generation c

        # Calculate pairwise comparisons
        pairwise_results <- pairwise_survdiff(Surv(OS.time, OS) ~ cluster, data = df_merge)

        # Extract p-values from pairwise comparisons
        pvals <- pairwise_results$p.value

        # Create a p-value matrix for display
        pval_matrix <- round(pvals, 3)
        pval_text <- apply(pval_matrix, 1, function(x) paste(names(x), x, sep = " = ", collapse = "\n"))

        colors <- rainbow(K)
        # Generate the survival plot
        surv_plot <- ggsurvplot(fit, 
                                pval = FALSE, conf.int = TRUE,
                                risk.table = TRUE,
                                risk.table.col = "strata",
                                linetype = "strata",
                                ggtheme = theme_bw(),
                                palette = colors,
                                xlab = "Survival Time (Days)")
        
        # Add title
        surv_plot$plot <- surv_plot$plot + 
                          labs(title = paste(method, "Clusters Survival Analysis")) + 
                          theme(plot.title = element_text(hjust = 0.5))

        # Create a ggtext plot for the p-value matrix
        pval_df <- as.data.frame(as.table(pval_matrix))
        colnames(pval_df) <- c("Cluster1", "Cluster2", "PValue")
        pval_df$PValue <- round(pval_df$PValue, 3)
        pval_df <- na.omit(pval_df)
        pval_df$PValue <- ifelse(pval_df$PValue < 0.001, "< 0.001", as.character(pval_df$PValue))
        
        # Ensure PValue is numeric for scale_fill_gradient
        pval_df$PValue_num <- as.numeric(ifelse(pval_df$PValue == "< 0.001", 0.001, pval_df$PValue))

        pval_matrix_plot <- ggplot(pval_df, aes(x = Cluster1, y = Cluster2, fill = PValue_num)) +
          geom_tile(color = "black", lwd = 0.5, linetype = 1) + # Add boundary around the matrix
          geom_text(aes(label = PValue), color = "black") +
          scale_fill_gradient2(low = "green", mid = "lightgreen", high = "red", midpoint = 0.05, space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill") +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white", color = "white"), # Set background to white
            panel.background = element_rect(fill = "white", color = "white"), # Set panel background to white
            panel.grid.major = element_line(color = "grey80", size = 0.5) # Add grid lines
          ) +
          labs(title = "Pairwise P-Values", fill = "P-Value")

        # Combine the plots
        combined_plot <- cowplot::plot_grid(surv_plot$plot, pval_matrix_plot, ncol = 1, rel_heights = c(2, 1))

        # Save the combined plot
        filename <- file.path("output", "survival_analysis", paste0(method, "_survival_plot.png"))
        ggsave(filename, plot = combined_plot, width = 7, height = 7, dpi = 300)

      }}
      
      # Store p-values for the current cancer type
      df_methods[paste("K_", K, "_", tumor_type, sep = "")] <- pval_list
    }}
    
    # Save p-values to file
    write.table(df_methods, file = file.path("output", "survival_analysis", "p_values.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    ''')

# Main function
def main():
    # Load the data from the .pkl file
    with open(os.path.join('.', 'modules', 'data.pkl'), 'rb') as f:
        X_scaled, pca_df, cancer_type, k, algorithm_choice = pickle.load(f)

    # Extract selected cancer and K values
    selected_cancer = cancer_type
    K = k

    # Debugging print statements to verify inputs
    # print("Selected Cancer:", cancer_type)
    # print("Number of Clusters (K):", k)
    # print("Algorithm Choice:", algorithm_choice)

    # Convert algorithm_choice from string to integer
    algorithm_choice = int(algorithm_choice)
    algorithms = ["SNF", "KMeans", "Hierarchical", "SpectralClustering", "FuzzyCMeans", "All"]
    # print(f"Choosen algorithm: {algorithms[algorithm_choice-1]}")
    # Call the survival analysis function
    # survival_analysis(selected_cancer, K, algorithm_choice)

    if algorithm_choice == 6:
        print(f"{BOLD}🧬 Running Survival Analysis for Clusters: {RESET}SNF, KMeans, Hierarchical, SpectralClustering, and FuzzyCMeans")
        for i, algorithm in enumerate(algorithms[:-1], 1):  # Iterate over all algorithms except "All"
            # print(f"Choosen algorithm: {algorithm}")
            survival_analysis(selected_cancer, K, i)
    else:
        algorithm = algorithms[algorithm_choice-1]
        print(f"{BOLD}🧬 Running Survival Analysis for Cluster: {RESET}{algorithm}")
        survival_analysis(selected_cancer, K, algorithm_choice)

    print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
    print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW}Survival analysis results saved to survival_analysis folder.{RESET}")

if __name__ == "__main__":
    main()
