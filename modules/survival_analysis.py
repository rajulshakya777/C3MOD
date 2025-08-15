import os
import pickle
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict

# Try R-based environment first
try:
  import rpy2.robjects as robjects  # type: ignore
  from rpy2.robjects import r  # type: ignore
  R_AVAILABLE = True
except Exception as e:  # noqa: E722
  robjects = None  # type: ignore
  r = None  # type: ignore
  R_AVAILABLE = False
  _R_IMPORT_ERROR = e

# Try Python lifelines fallback
try:
  from lifelines import KaplanMeierFitter
  from lifelines.statistics import logrank_test
  LIFELINES_AVAILABLE = True
except Exception as e:  # noqa: E722
  KaplanMeierFitter = None  # type: ignore
  logrank_test = None  # type: ignore
  LIFELINES_AVAILABLE = False
  _LIFELINES_IMPORT_ERROR = e

SURVIVAL_ENABLED = R_AVAILABLE or LIFELINES_AVAILABLE


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

def _survival_analysis_r(selected_cancer, K, algorithm_choice):
  """Original R-based survival analysis."""
  r(f'cancers <- list("{selected_cancer}")')
  r(f'K <- {K}')
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


def _find_survival_file(cancer: str) -> str:
  """Locate survival file (supports both root and survival_data subdir)."""
  candidates = [
    os.path.join('data', 'input_data', 'TCGA_data', cancer, f'{cancer}.survival_UCal.tsv'),
    os.path.join('data', 'input_data', 'TCGA_data', cancer, 'survival_data', f'{cancer}.survival_UCal.tsv'),
  ]
  for path in candidates:
    if os.path.exists(path):
      return path
  raise FileNotFoundError(f'Survival file not found for {cancer}. Looked in: {candidates}')


def _python_survival_for_algorithm(cancer: str, K: int, algorithm: str) -> Dict[str, str]:
  """Run survival analysis using lifelines for one algorithm; returns produced file paths."""
  if not LIFELINES_AVAILABLE:
    raise RuntimeError(f'lifelines not installed: {_LIFELINES_IMPORT_ERROR}')
  survival_file = _find_survival_file(cancer)
  surv_df = pd.read_csv(survival_file, sep='\t', index_col=0)
  # Normalize sample IDs similar to R code (replace '-' with '.')
  surv_df.index = surv_df.index.str.replace('-', '.')
  # Needed columns
  required_cols = {'OS', 'OS.time'}
  if not required_cols.issubset(set(surv_df.columns)):
    raise ValueError(f'Survival file missing required columns {required_cols}')

  class_file = os.path.join('output', 'clustering_results', f'{cancer}_classification_{algorithm}.txt')
  if not os.path.exists(class_file):
    warnings.warn(f'Classification file missing: {class_file}', RuntimeWarning)
    return {}
  cls_df = pd.read_csv(class_file, sep='\t', index_col=0)
  # Ensure consistent naming
  if 'samples' not in cls_df.columns or 'cluster' not in cls_df.columns or 'K' not in cls_df.columns:
    raise ValueError('Classification file missing required columns (samples, cluster, K)')
  cls_df = cls_df[cls_df['K'] == K]
  cls_df = cls_df.set_index('samples')
  # Intersect
  common = surv_df.index.intersection(cls_df.index)
  if len(common) == 0:
    warnings.warn(f'No overlapping samples for survival analysis ({algorithm}).', RuntimeWarning)
    return {}
  merged = pd.DataFrame({
    'OS': surv_df.loc[common, 'OS'],
    'OS.time': surv_df.loc[common, 'OS.time'],
    'cluster': cls_df.loc[common, 'cluster']
  }).dropna()
  # Plot KM curves
  km_fig, ax = plt.subplots(figsize=(7, 5))
  kmf = KaplanMeierFitter()
  p_values = {}
  clusters = sorted(merged['cluster'].unique())
  for c in clusters:
    mask = merged['cluster'] == c
    kmf.fit(durations=merged.loc[mask, 'OS.time'], event_observed=merged.loc[mask, 'OS'], label=f'Cluster {c}')
    kmf.plot(ax=ax)
  ax.set_title(f'{algorithm} Clusters Survival (Python)')
  ax.set_xlabel('Time (Days)')
  ax.set_ylabel('Survival Probability')
  km_plot_path = os.path.join(output_dir, f'{algorithm}_survival_plot_python.png')
  km_fig.tight_layout()
  km_fig.savefig(km_plot_path, dpi=150)
  plt.close(km_fig)
  # Pairwise log-rank tests
  import numpy as np
  p_matrix = np.full((len(clusters), len(clusters)), np.nan)
  for i, ci in enumerate(clusters):
    for j, cj in enumerate(clusters):
      if j <= i:
        continue
      r1 = merged['cluster'] == ci
      r2 = merged['cluster'] == cj
      res = logrank_test(merged.loc[r1, 'OS.time'], merged.loc[r2, 'OS.time'],
                 event_observed_A=merged.loc[r1, 'OS'], event_observed_B=merged.loc[r2, 'OS'])
      p_matrix[i, j] = res.p_value
  # Heatmap
  fig_h, ax_h = plt.subplots(figsize=(5, 4))
  im = ax_h.imshow(p_matrix, cmap='viridis', vmin=0, vmax=0.1)
  ax_h.set_xticks(range(len(clusters)))
  ax_h.set_yticks(range(len(clusters)))
  ax_h.set_xticklabels(clusters)
  ax_h.set_yticklabels(clusters)
  ax_h.set_title('Pairwise Log-rank p-values')
  for i in range(len(clusters)):
    for j in range(len(clusters)):
      if not (j > i):
        continue
      val = p_matrix[i, j]
      if not (val is None or pd.isna(val)):
        ax_h.text(j, i, f'{val:.3g}', ha='center', va='center', color='white', fontsize=9)
  fig_h.colorbar(im, ax=ax_h, fraction=0.046, pad=0.04)
  heatmap_path = os.path.join(output_dir, f'{algorithm}_survival_pvalues_python.png')
  fig_h.tight_layout()
  fig_h.savefig(heatmap_path, dpi=150)
  plt.close(fig_h)
  return {'km_plot': km_plot_path, 'pvalue_heatmap': heatmap_path}


def survival_analysis(selected_cancer: str, K: int, algorithm_choice: int):
  """Dispatch to R or Python implementation based on availability.

  algorithm_choice follows original convention (1..5 or 6 for all).
  """
  algorithms = ["SNF", "KMeans", "Hierarchical", "SpectralClustering", "FuzzyCMeans", "All"]
  if R_AVAILABLE:
    _survival_analysis_r(selected_cancer, K, algorithm_choice)
    return
  if not LIFELINES_AVAILABLE:
    print(f"[WARN] Survival analysis skipped: neither R nor lifelines available (R error: {_R_IMPORT_ERROR}, lifelines error: {_LIFELINES_IMPORT_ERROR})")
    return
  # Python fallback
  if algorithm_choice == 6:
    algos_to_run = algorithms[:-1]  # all but 'All'
  else:
    algos_to_run = [algorithms[algorithm_choice - 1]]
  print(f"[INFO] Python survival analysis running for: {', '.join(algos_to_run)}")
  for algo in algos_to_run:
    try:
      _python_survival_for_algorithm(selected_cancer, K, algo)
    except Exception as e:
      print(f"[WARN] Survival analysis failed for {algo}: {e}")

# Main function
def main():
  if not SURVIVAL_ENABLED:
    print("[WARN] No survival backend available (install R + packages or 'lifelines').")
    return
  with open(os.path.join('.', 'modules', 'data.pkl'), 'rb') as f:
    X_scaled, pca_df, cancer_type, k, algorithm_choice = pickle.load(f)
  selected_cancer = cancer_type
  K = k
  algorithm_choice = int(algorithm_choice)
  # Run analysis
  survival_analysis(selected_cancer, K, algorithm_choice)
  print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
  print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW}Survival analysis results saved to survival_analysis folder.{RESET}")

if __name__ == "__main__":
    main()
