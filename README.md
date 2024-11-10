# C3MOD - Cancer Clustering and Characterization using Multi-Omics Data

**C3MOD** is a Python tool that identifies cancer subtypes in multi-omics data using unsupervised clustering algorithms and enables a variety of analyses on the identified subtypes. It provides flexibility for customized analysis, allowing users to select the cancer type, algorithm, parameters, and analyses that suit their research focus.

---

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Input Data](#input-data)
- [Running C3MOD](#running-c3mod)
- [Detailed Analysis Workflow](#detailed-analysis-workflow)
  - [Survival Analysis](#survival-analysis)
  - [Mutation Analysis](#mutation-analysis)
  - [Stage Analysis](#stage-analysis)
  - [Immune Analysis](#immune-analysis)
- [License](#license)

---

## Requirements
- **Python Version**: >=3.10

## Installation

1. Clone the GitHub repository:
   ```bash
   git clone https://github.com/rajulshakya777/C3MOD.git
2. Create a virtual python environment:
   ```bash
   python3 -m venv <env>
3. Activate the virtual environment:
   ```bash
   source <env>/bin/activate
4. Run the main script:
   ```bash
   python c3mod_lib/main.py

## Input Data
All input data and analysis scripts are included in the cloned repository.

# Running C3MOD

## Example 1: Run `main.py`script

### Step 1: Python and R Packages Installation
- Install all required Python and R packages necessary for the analysis.

### Step 2: Data Pre-Processing
- **User Input**: Select cancer data type from the menu.
- **Available Cancer Data**: 
  `['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LAML', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']`
- **Processing Steps**: Data cleaning, standardization, PCA for dimensionality reduction, and data integration.
- **Output**: Summary message confirming pre-processing completion.

### Step 3: K Recommendation (Cluster Count)
- **User Input**: Choose the algorithm for which to find the optimal K value:
  1. SNF
  2. KMeans
  3. Hierarchical Clustering
  4. Spectral Clustering
  5. Fuzzy C-Means
  6. Run All

- **Recommendation Method**: Based on silhouette scores and WCSS (within-cluster sum of squares) methods. Plots are saved in the `output/recommended_K` directory.

  **Example**: Based on the silhouette score and WCSS plot for KMeans, the tool recommends an optimal K (cluster count) to use.  
  <img src="output/recommended_K/KMeans_recommended_k_plots.png" width="500"/>

### Step 4: Algorithm Selection and Clustering
- **User Input**: Select the clustering algorithm and specify K (cluster count).
  - **Available Algorithms**:
    1. SNF
    2. KMeans
    3. Hierarchical Clustering
    4. Spectral Clustering
    5. Fuzzy C-Means
    6. Run All

- **Optional Parameters for SNF**:
  1. Number of neighbors for constructing affinity matrix (Default: 20)
  2. Hyperparameter alpha for affinity matrix (Default: 0.5)
  3. Number of SNF iterations (Default: 15)
  4. Keep default values

- **Output**: Results are saved in the `output/clustering_results` directory.
  1. **Classification File**
     <img src="https://github.com/user-attachments/assets/dac1096d-e61a-40c3-8c88-3ab067632922"/>
  2. **PCA Plots**
     <img src="https://github.com/user-attachments/assets/4058aeb4-b72f-4ac3-819c-3dbe41bea426"/>

### Step 5: Further Analysis Selection
- Choose the type of analysis to perform:
  1. Survival Analysis
  2. Mutation Analysis
  3. Stage Analysis
  4. Immune Analysis
  5. Run All

- **Example**: Selecting option 5 runs all analyses sequentially.

---

## Detailed Analysis Workflow

### Survival Analysis
- **Method**: Cox Proportional-Hazards Model (coxph).
- **Output**: Results are saved in the `output/survival_analysis` directory.

  **KMeans Survival Plot**:
  <img src="https://github.com/user-attachments/assets/3a4d8816-a7cb-4851-9d90-21dd24b43b41" width="500"/>

  **Interpretation**:
  1. **Cluster-Specific Survival Trends**: Each cluster shows unique survival curves with varying survival probabilities.
  2. **Variability**: Shaded areas represent the range of survival times within each cluster.
  3. **Significance Testing**: Pairwise p-values show statistical significance (low values) or similarity (high values) between clusters.
  4. **Cluster Comparison**: Survival curves and p-values highlight differences in survival outcomes.
  5. **Potential Clinical Insights**: Variations in survival suggest possible biological or clinical differences among clusters.

### Mutation Analysis
- **Input Data**: Identified clusters and gene mutation data.
- **Analyses Performed**:
  1. **Chi-Square Test**: Significance test for each gene between pairs of clusters.
     <img src="https://github.com/user-attachments/assets/df9509f8-8255-4b23-ae79-98570b8f8038"/>

     **Interpretation**: Example - For gene CPHL1P, a p-value of 0.99941 between clusters 1 and 2 indicates a low significance difference.

  2. **Top 5 Significant Genes**: Genes with the lowest p-values between cluster pairs.
     <img src="https://github.com/user-attachments/assets/341c9165-eae5-41c0-bc94-6eb21ac6b4f0"/>

     **Interpretation**: Genes such as BMP8A, FCGR3B, REEP6, IGKV1D-12, and GNG4 show the lowest p-values, meaning they are highly significant between cluster pairs 1 and 3.

  3. **Cluster Similarity**: Correlation matrix shows mutation similarity between clusters.
     <img src="https://github.com/user-attachments/assets/d664de49-33df-433e-835d-eaa9f5645762" width="500"/>

     **Interpretation**:
     - **High Correlation Values**: Correlation coefficients close to 1 suggest strong mutation pattern similarities.
     - **Statistical Significance**: Extremely low p-values (e.g., 1.000e-100) indicate high statistical significance.
     - **Cluster Similarities**: Clusters 1 and 3, as well as clusters 2 and 3, show high correlation values (0.871 and 0.749), suggesting similar mutation patterns.

  4. **Total Mutations per Cluster**: Shows total mutations for each cluster.
     <img src="https://github.com/user-attachments/assets/beb06da6-211a-4c5e-9623-f7d559f8edf1" width="500"/>

     **Interpretation**: Distinct mutation counts suggest significant differences between clusters based on mutation frequency.

  5. **Genes vs. P-Values**: Displays gene significance across cluster pairs.
     <img src="https://github.com/user-attachments/assets/f0950ee1-f3e6-4ead-8a1c-8cc577327417" width="500"/>

     **Interpretation**: Few genes in clusters 1 and 2 (shown in red) have a p-value < 0.05, indicating low evidence for significant differences between these clusters.

### Stage Analysis
- **Objective**: Show stage distribution among patients across clusters.
- **Input Data**: Stage data and identified clusters.
  
  <img src="https://github.com/user-attachments/assets/8f6e7999-3961-4a60-808c-24177edc768f" width="450"/>

  **Interpretation**: The plot shows how the stage distribution varies across clusters. Such insights can inform stage-based cluster analysis.

### Immune Analysis
- **Input Data**: Immune cell data and cluster IDs.
- **Output**: Bar plots showing immune-related characteristics of each cluster.
  <img src="https://github.com/user-attachments/assets/c7982f96-8a10-4ca9-a299-cb1c5a32fd33" width="500"/>

---

## License
C3MOD is licensed under the MIT License. For full terms and conditions, please refer to the [LICENSE](LICENSE) file in the repository.

---

## Contributing
We welcome contributions from the community! Please fork the repository, make your changes, and submit a pull request.

## Acknowledgments
C3MOD leverages cutting-edge clustering techniques, genomic data processing, and statistical analysis to identify clinically relevant patterns in cancer data. Special thanks to all contributors and collaborators.

---  

## Report Issue
In case of errors or improvements, please raise an issue via GitHub.

