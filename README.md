# C3MOD - Cancer Clustering and Characterization using Multi Omics Data
C3MOD is a python tool that identify cancer subtypes in multi omics data using unsupervised clustering algorithms and perform the various analysis on the identified subtypes.
The main objective of C3MOD is to provide user a centralized tool with complete flexibility of customized analysis to decide the type of cancer, algorithm, algorithm parameters and analysis so that, user can decide
based on it's area of interest user can draw the conclusion based on the analysis it has performed.

# Requirements
This tool can be download from https://github.com/rajulshakya777/C3MOD.git 
Recommended python version is : >=3.10

# Steps for Running C3MOD
1. Clone the github repositry using 'git clone https://github.com/rajulshakya777/C3MOD.git' in your Local machine.
2. Create a virtual python environment : python3 -m venv <env>
3. Activate the environemtn : source <env>/bin/activate
4. Run main.py file present inside c3mod_lib directory : python main.py

# Running C3MOD 
## Example 1: Run main.py

Step 1: Python and R packages installation
- Install all the required Python and R packages required for analysis

Step 2: Data Pre-Processing
- User input : Select the cancer data from the menu 
- Cancer data : ['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LAML', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']
- Data Pre-Processing includes Data cleaning, standardization, Dimensionality reduction using PCA and Data Integration.
- Prints Analysis Summary and the Pre-Processing completion message

Step 3: Recommend K as per your choice of data and algorithm
- User input : Select the algorithm for which you want to find the recommended K from the menu

1. SNF
2. KMeans
3. Hierarchical Clustering
4. Spectral Clustering
5. Fuzzy C-Means
6. Run All

Recommend K based on Silhouti score and WCSS methods.
Silhouti Scores and WCSS Plots Will be saved to recommended_K directory inside output directory
Example : Based on generated Silhouti score and WCSS plot for KMeans algorithm, tool will recommend you the best K (number of cluster) to choose else
you can decide by your own choice also
  
<img src="output/recommended_K/KMeans_recommended_k_plots.png" width="450"/>

<img src="https://github.com/user-attachments/assets/4523b77f-0ed8-4e1f-a09c-a69aeb51f45c" width="300" height="300"/>

Step 4: User Input: Algorithm Selection and K (number of cluster)

1. SNF
2. KMeans
3. Hierarchical Clustering
4. Spectral Clustering
5. Fuzzy C-Means
6. Run All

You need to select the algorithm from the menu and enter any integer k (between 2 and 7)
If you choose SNF, you can also change some other algorithm parameters from the below menu or you can keep them default by selecting option 4 

1. Number of neighbors for constructing affinity matrix (Default is 20)
2. Hyperparameter alpha for constructing affinity matrix (Default is 0.5)
3. Number of iterations for SNF (Default is 15)
4. No (Keep the default values)

Based on your choice, algorithm will run and all the results will be saved in clustering_resutls directory inside output directory.



Based on your selection, algorithm will run 

To predict 3 clusters for cancer type “CHOL” with n_iter (Number of iterations) = 7, use the following command-

<img src="iCluF_modules/CHOL_Readme.png" width="450"/>

## Example 2: 
To predict 4 clusters for cancer type “ACC” with n_iter (Number of iterations) = 6, use the following command-

<img src="iCluF_modules/ACC_Readme.png" width="450"/>

The predicted clusters are saved in the output folder .data/output/Clusters

Please choose a cancer type from the following list of 30 cancers-
['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS'
, 'UVM']

# Running iCluF on different dataset
For running the algorithm on a dataset that is not in the list, the user needs to create a folder with the proper name in “./data/input_data/TCGA_data”.  The format of the data types should be the same as given in the different cancer types. 
