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

# Running iCluF 
## Example 1: 
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
