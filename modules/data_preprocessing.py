import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pickle
import os

# Color and style codes
YELLOW = "\033[1;33m"
CYAN = "\033[1;36m"
BLUE = "\033[1;34m"
BOLD = "\033[1m"
RESET = "\033[0m"
GREEN = "\033[1;32m"


CHECK_EMOJI = f"{GREEN}✔{RESET}"
# Global variables
X_scaled = None
pca_df = None

# Function to load and process data
def read_and_process_data(file_path):
    data = pd.read_csv(file_path, sep='\t')
    data = data.rename(columns={'Unnamed: 0': 'PatId'}).transpose().reset_index()
    data.columns = data.iloc[0]
    data = data[1:]
    return data

# Function to merge data
def merge_data(mRNA_data, miRNA_data, meth_data):
    return mRNA_data.merge(miRNA_data, on='PatId').merge(meth_data, on='PatId')

# Function to print dataframe dimensions
def print_data_info(df, name):
    total_rows, total_columns = df.shape
    print(f"{BLUE}{BOLD}{name}:{RESET}")
    print(f"  {YELLOW} • Total number of samples:{RESET} {total_rows}")
    print(f"  {YELLOW} • Total features:{RESET} {total_columns}")

# Function to standardize data
def standardize_data(features):
    print(f"{YELLOW}⚙️  Standardizing data...{RESET}")
    scaler = StandardScaler()
    return scaler.fit_transform(features)

# Function to perform PCA
def perform_pca(X_scaled):
    print(f"{YELLOW}⚙️  Performing Principal Component Analysis (PCA)...{RESET}")
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(X_scaled)
    return pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

# Function to get cancer dataset type from user input
def get_cancer_type():
    cancer_types = [
        'ACC', 'BRCA', 'BLCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 
        'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 
        'MESO', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'STAD', 'TGCT', 'THCA', 
        'THYM', 'UCEC', 'UCS', 'UVM'
    ]
    
    print(f"{BOLD}Select the cancer dataset:{RESET}")
    for idx, cancer_type in enumerate(cancer_types, 1):
        print(f"{idx}. TCGA-{cancer_type}")

    while True:
        try:
            choice = int(input("Enter the number corresponding to your cancer type choice: "))
            if 1 <= choice <= len(cancer_types):
                selected_cancer_type = cancer_types[choice - 1]
                print(f"{BOLD}🔍 Selected Cancer Type:{RESET} TCGA-{selected_cancer_type}")
                return selected_cancer_type
            else:
                print(f"{YELLOW}⚠️ Invalid choice. Please select a number between 1 and {len(cancer_types)}.{RESET}")
        except ValueError:
            print(f"{YELLOW}⚠️ Invalid input. Please enter a valid number.{RESET}")

# Main script logic
def main():
    global X_scaled, merged_data  # Declare global variables

    cancer_type = get_cancer_type()
    
    
    # Generate dynamic file paths based on cancer type input
    base_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'input_data', 'TCGA_data')
    mRNA_file_path = os.path.join(base_path, cancer_type, f'{cancer_type}_mRNA_GANfeature.txt')
    miRNA_file_path = os.path.join(base_path, cancer_type, f'{cancer_type}_miRNA_GANfeature.txt')
    meth_file_path = os.path.join(base_path, cancer_type, f'{cancer_type}_methyl_GANfeature.txt')
    
    # Load and process data
    print(f"  {CYAN}{BOLD} Reading and processing data...")
    mRNA_data = read_and_process_data(mRNA_file_path)
    miRNA_data = read_and_process_data(miRNA_file_path)
    meth_data = read_and_process_data(meth_file_path)

    # Merge data
    merged_data = merge_data(mRNA_data, miRNA_data, meth_data)

    # Print data info
    print(f"{BOLD}{BLUE}📊 Analysis Summary:{RESET}")
    print_data_info(merged_data, f"{YELLOW}🔗{BOLD} Integrated mRNA, miRNA, and methylation data{RESET}")

    # Drop 'PatId' and standardize data
    features = merged_data.drop(columns=['PatId'])
    X_scaled = standardize_data(features)
    print(f" {YELLOW}{CHECK_EMOJI} Data Standardization completed.{RESET}")

    # Perform PCA
    pca_df = perform_pca(X_scaled)
    pca_df['PatId'] = merged_data['PatId']
    print(f" {YELLOW}{CHECK_EMOJI} Principal Component Analysis (PCA) completed.{RESET}")

    # Save X_scaled, pca_df, and cancer_type

    output_path = os.path.join('.', 'modules', 'data.pkl')
   
    # Delete the file if it exists
    if os.path.exists(output_path):
        os.remove(output_path)

    # Save X_scaled, pca_df, and cancer_type
    with open(output_path, 'wb') as f:
        pickle.dump((X_scaled, pca_df, cancer_type), f)
    print(f"{CYAN}📈 {BOLD}Results:{RESET} {YELLOW}All results and plots saved to recommended_K folder.{RESET}")

if __name__ == "__main__":
    main()