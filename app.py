import os
import pickle
import streamlit as st
import pandas as pd
from modules import data_preprocessing
from modules import clustering as clustering_mod
from modules import survival_analysis as surv_mod
from modules import mutation_analysis as mut_mod
from modules import stage_analysis as stage_mod
from modules import immune_analysis as immune_mod

# NOTE: Existing modules are very file-output oriented. For first iteration we trigger them and then surface produced files.
# Later we can refactor each analysis to return in-memory figures / dataframes instead of only writing to disk.

DATA_PKL_PATH = os.path.join('modules', 'data.pkl')
OUTPUT_BASE = os.path.join('output')

ALGORITHMS = ["SNF", "KMeans", "Hierarchical", "SpectralClustering", "FuzzyCMeans"]
ANALYSES = {
    'Survival Analysis': 'survival',
    'Mutation Analysis': 'mutation',
    'Stage Analysis': 'stage',
    'Immune Analysis': 'immune'
}

def load_pickle():
    if os.path.exists(DATA_PKL_PATH):
        with open(DATA_PKL_PATH, 'rb') as f:
            return pickle.load(f)
    return None

def run_preprocessing(cancer_type: str):
    # Temporarily replace interactive get_cancer_type
    orig = data_preprocessing.get_cancer_type
    data_preprocessing.get_cancer_type = lambda: cancer_type
    try:
        data_preprocessing.main()
    finally:
        data_preprocessing.get_cancer_type = orig


def run_clustering(cancer_type: str, k: int, selected_algos, snf_params):
    # Load data produced by preprocessing
    with open(DATA_PKL_PATH, 'rb') as f:
        X_scaled, pca_df, cancer_type_loaded = pickle.load(f)[:3]
    # Save expanded tuple including k + algorithm choice placeholder (we keep algorithm choice for downstream modules)
    algo_choice_number = 6 if len(selected_algos) > 1 else (ALGORITHMS.index(selected_algos[0]) + 1)
    with open(DATA_PKL_PATH, 'wb') as f:
        pickle.dump((X_scaled, pca_df, cancer_type_loaded, k, algo_choice_number), f)
    result_files = clustering_mod.run_non_interactive_clustering(
        X_scaled, pca_df, cancer_type_loaded, k, selected_algos, snf_params=snf_params
    )
    # Collect PCA plot paths
    pca_plots = {}
    for algo in selected_algos:
        if algo == 'SNF':
            p = os.path.join('output', 'clustering_results', f'{cancer_type_loaded}_SNF_pca.png')
        else:
            p = os.path.join('output', 'clustering_results', f'pca_{algo}.png')
        if os.path.exists(p):
            pca_plots[algo] = p
    return result_files, pca_plots


def run_analysis(kind: str, algo_choice_number: int):
    # Each module reads modules/data.pkl for needed metadata; we ensure algo_choice is stored already.
    if kind == 'survival':
        surv_mod.main()
    elif kind == 'mutation':
        mut_mod.main()
    elif kind == 'stage':
        stage_mod.main()
    elif kind == 'immune':
        immune_mod.main()


def list_generated_files(subfolder):
    folder = os.path.join(OUTPUT_BASE, subfolder)
    paths = []
    for root, _, files in os.walk(folder):
        for f in files:
            paths.append(os.path.join(root, f))
    return sorted(paths)

st.set_page_config(page_title="C3MOD Multi-Omics Clustering", layout='wide')

st.title("C3MOD: Multi-Omics Clustering and Downstream Analyses")
st.markdown("This Streamlit interface wraps the existing pipeline. Select options, run steps, and view or download results.")

with st.expander("1. Data Preprocessing", expanded=True):
    cancer_type = st.selectbox("Select Cancer Type", [
        'ACC','BRCA','BLCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','PAAD','PCPG','PRAD','READ','SARC','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'
    ], index=16)  # default LUAD
    if st.button("Run Preprocessing", key="preprocess_btn"):
        with st.spinner("Preprocessing data..."):
            run_preprocessing(cancer_type)
        st.success("Preprocessing complete. Data stored.")
        st.session_state['preprocessed'] = True

preprocessed = st.session_state.get('preprocessed', os.path.exists(DATA_PKL_PATH))

with st.expander("2. Clustering", expanded=preprocessed):
    if not preprocessed:
        st.info("Run preprocessing first.")
    else:
        k = st.slider("Number of Clusters (K)", 2, 7, 3)
        selected_algos = st.multiselect("Algorithms", ALGORITHMS, default=['KMeans'])
        use_snf = 'SNF' in selected_algos
        snf_params = {}
        if use_snf:
            with st.popover("SNF Parameters"):
                snf_params['n_neighbour'] = st.number_input("Number of neighbors", 5, 25, 20)
                snf_params['alpha'] = st.slider("Alpha", 0.1, 0.9, 0.5, 0.1)
                snf_params['T'] = st.number_input("Iterations (T)", 5, 20, 15)
        if st.button("Run Clustering", key="cluster_btn"):
            with st.spinner("Running clustering algorithms..."):
                result_files, pca_plots = run_clustering(cancer_type, k, selected_algos, snf_params)
            st.session_state['clustering_done'] = True
            st.session_state['pca_plots'] = pca_plots
            st.session_state['classification_files'] = result_files
            st.success("Clustering complete.")
            for algo, path in pca_plots.items():
                st.image(path, caption=f"PCA Plot - {algo}")

clustering_done = st.session_state.get('clustering_done', False)

with st.expander("3. Downstream Analyses", expanded=clustering_done):
    if not clustering_done:
        st.info("Run clustering first.")
    else:
        analyses_to_run = st.multiselect("Select analyses to run", list(ANALYSES.keys()))
        if st.button("Run Selected Analyses", key="analyses_btn"):
            # Retrieve algorithm choice number from pickle (already set during clustering)
            with st.spinner("Running analyses..."):
                data_tuple = load_pickle()
                if data_tuple and len(data_tuple) >= 5:
                    _, _, _, _, algo_choice_number = data_tuple
                    for label in analyses_to_run:
                        run_analysis(ANALYSES[label], algo_choice_number)
            st.success("Analyses complete.")

with st.expander("4. Results Browser", expanded=False):
    section = st.selectbox("Section", ["clustering_results", "survival_analysis", "mutation_analysis", "stage_analysis", "immune_analysis"])    
    files = list_generated_files(section)
    if not files:
        st.info("No files yet. Run steps above.")
    else:
        for fpath in files:
            rel = os.path.relpath(fpath)
            st.write(rel)
            if any(rel.lower().endswith(ext) for ext in ['.png', '.jpg', '.jpeg']):
                st.image(rel)
            if rel.lower().endswith(('.txt', '.tsv', '.csv')):
                try:
                    if rel.endswith('.tsv') or rel.endswith('.txt'):
                        df = pd.read_csv(rel, sep='\t')
                    else:
                        df = pd.read_csv(rel)
                    st.dataframe(df.head(100))
                except Exception as e:
                    st.caption(f"Could not display table: {e}")
            with open(rel, 'rb') as f:
                st.download_button(label="Download", data=f, file_name=os.path.basename(rel))

st.markdown("---")
st.caption("C3MOD Streamlit interface (initial migration). Further refactoring can return in-memory objects instead of reading from disk.")
