import os
import sys
import io
import zipfile
import tempfile
import datetime
import time
import uuid
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
# Detect R availability via clustering module flag if present
try:
    from modules.clustering import R_AVAILABLE as R_OK, _R_IMPORT_ERROR  # type: ignore
except Exception:
    R_OK = False
    _R_IMPORT_ERROR = None
try:
    from modules.survival_analysis import SURVIVAL_ENABLED as SURV_OK  # type: ignore
except Exception:
    SURV_OK = False
ANALYSES = {
    'Survival Analysis': 'survival',
    'Mutation Analysis': 'mutation',
    'Stage Analysis': 'stage',
    'Immune Analysis': 'immune'
}
if not SURV_OK:
    ANALYSES = {k: v for k, v in ANALYSES.items() if v != 'survival'}

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
        # Set session state after successful preprocessing
        st.session_state['preprocessing_completed'] = True
        st.session_state['show_preprocessing_success'] = True
        print('[INFO] Preprocessing completed successfully, session state updated')
    finally:
        data_preprocessing.get_cancer_type = orig


def run_clustering(cancer_type: str, k: int, selected_algos, snf_params):
    # Validate input
    if not selected_algos:
        raise ValueError("No algorithms selected for clustering")
    
    print(f"[DEBUG] Selected algorithms: {selected_algos}")
    print(f"[DEBUG] Available algorithms: {ALGORITHMS}")
    
    # Load data produced by preprocessing
    with open(DATA_PKL_PATH, 'rb') as f:
        X_scaled, pca_df, cancer_type_loaded = pickle.load(f)[:3]
    # Save expanded tuple including k + algorithm choice placeholder (we keep algorithm choice for downstream modules)
    if len(selected_algos) > 1:
        algo_choice_number = 6
    else:
        try:
            algo_choice_number = ALGORITHMS.index(selected_algos[0]) + 1
        except (ValueError, IndexError) as e:
            print(f"[ERROR] Algorithm '{selected_algos[0]}' not found in ALGORITHMS list. Error: {e}")
            # Fallback if algorithm not found in list
            algo_choice_number = 1
    with open(DATA_PKL_PATH, 'wb') as f:
        pickle.dump((X_scaled, pca_df, cancer_type_loaded, k, algo_choice_number), f)
    result_files, metrics = clustering_mod.run_non_interactive_clustering(
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
    return result_files, pca_plots, metrics


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

st.set_page_config(page_title="C3MOD - Cancer Clustering and Characterization", layout='wide', page_icon="🧬")

# Initialize session state for cloud deployment compatibility
if 'app_initialized' not in st.session_state:
    st.session_state['app_initialized'] = True
    st.session_state['preprocessing_completed'] = False
    st.session_state['clustering_done'] = False
    st.session_state['analyses_done'] = False
    st.session_state['show_preprocessing_success'] = False
    st.session_state['show_clustering_success'] = False
    st.session_state['show_analyses_success'] = False

# ----------------- Theming / CSS -----------------
CUSTOM_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Poppins:wght@300;400;500;600;700;800&family=Source+Sans+3:wght@300;400;500;600;700&display=swap');
html, body, .stApp {
  font-family: 'Source Sans 3', -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
  font-weight: 400;
  line-height: 1.6;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
/* Beautiful headings */
h1, h2, h3, h4, h5, h6, .stMarkdown h1, .stMarkdown h2, .stMarkdown h3, .stMarkdown h4, .stMarkdown h5, .stMarkdown h6 {
  font-family: 'Poppins', sans-serif !important;
  font-weight: 600 !important;
  letter-spacing: -0.025em !important;
  line-height: 1.3 !important;
}
/* Enhanced button typography */
button, .stButton > button {
  font-family: 'Poppins', sans-serif !important;
  font-weight: 600 !important;
  letter-spacing: 0.025em !important;
}
/* Tab and label improvements */
button[role="tab"], label {
  font-family: 'Poppins', sans-serif !important;
  font-weight: 500 !important;
}
body, .main, .stApp {background:linear-gradient(180deg,#ffffff 0%, #f5f9fc 70%); color:#12394d;}
/* Design Tokens */
:root {--c3-accent:#0092b9; --c3-accent-alt:#00b28e; --c3-accent-gradient:linear-gradient(135deg,var(--c3-accent),var(--c3-accent-alt)); --c3-radius:14px; --c3-radius-pill:999px; --c3-border:#c2dde8; --c3-bg-soft:#f0f6fa; --c3-shadow:0 4px 14px -4px rgba(20,60,90,0.15),0 2px 4px rgba(30,70,100,0.08);}
/* Smooth entrance */
.c3-fade-in {animation:c3fade .55s ease-out both;}
@keyframes c3fade {0%{opacity:0; transform:translateY(10px);}100%{opacity:1; transform:translateY(0);}}
/* Header */
.c3mod-header {padding:1.6rem 1.5rem 1.25rem 1.5rem; margin-bottom:1rem; background:linear-gradient(135deg,#ffffff,#f0f7fb); border:1px solid #d9e6ef; border-radius:22px; box-shadow:0 8px 24px -10px rgba(20,60,90,0.18), 0 2px 4px rgba(30,70,100,0.08); position:relative; overflow:hidden;}
.c3mod-header:before {content:""; position:absolute; inset:0; background:radial-gradient(circle at 18% 22%,rgba(90,160,200,0.15),rgba(255,255,255,0) 60%); opacity:0.9; pointer-events:none;}
.c3mod-header:after {content:""; position:absolute; top:-60%; left:-30%; width:160%; height:220%; background:conic-gradient(from 140deg,rgba(90,160,200,0.15),rgba(255,255,255,0) 55%); animation:shine 12s linear infinite; mix-blend-mode:overlay; opacity:0.55;}
@keyframes shine {0%{transform:rotate(0deg);}100%{transform:rotate(360deg);}}
.c3mod-header h1 {
  background:linear-gradient(90deg,#004d7a,#008793,#00bf72,#00bf72 70%,#008793 90%); 
  -webkit-background-clip:text; 
  color:transparent; 
  font-size:1.8rem; 
  margin:0; 
  letter-spacing:-0.02em; 
  font-weight:700; 
  font-family: 'Poppins', sans-serif;
  text-shadow:0 1px 2px rgba(0,0,0,0.08);
  line-height: 1.1;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
} 
.c3mod-sub {
  color:#2e657c; 
  font-size:1rem; 
  margin-top:0.6rem; 
  font-weight:400; 
  letter-spacing:0.01em; 
  display:flex; 
  flex-wrap:wrap; 
  gap:0.6rem; 
  align-items:center;
  font-family: 'Source Sans 3', sans-serif;
  line-height: 1.5;
}
.c3mod-badge {
  background:linear-gradient(120deg,#e2f4fb,#d3edf6); 
  padding:0.35rem 0.8rem; 
  border-radius:999px; 
  font-size:0.7rem; 
  letter-spacing:0.05em; 
  font-weight:600; 
  color:#0d4a63; 
  border:1px solid #c2dde8; 
  box-shadow:0 2px 4px -1px rgba(0,0,0,0.08);
  font-family: 'Poppins', sans-serif;
  text-transform: uppercase;
} 
.c3mod-badge.dim {
  background:#eef7fb; 
  color:#4c7e92;
  font-weight: 500;
}
.c3mod-pulse {position:relative;}
.c3mod-pulse:after {content:""; position:absolute; inset:0; border-radius:inherit; box-shadow:0 0 0 0 rgba(0,140,180,0.4); animation:pulse 3.2s ease-in-out infinite;}
@keyframes pulse {0%{box-shadow:0 0 0 0 rgba(0,140,180,0.4);}70%{box-shadow:0 0 0 14px rgba(0,140,180,0);}100%{box-shadow:0 0 0 0 rgba(0,140,180,0);}}
/* Tabs */
.stTabs [data-baseweb="tab-list"] {gap:0.5rem; margin-bottom:1.5rem; padding:0.4rem; background:rgba(255,255,255,0.8); border-radius:16px; border:1px solid #e2eff5; box-shadow:0 2px 8px -2px rgba(20,60,90,0.08);}
button[role="tab"] {
  background:linear-gradient(145deg,#ffffff,#f8fcfe) !important; 
  color:#1f5b72 !important; 
  border-radius:12px !important; 
  border:1px solid #d4e5ed !important; 
  font-weight:600 !important;
  font-size:0.95rem !important;
  padding:0.75rem 1.25rem !important;
  letter-spacing:0.3px !important;
  position:relative !important;
  transition:all 0.3s cubic-bezier(0.4,0.2,0.2,1) !important;
  min-height:48px !important;
  display:flex !important;
  align-items:center !important;
  justify-content:center !important;
  box-shadow:0 2px 4px -1px rgba(0,0,0,0.06) !important;
}
button[role="tab"]:hover {
  background:linear-gradient(145deg,#f1f9fc,#e8f4f9) !important;
  border-color:#b9d6e2 !important;
  transform:translateY(-2px) !important;
  box-shadow:0 4px 12px -3px rgba(20,60,90,0.15) !important;
}
button[role="tab"][aria-selected="true"] {
  background:linear-gradient(135deg,#0092b9,#00b28e) !important; 
  color:#ffffff !important; 
  border-color:#009fb3 !important;
  box-shadow:0 0 0 1px #5ed3ec inset, 0 6px 16px -4px rgba(0,120,140,0.4), 0 2px 4px rgba(0,0,0,0.08) !important;
  transform:translateY(-1px) !important;
}
button[role="tab"][aria-selected="true"]:hover {
  filter:brightness(1.05) !important;
  transform:translateY(-3px) !important;
}
/* Tab content area */
.stTabs [data-baseweb="tab-panel"] {
  padding:1.5rem 0 0 0 !important;
  animation:tabFadeIn 0.4s ease-out both;
}
@keyframes tabFadeIn {
  0% {opacity:0; transform:translateY(8px);}
  100% {opacity:1; transform:translateY(0);}
} 
/* Cards */
.c3-card, .glass {background:rgba(255,255,255,0.7); backdrop-filter:blur(8px); -webkit-backdrop-filter:blur(8px); padding:1rem 1.15rem; border-radius:16px; border:1px solid #d4e5ed; box-shadow:0 4px 14px -4px rgba(20,60,90,0.15);} 
.mini-card {display:inline-block; padding:0.35rem 0.65rem; background:#e3f4fa; border:1px solid #c2dde8; border-radius:999px; font-size:0.72rem; font-weight:600; letter-spacing:0.5px; margin:0.25rem 0.35rem 0.25rem 0; color:#0f4d63;}
.algo-badge {cursor:pointer; transition:all .18s ease;}
.algo-badge:hover {background:#d5ecf5; transform:translateY(-2px); box-shadow:0 4px 10px -4px rgba(20,60,90,0.15);} 
.algo-selected {background:linear-gradient(125deg,#00a4b7,#00bf72); border-color:#00b38a; color:#fff; box-shadow:0 0 0 1px #6ae3ff inset, 0 0 0 2px rgba(255,255,255,0.5);} 
/* Metrics */
.metrics-grid {display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:0.8rem; margin-top:0.6rem;}
.metric-box {padding:0.7rem 0.85rem; background:#f2f9fc; border:1px solid #c2dde8; border-radius:14px; position:relative;}
.metric-box.good {border-color:#41b07a; background:linear-gradient(145deg,#e9f8f2,#d5f2e6);} 
.metric-box.ok {border-color:#d1b34b; background:linear-gradient(145deg,#faf6e5,#f2e6c5);} 
.metric-box.bad {border-color:#d46a6a; background:linear-gradient(145deg,#fdecec,#f9d6d6);} 
.metric-box h5 {margin:0 0 0.25rem 0; font-size:0.8rem; text-transform:uppercase; letter-spacing:0.5px; font-weight:600; opacity:0.75;}
.metric-val {font-size:1.15rem; font-weight:600; color:#0d4a63;}
/* Dataframe */
div[data-testid="stDataFrame"] {background:#ffffff; border:1px solid #d4e5ed; border-radius:14px;}
/* Sidebar */
section[data-testid="stSidebar"] {background:linear-gradient(180deg,#ffffff,#f1f7fa); border-right:1px solid #d9e6ef;}
/* Download buttons spacing */
div[data-testid="stDownloadButton"] {margin-bottom:0.55rem;}
/* Log */
.log-box {font-family:monospace; white-space:pre-wrap; background:#082b38; color:#e8f8fd; padding:0.85rem 1rem; border:1px solid #0d4356; border-radius:14px; max-height:420px; overflow:auto; font-size:0.72rem; line-height:1.2;}
/* Footer */
.footer-note {font-size:0.65rem; letter-spacing:0.4px; opacity:0.55;}
/* Primary Buttons */
div.stButton > button {background:var(--c3-accent-gradient); color:#fff; border:1px solid #009fb3; padding:0.55rem 1.15rem; font-weight:600; letter-spacing:0.3px; border-radius:12px; box-shadow:0 4px 10px -2px rgba(0,120,140,0.35),0 2px 4px rgba(0,0,0,0.08); transition:all .25s ease;}
div.stButton > button:hover {filter:brightness(1.05); transform:translateY(-2px); box-shadow:0 6px 16px -4px rgba(0,120,140,0.4);}
div.stButton > button:active {transform:translateY(0); box-shadow:0 2px 6px -2px rgba(0,120,140,0.3);}
/* Secondary Buttons / Download */
div[data-testid="stDownloadButton"] button, div.stButton.secondary > button {background:#ffffff; color:#0d4a63; border:1px solid #b9d6e2; border-radius:10px; font-weight:500; box-shadow:0 2px 6px -2px rgba(0,0,0,0.08); transition:all .25s ease;}
div[data-testid="stDownloadButton"] button:hover, div.stButton.secondary > button:hover {background:#f1f9fc; border-color:#86bfd2;}
/* Checkbox / Slider refinement */
div.stCheckbox > label, div.stSlider > label {font-weight:600; color:#0d4a63;}
/* Inputs */
div[data-baseweb="select"] > div {border-radius:10px !important;}
input, textarea {border-radius:10px !important;}
/* Expander */
details[data-testid="stExpander"] {border:1px solid #d4e5ed; border-radius:12px; background:#ffffff;}
/* Animations for algorithm badges */
.algo-selected {animation:algopulse 2.4s ease-in-out infinite alternate;}
@keyframes algopulse {0%{filter:brightness(1);}100%{filter:brightness(1.1);}}
/* Refresh/Clear Button Styling */
section[data-testid="stSidebar"] button[data-testid="baseButton-secondary"],
section[data-testid="stSidebar"] div.stButton button[kind="secondary"],
section[data-testid="stSidebar"] button:has-text("Refresh & Clear Analysis") {
  background: linear-gradient(135deg, #fbbf24, #f59e0b) !important;
  color: #92400e !important;
  border: 1px solid #f59e0b !important;
  border-radius: 12px !important;
  font-weight: 600 !important;
  font-family: 'Poppins', sans-serif !important;
  letter-spacing: 0.025em !important;
  padding: 0.7rem 1.2rem !important;
  transition: all 0.3s cubic-bezier(0.4, 0.2, 0.2, 1) !important;
  box-shadow: 0 4px 12px -2px rgba(251, 191, 36, 0.3) !important;
  position: relative !important;
  overflow: hidden !important;
  width: 100% !important;
}
section[data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:before,
section[data-testid="stSidebar"] div.stButton button[kind="secondary"]:before {
  content: "";
  position: absolute;
  top: 0;
  left: -100%;
  width: 100%;
  height: 100%;
  background: linear-gradient(90deg, transparent, rgba(255,255,255,0.3), transparent);
  transition: left 0.6s ease;
}
section[data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:after,
section[data-testid="stSidebar"] div.stButton button[kind="secondary"]:after {
  content: "⚠️ Reset Pipeline\AThis will clear all analysis data and reset the pipeline to start fresh.";
  white-space: pre-line;
  position: absolute;
  top: 100%;
  left: 50%;
  transform: translateX(-50%);
  background: linear-gradient(135deg, #1f2937, #374151);
  color: #f9fafb;
  padding: 0.8rem 1rem;
  border-radius: 8px;
  font-size: 0.75rem;
  line-height: 1.4;
  min-width: 200px;
  text-align: center;
  box-shadow: 0 8px 20px -4px rgba(0,0,0,0.4);
  opacity: 0;
  visibility: hidden;
  transition: all 0.3s ease;
  z-index: 1000;
  border: 1px solid #4b5563;
}
section[data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:hover,
section[data-testid="stSidebar"] div.stButton button[kind="secondary"]:hover {
  background: linear-gradient(135deg, #d97706, #b45309) !important;
  color: #ffffff !important;
  border-color: #d97706 !important;
  transform: translateY(-3px) !important;
  box-shadow: 0 6px 18px -2px rgba(217, 119, 6, 0.4) !important;
}
section[data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:hover:before,
section[data-testid="stSidebar"] div.stButton button[kind="secondary"]:hover:before {
  left: 100%;
}
section[data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:hover:after,
section[data-testid="stSidebar"] div.stButton button[kind="secondary"]:hover:after {
  opacity: 1;
  visibility: visible;
  top: calc(100% + 8px);
}
/* Pipeline Status Styling */
section[data-testid="stSidebar"] h2 {
  background: linear-gradient(135deg, #ffffff, #f8fafc) !important;
  color: #1e293b !important;
  padding: 1rem 1.2rem !important;
  margin: -1rem -1rem 1.5rem -1rem !important;
  border-radius: 0 0 16px 16px !important;
  font-family: 'Poppins', sans-serif !important;
  font-weight: 600 !important;
  font-size: 1.1rem !important;
  letter-spacing: 0.025em !important;
  box-shadow: 0 4px 12px -4px rgba(30, 41, 59, 0.15) !important;
  border-bottom: 2px solid #3b82f6 !important;
  border: 1px solid #e2e8f0 !important;
  border-top: none !important;
}
.status-card {
  background: linear-gradient(135deg, #ffffff, #f8fafc);
  border: 1px solid #e2e8f0;
  border-radius: 12px;
  padding: 0.8rem 1rem;
  margin: 0.5rem 0;
  box-shadow: 0 2px 8px -2px rgba(0,0,0,0.08);
  transition: all 0.3s ease;
  position: relative;
  overflow: hidden;
}
.status-card:before {
  content: "";
  position: absolute;
  left: 0;
  top: 0;
  height: 100%;
  width: 4px;
  background: linear-gradient(180deg, #94a3b8, #cbd5e1);
  transition: all 0.3s ease;
}
.status-card.complete:before {
  background: linear-gradient(180deg, #22c55e, #16a34a);
}
.status-card.pending:before {
  background: linear-gradient(180deg, #f59e0b, #d97706);
}
.status-card:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 12px -2px rgba(0,0,0,0.12);
}
.status-label {
  font-family: 'Poppins', sans-serif;
  font-weight: 600;
  font-size: 0.85rem;
  color: #334155;
  margin-bottom: 0.2rem;
  letter-spacing: 0.025em;
}
.status-value {
  font-family: 'Source Sans 3', sans-serif;
  font-weight: 500;
  font-size: 0.8rem;
  color: #64748b;
  display: flex;
  align-items: center;
  gap: 0.4rem;
}
.status-value.complete {
  color: #16a34a;
}
.status-value.pending {
  color: #d97706;
}
.status-icon {
  font-size: 1rem;
  animation: statusPulse 2s ease-in-out infinite;
}
@keyframes statusPulse {
  0%, 100% { opacity: 1; }
  50% { opacity: 0.6; }
}
section[data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:active,
section[data-testid="stSidebar"] div.stButton button[kind="secondary"]:active {
  transform: translateY(-1px) !important;
  box-shadow: 0 4px 10px -2px rgba(217, 119, 6, 0.35) !important;
}
/* Success Messages */
.success-message {
  background: linear-gradient(135deg, #dcfce7, #bbf7d0);
  border: 1px solid #86efac;
  border-radius: 12px;
  padding: 1.2rem 1.4rem;
  margin: 1.2rem 0;
  box-shadow: 0 6px 16px -4px rgba(34, 197, 94, 0.3), 0 2px 4px rgba(34, 197, 94, 0.1);
  animation: successSlideIn 0.6s ease-out both, successPulse 1.0s ease-in-out 0.6s;
  position: relative;
  overflow: hidden;
}
.success-message::before {
  content: "";
  position: absolute;
  top: 0;
  left: -100%;
  width: 100%;
  height: 100%;
  background: linear-gradient(90deg, transparent, rgba(255,255,255,0.3), transparent);
  animation: successShine 3s ease-in-out 0.8s;
}
.success-message .icon {
  font-size: 1.3rem;
  margin-right: 0.6rem;
  line-height: 1;
  margin-top: 0.1rem;
  animation: iconBounce 0.8s ease-out 0.3s both;
}
@keyframes successSlideIn {
  0% {
    opacity: 0;
    transform: translateY(-20px) scale(0.95);
  }
  50% {
    opacity: 0.7;
    transform: translateY(-5px) scale(0.98);
  }
  100% {
    opacity: 1;
    transform: translateY(0) scale(1);
  }
}
@keyframes successPulse {
  0%, 100% {
    box-shadow: 0 6px 16px -4px rgba(34, 197, 94, 0.3), 0 2px 4px rgba(34, 197, 94, 0.1);
  }
  50% {
    box-shadow: 0 8px 20px -4px rgba(34, 197, 94, 0.4), 0 4px 8px rgba(34, 197, 94, 0.2);
  }
}
@keyframes successFadeOut {
  0% {
    opacity: 1;
    transform: translateY(0) scale(1);
  }
  100% {
    opacity: 0;
    transform: translateY(-10px) scale(0.98);
  }
}
@keyframes successShine {
  0% {
    left: -100%;
  }
  100% {
    left: 100%;
  }
}
@keyframes iconBounce {
  0% {
    transform: scale(0.8);
  }
  50% {
    transform: scale(1.1);
  }
  100% {
    transform: scale(1);
  }
}
.success-message .title {
  font-weight: 600;
  color: #15803d;
  font-size: 1.1rem;
  margin: 0;
  font-family: 'Poppins', sans-serif;
  letter-spacing: -0.01em;
}
.success-message .description {
  color: #166534;
  font-size: 0.95rem;
  margin: 0.35rem 0 0 0;
  opacity: 0.9;
  font-family: 'Source Sans 3', sans-serif;
  line-height: 1.5;
  font-weight: 400;
}

/* Increase font size for form elements */
.stSelectbox > div > div > div {
  font-size: 1.1rem !important;
  font-weight: 500 !important;
}

.stSelectbox label {
  font-size: 1.15rem !important;
  font-weight: 600 !important;
}

.stSlider > div > div > div > div {
  font-size: 1.1rem !important;
  font-weight: 500 !important;
}

.stSlider label {
  font-size: 1.15rem !important;
  font-weight: 600 !important;
}

.stMultiSelect > div > div > div {
  font-size: 1.1rem !important;
  font-weight: 500 !important;
}

.stMultiSelect label {
  font-size: 1.15rem !important;
  font-weight: 600 !important;
}

/* Style the dropdown options */
.stSelectbox > div > div > div > div {
  font-size: 1.1rem !important;
}

/* Style multiselect options */
.stMultiSelect > div > div > div > div {
  font-size: 1.1rem !important;
}
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# Hard refresh status calculation for header badge
_preprocessing_done = bool(os.path.exists(DATA_PKL_PATH)) or bool(st.session_state.get('preprocessing_completed', False))
_clustering_done = bool(st.session_state.get('clustering_done'))
_analyses_done = bool(st.session_state.get('analyses_done'))
_active_states = [int(_preprocessing_done), int(_clustering_done), int(_analyses_done)]
_completed = sum(_active_states)

# Debug logging for hard refresh
print(f'[HEADER] Hard refresh status: Preprocessing={_preprocessing_done}, Clustering={_clustering_done}, Analyses={_analyses_done} (Total: {_completed}/3)')
print(f'[DEBUG] Session state: preprocessing_completed={st.session_state.get("preprocessing_completed", False)}, clustering_done={st.session_state.get("clustering_done", False)}, analyses_done={st.session_state.get("analyses_done", False)}')
print(f'[DEBUG] File exists: {os.path.exists(DATA_PKL_PATH)}')

_badge_html = f"<span class='c3mod-badge c3mod-pulse'>STEP {_completed}/3</span>" if _completed < 3 else "<span class='c3mod-badge'>ALL STEPS COMPLETE</span>"
st.markdown(f"""
<div class="c3mod-header c3-fade-in">
    <h1>🧬 C3MOD - CANCER CLUSTERING AND CHARACTERIZATION USING MULTIOMICS DATA</h1>
    <div class="c3mod-sub">
        <span>Interactive pipeline for subtype discovery & characterization</span>
        {_badge_html}
        <span class='c3mod-badge dim'>Unsupervised Clustering</span>
    </div>
</div>
""", unsafe_allow_html=True)

# Stepper removed per user request

# ----------------- Logging Helper -----------------
class TeeLogger(io.TextIOBase):
    """Duplicates writes to stdout and an internal buffer for UI display."""
    def __init__(self, original):
        self.original = original
        self.buffer = io.StringIO()
    def write(self, s):
        self.original.write(s)
        self.buffer.write(s)
    def flush(self):
        self.original.flush()
    def get_value(self):
        return self.buffer.getvalue()

if 'tee_logger' not in st.session_state:
    st.session_state['tee_logger'] = TeeLogger(sys.stdout)
    sys.stdout = st.session_state['tee_logger']

def get_logs():
    return st.session_state['tee_logger'].get_value()

# ----------------- Output Folder Cleanup -----------------
def clear_output_files():
    base = OUTPUT_BASE
    if not os.path.isdir(base):
        return 0
    removed = 0
    for root, _, files in os.walk(base):
        for fname in files:
            fpath = os.path.join(root, fname)
            try:
                os.remove(fpath)
                removed += 1
            except Exception as e:  # noqa: E722
                print(f"[WARN] Could not remove {fpath}: {e}")
    print(f"[INFO] Cleared {removed} existing output files.")
    return removed

SESSION_TIMEOUT_SECONDS = int(os.environ.get('C3MOD_SESSION_TIMEOUT_SECONDS', '1800'))  # 30 min default
SESSION_MARKER = os.path.join(OUTPUT_BASE, '.session_marker')

def _session_expired() -> bool:
    if not os.path.exists(SESSION_MARKER):
        return False
    try:
        mtime = os.path.getmtime(SESSION_MARKER)
        return (time.time() - mtime) > SESSION_TIMEOUT_SECONDS
    except Exception:
        return False

def _touch_session_marker():
    os.makedirs(OUTPUT_BASE, exist_ok=True)
    try:
        with open(SESSION_MARKER, 'w') as f:
            f.write(str(time.time()))
    except Exception as e:
        print(f"[WARN] Could not write session marker: {e}")

def initialize_session():
    if _session_expired():
        print('[INFO] Previous session expired. Cleaning outputs...')
        clear_output_files()
        # Clear data.pkl on session expiry
        if os.path.exists(DATA_PKL_PATH):
            try:
                os.remove(DATA_PKL_PATH)
                print('[INFO] Removed expired data.pkl file.')
            except Exception as e:
                print(f'[WARN] Could not remove data.pkl: {e}')
        # Reset session state
        for key in ['preprocessed','clustering_done','analyses_done','pca_plots','classification_files','clustering_metrics',
                    'show_preprocessing_success','show_clustering_success','show_analyses_success']:
            if key in st.session_state:
                del st.session_state[key]
    if 'session_id' not in st.session_state:
        st.session_state['session_id'] = str(uuid.uuid4())
    _touch_session_marker()

initialize_session()

# Clean start: if this is a fresh app load, ensure clean state
if 'app_initialized' not in st.session_state:
    # Remove existing data.pkl to force fresh start
    if os.path.exists(DATA_PKL_PATH):
        try:
            os.remove(DATA_PKL_PATH)
            print('[INFO] Removed existing data.pkl for fresh start.')
        except Exception as e:
            print(f'[WARN] Could not remove data.pkl: {e}')
    
    # Clear any existing output files
    clear_output_files()
    
    # Reset all pipeline states
    for key in ['preprocessed','clustering_done','analyses_done','pca_plots','classification_files','clustering_metrics',
                'show_preprocessing_success','show_clustering_success','show_analyses_success']:
        if key in st.session_state:
            del st.session_state[key]
    
    st.session_state['app_initialized'] = True
    print('[INFO] App initialized with clean state.')

# Sidebar status panel
with st.sidebar:
    st.header("⚙️ Pipeline Status")
    
    # Force refresh status by re-checking session state as primary source
    preprocessing_done = bool(st.session_state.get('preprocessing_completed', False))
    clustering_done = bool(st.session_state.get('clustering_done', False))
    analyses_done = bool(st.session_state.get('analyses_done', False))
    
    # Hard refresh status display with detailed logging
    if preprocessing_done:
        print('[STATUS] Preprocessing: Session state set - showing Complete')
    else:
        print('[STATUS] Preprocessing: No session state - showing Pending')
        
    if clustering_done:
        print('[STATUS] Clustering: Session state set - showing Complete')
    else:
        print('[STATUS] Clustering: No session state - showing Pending')
        
    if analyses_done:
        print('[STATUS] Analyses: Session state set - showing Complete')
    else:
        print('[STATUS] Analyses: No session state - showing Pending')
    
    # Beautiful status cards
    st.markdown(f"""
    <div class="status-card {'complete' if preprocessing_done else 'pending'}">
        <div class="status-label">📊 Data Preprocessing</div>
        <div class="status-value {'complete' if preprocessing_done else 'pending'}">
            <span class="status-icon">{'✅' if preprocessing_done else '⏳'}</span>
            {'Complete' if preprocessing_done else 'Pending'}
        </div>
    </div>
    
    <div class="status-card {'complete' if clustering_done else 'pending'}">
        <div class="status-label">🎯 Clustering Analysis</div>
        <div class="status-value {'complete' if clustering_done else 'pending'}">
            <span class="status-icon">{'✅' if clustering_done else '⏳'}</span>
            {'Complete' if clustering_done else 'Pending'}
        </div>
    </div>
    
    <div class="status-card {'complete' if analyses_done else 'pending'}">
        <div class="status-label">🧬 Downstream Analysis</div>
        <div class="status-value {'complete' if analyses_done else 'pending'}">
            <span class="status-icon">{'✅' if analyses_done else '⏳'}</span>
            {'Complete' if analyses_done else 'Pending'}
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    st.divider()
    
    # Enhanced Refresh & Clear button
    if st.button("🔄 Refresh & Clear Analysis", key="clear_btn", type="secondary"):
        cleared = clear_output_files()
        # Remove data.pkl to reset preprocessing
        if os.path.exists(DATA_PKL_PATH):
            try:
                os.remove(DATA_PKL_PATH)
                print('[INFO] Removed data.pkl file.')
            except Exception as e:
                print(f'[WARN] Could not remove data.pkl: {e}')
        # Reset relevant session state keys
        for key in ['preprocessed','clustering_done','analyses_done','pca_plots','classification_files','clustering_metrics']:
            if key in st.session_state:
                del st.session_state[key]
        # Optionally clear cache
        try:
            st.cache_data.clear()  # type: ignore
        except Exception:
            pass
        st.success(f"Cleared {cleared} files. Pipeline reset to initial state.")
        st.rerun()
    
    st.divider()
    
    # Developer Info Section
    with st.expander("👨‍💻 Developer Info", expanded=False):
        st.markdown("""
        <div style="text-align: center; margin-bottom: 1rem;">
        """, unsafe_allow_html=True)
        
        # Display developer image
        st.image("images/developer-image.png", width=100)
        
        st.markdown("""
            <h3 style="margin: 0; color: #2c3e50;">Rajul Shakywar</h3>
            <p style="margin: 0.5rem 0; color: #7f8c8d; font-style: italic;">Software Engineer | GenAI</p>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("**About C3MOD:**")
        st.markdown("""
        C3MOD is a Python tool that identifies cancer subtypes in multi-omics data using unsupervised clustering algorithms and enables a variety of analyses on the identified subtypes.
        """)
        
        st.markdown("**Tech Stack:**")
        st.markdown("""
        • Unsupervised clustering  
        • Data analysis  
        • Machine Learning  
        • Python, R  
        • Bioinformatics  
        • Statistics  
        """)
        
        st.markdown("**Connect with me:**")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown("[![GitHub](https://img.shields.io/badge/GitHub-181717?style=for-the-badge&logo=github&logoColor=white)](https://github.com/rajulshakya777/C3MOD)")
        with col2:
            st.markdown("[![LinkedIn](https://img.shields.io/badge/LinkedIn-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/rajul-shakywar)")
        with col3:
            st.markdown("[![Portfolio](https://img.shields.io/badge/Portfolio-FF5722?style=for-the-badge&logo=web&logoColor=white)](https://rajulshakya777.github.io/portfolio/)")
    
    st.divider()
    
    # Data Download Section
    with st.expander("📁 Download Data", expanded=False):
        st.markdown("**Available Cancer Types:**")
        
        data_path = "data/input_data/TCGA_data"
        
        # Get list of available cancer types
        cancer_types = []
        try:
            for item in os.listdir(data_path):
                if os.path.isdir(os.path.join(data_path, item)) and not item.startswith('.'):
                    cancer_types.append(item)
            cancer_types.sort()
        except Exception as e:
            st.error(f"Error accessing data folder: {e}")
            cancer_types = []
        
        if cancer_types:
            selected_cancer = st.selectbox("Select Cancer Type for Download:", cancer_types, key="download_cancer")
            
            # Show files for selected cancer type
            if selected_cancer:
                cancer_data_path = os.path.join(data_path, selected_cancer)
                try:
                    files = []
                    for item in os.listdir(cancer_data_path):
                        file_path = os.path.join(cancer_data_path, item)
                        if os.path.isfile(file_path) and not item.startswith('.'):
                            files.append((item, file_path))
                    
                    if files:
                        st.markdown(f"**Files in {selected_cancer}:**")
                        for file_name, file_path in files:
                            try:
                                with open(file_path, 'rb') as f:
                                    file_data = f.read()
                                st.download_button(
                                    label=f"⬇️ {file_name}",
                                    data=file_data,
                                    file_name=file_name,
                                    mime="application/octet-stream",
                                    key=f"download_{selected_cancer}_{file_name}"
                                )
                            except Exception as e:
                                st.error(f"Error reading {file_name}: {e}")
                    else:
                        st.info(f"No files found in {selected_cancer} folder.")
                        
                except Exception as e:
                    st.error(f"Error accessing {selected_cancer} folder: {e}")
        else:
            st.info("No cancer type data available for download.")
    
    # (Session ID & log info intentionally hidden per user request)

tab1, tab2, tab3, tab4 = st.tabs([
    "🔧 Preprocess Data",
    "🎯 Cluster Analysis", 
    "🧬 Downstream Analysis",
    "📊 Results & Insights"
])  # Enhanced tab labels with icons

with tab1:
    st.markdown("""
    <div class="c3-fade-in">
        <div style="background:linear-gradient(135deg,#f0f9ff,#e0f2fe); padding:1.5rem; border-radius:16px; border:1px solid #bae6fd; margin-bottom:1.5rem;">
            <h2 style="color:#0c4a6e; margin:0 0 0.5rem 0; font-size:1.6rem; font-weight:700; font-family:'Poppins',sans-serif; letter-spacing:-0.02em;">🔧 Data Preprocessing</h2>
            <p style="color:#075985; margin:0; font-size:1.05rem; line-height:1.6; font-family:'Source Sans 3',sans-serif; font-weight:400;">Select a cancer type and start preprocessing. This step integrates and scales omics inputs and computes PCA for downstream analysis.</p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    cancer_type = st.selectbox("Select Cancer Type", [
        'ACC','BRCA','BLCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','PAAD','PCPG','PRAD','READ','SARC','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'
    ], index=16)  # default LUAD
    if st.button("Run Preprocessing", key="preprocess_btn"):
        with st.spinner("Preprocessing data..."):
            run_preprocessing(cancer_type)
            print('[INFO] Hard refresh: Updating preprocessing status to Complete')
            st.rerun()
    
    # Show success message if preprocessing just completed
    if st.session_state.get('show_preprocessing_success', False):
        st.markdown("""
        <div class="success-message">
            <div style="display: flex; align-items: flex-start;">
                <span class="icon">✅</span>
                <div>
                    <div class="title">Preprocessing Complete!</div>
                    <div class="description">Omics data has been successfully integrated, scaled, and prepared for clustering analysis.</div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        # Clear the success flag after showing it once
        st.session_state['show_preprocessing_success'] = False

preprocessed = bool(os.path.exists(DATA_PKL_PATH)) or bool(st.session_state.get('preprocessing_completed', False))

with tab2:
    st.markdown("""
    <div class="c3-fade-in">
        <div style="background:linear-gradient(135deg,#f0fdf4,#ecfdf5); padding:1.5rem; border-radius:16px; border:1px solid #bbf7d0; margin-bottom:1.5rem;">
            <h2 style="color:#14532d; margin:0 0 0.5rem 0; font-size:1.6rem; font-weight:700; font-family:'Poppins',sans-serif; letter-spacing:-0.02em;">🎯 Clustering Analysis</h2>
            <p style="color:#166534; margin:0; font-size:1.05rem; line-height:1.6; font-family:'Source Sans 3',sans-serif; font-weight:400;">Configure clustering parameters and run multiple algorithms to discover molecular subtypes in your data.</p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    if not preprocessed:
        st.info("Run preprocessing first.")
    else:
        top_cols = st.columns([1,1,2])
        with top_cols[0]:
            k = st.slider("Number of Clusters (K)", 2, 7, 3)
        with top_cols[1]:
            use_all = st.toggle("Select All Algos", value=False)
        # Algorithm badge selector
        if 'algo_selection' not in st.session_state:
            st.session_state['algo_selection'] = ['KMeans']
        if use_all:
            st.session_state['algo_selection'] = ALGORITHMS.copy()
        badge_cols = st.columns(len(ALGORITHMS))
        for algo, col in zip(ALGORITHMS, badge_cols):
            sel = algo in st.session_state['algo_selection']
            style = "algo-badge algo-selected mini-card" if sel else "algo-badge mini-card"
            if col.button(algo + (" ✓" if sel else ""), key=f"badge_{algo}"):
                if sel:
                    st.session_state['algo_selection'] = [a for a in st.session_state['algo_selection'] if a != algo]
                else:
                    st.session_state['algo_selection'].append(algo)
        st.caption("Click badges to toggle algorithms. SNF requires an R environment. ✓ = selected.")
        selected_algos = st.session_state['algo_selection'] if not use_all else ALGORITHMS
        snf_params = {}
        if 'SNF' in selected_algos:
            if not R_OK:
                st.warning("To run SNF analysis, please use CLI tool - [GitHub](https://github.com/rajulshakya777/C3MOD)")
                selected_algos = [a for a in selected_algos if a != 'SNF']
            else:
                with st.expander("SNF Parameters", expanded=False):
                    c1, c2, c3 = st.columns(3)
                    snf_params['n_neighbour'] = c1.number_input("Neighbors", 5, 25, 20)
                    snf_params['alpha'] = c2.slider("Alpha", 0.1, 0.9, 0.5, 0.1)
                    snf_params['T'] = c3.number_input("Iterations T", 5, 20, 15)
        
        # Ensure at least one algorithm is selected
        if not selected_algos:
            st.error("Please select at least one clustering algorithm.")
            st.stop()
            
        if st.button("🚀 Run Clustering", key="cluster_btn", type="primary"):
            with st.spinner("Running clustering algorithms..."):
                # Clear previous outputs prior to new clustering run
                clear_output_files()
                result_files, pca_plots, metrics = run_clustering(cancer_type, k, selected_algos, snf_params)
            st.session_state['clustering_done'] = True
            st.session_state['show_clustering_success'] = True
            st.session_state['pca_plots'] = pca_plots
            st.session_state['classification_files'] = result_files
            st.session_state['clustering_metrics'] = metrics
            print('[INFO] Hard refresh: Updating clustering status to Complete')
            st.rerun()
        
        # Show success message if clustering just completed
        if st.session_state.get('show_clustering_success', False):
            st.markdown("""
            <div class="success-message">
                <div style="display: flex; align-items: flex-start;">
                    <span class="icon">🎯</span>
                    <div>
                        <div class="title">Clustering Analysis Complete!</div>
                        <div class="description">Molecular subtypes have been successfully identified using selected algorithms. Review PCA plots and quality metrics below.</div>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            # Clear the success flag after showing it once
            st.session_state['show_clustering_success'] = False
        if st.session_state.get('pca_plots'):
            st.markdown("### PCA Plots")
            cols = st.columns(len(st.session_state['pca_plots']))
            for (algo, path), col in zip(st.session_state['pca_plots'].items(), cols):
                with col:
                    st.image(path, caption=algo, use_container_width=True)
        if st.session_state.get('clustering_metrics'):
            st.markdown("### Quality Metrics Dashboard")
            metrics_dict = st.session_state['clustering_metrics']
            all_keys = sorted({k for d in metrics_dict.values() for k in d.keys()})
            if all_keys:
                for algo, vals in metrics_dict.items():
                    st.markdown(f"#### {algo}")
                    boxes = []
                    for mk in all_keys:
                        val = vals.get(mk)
                        if val is None:
                            continue
                        cls = 'ok'
                        if mk == 'silhouette':
                            cls = 'good' if val >= 0.5 else ('ok' if val >= 0.25 else 'bad')
                        elif mk == 'davies_bouldin':
                            cls = 'good' if val < 0.8 else ('ok' if val < 1.5 else 'bad')
                        elif mk == 'calinski_harabasz':
                            cls = 'good' if val > 500 else ('ok' if val > 200 else 'bad')
                        boxes.append(f"<div class='metric-box {cls}'><h5>{mk.replace('_',' ')}</h5><div class='metric-val'>{val:.3f}</div></div>")
                    if boxes:
                        st.markdown(f"<div class='metrics-grid'>{''.join(boxes)}</div>", unsafe_allow_html=True)
                import pandas as _pd
                rows = []
                for algo, vals in metrics_dict.items():
                    row = {'algorithm': algo}
                    row.update(vals)
                    rows.append(row)
                dfm = _pd.DataFrame(rows)
                st.download_button("⬇️ Download Metrics CSV", data=dfm.to_csv(index=False).encode(), file_name="clustering_metrics.csv", mime="text/csv")
            else:
                st.info("Metrics not computed (insufficient clusters or data).")

clustering_done = st.session_state.get('clustering_done', False)

with tab3:
    st.markdown("""
    <div class="c3-fade-in">
        <div style="background:linear-gradient(135deg,#fefce8,#fef3c7); padding:1.5rem; border-radius:16px; border:1px solid #fde68a; margin-bottom:1.5rem;">
            <h2 style="color:#92400e; margin:0 0 0.5rem 0; font-size:1.6rem; font-weight:700; font-family:'Poppins',sans-serif; letter-spacing:-0.02em;">🧬 Downstream Analysis</h2>
            <p style="color:#a16207; margin:0; font-size:1.05rem; line-height:1.6; font-family:'Source Sans 3',sans-serif; font-weight:400;">Characterize discovered subtypes through survival, mutation, stage, and immune analysis to understand biological differences.</p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    if not clustering_done:
        st.info("Run clustering first.")
    else:
        analyses_to_run = st.multiselect("Analyses", list(ANALYSES.keys()))
        run_all = st.checkbox("Run All", value=False)
        if run_all:
            analyses_to_run = list(ANALYSES.keys())
        if st.button("🧪 Run Analyses", key="analyses_btn"):
            with st.spinner("Running selected analyses..."):
                data_tuple = load_pickle()
                if data_tuple and len(data_tuple) >= 5:
                    _, _, _, _, algo_choice_number = data_tuple
                    for label in analyses_to_run:
                        run_analysis(ANALYSES[label], algo_choice_number)
            st.session_state['analyses_done'] = True
            st.session_state['show_analyses_success'] = True
            print('[INFO] Hard refresh: Updating analyses status to Complete')
            st.rerun()
        
        # Show success message if analyses just completed
        if st.session_state.get('show_analyses_success', False):
            st.markdown("""
            <div class="success-message">
                <div style="display: flex; align-items: flex-start;">
                    <span class="icon">🧬</span>
                    <div>
                        <div class="title">Downstream Analysis Complete!</div>
                        <div class="description">Subtype characterization has been completed. Navigate to Results & Insights to explore findings and download outputs.</div>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            # Clear the success flag after showing it once
            st.session_state['show_analyses_success'] = False

with tab4:
    st.markdown("""
    <div class="c3-fade-in">
        <div style="background:linear-gradient(135deg,#f8fafc,#f1f5f9); padding:1.5rem; border-radius:16px; border:1px solid #cbd5e1; margin-bottom:1.5rem;">
            <h2 style="color:#334155; margin:0 0 0.5rem 0; font-size:1.6rem; font-weight:700; font-family:'Poppins',sans-serif; letter-spacing:-0.02em;">📊 Results & Insights</h2>
            <p style="color:#475569; margin:0; font-size:1.05rem; line-height:1.6; font-family:'Source Sans 3',sans-serif; font-weight:400;">Explore generated outputs, download results, and review analysis logs. Access visualizations and data tables from all pipeline steps.</p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    colA, colB = st.columns([2,1])
    with colA:
        section = st.selectbox("Result Section", ["clustering_results", "survival_analysis", "mutation_analysis", "stage_analysis", "immune_analysis"])    
        files = list_generated_files(section)
        if not files:
            st.info("No files yet for this section.")
        else:
            st.markdown(f"**{len(files)} files** in `{section}`")
            # Zip download
            with tempfile.TemporaryDirectory() as tmpd:
                zip_path = os.path.join(tmpd, f"{section}_results_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.zip")
                with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
                    for fpath in files:
                        zf.write(fpath, arcname=os.path.relpath(fpath, os.path.join(OUTPUT_BASE, section)))
                with open(zip_path, 'rb') as zf_data:
                    st.download_button("⬇️ Download All as ZIP", data=zf_data.read(), file_name=os.path.basename(zip_path), mime='application/zip')
            for fpath in files:
                rel = os.path.relpath(fpath)
                st.markdown(f"**{os.path.basename(rel)}**")
                if any(rel.lower().endswith(ext) for ext in ['.png', '.jpg', '.jpeg']):
                    st.image(rel, use_container_width=True)
                if rel.lower().endswith(('.txt', '.tsv', '.csv')):
                    try:
                        if rel.endswith('.tsv') or rel.endswith('.txt'):
                            df = pd.read_csv(rel, sep='\t')
                        else:
                            df = pd.read_csv(rel)
                        st.dataframe(df.head(250))
                    except Exception as e:
                        st.caption(f"Could not display table: {e}")
                with open(rel, 'rb') as f:
                    st.download_button(label="Download", data=f, file_name=os.path.basename(rel), key=rel)
                st.divider()
    with colB:
        st.markdown("### Live Logs")
        logs = get_logs()
        if st.button("🔄 Refresh Logs"):
            logs = get_logs()  # Re-pull
        st.markdown(f"<div class='log-box'>{logs if logs else 'No logs yet.'}</div>", unsafe_allow_html=True)

st.markdown("---")
## Footer removed per user request
