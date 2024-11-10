import subprocess
import sys

# ANSI color codes for the log
BOLD = "\033[1m"
RESET = "\033[0m"
CYAN = "\033[1;36m"
GREEN = "\033[1;32m"
YELLOW = "\033[1;33m"
RED = "\033[1;31m"

# Emojis for better visualization
PACKAGE_EMOJI = "📦"
ERROR_EMOJI = "❌"
CHECK_EMOJI = f"{GREEN}✔{RESET}"  # Green check mark for success

# Function to install a Python package if not already installed
def install_package(package_name, pip_name=None):
    if pip_name is None:
        pip_name = package_name
    try:
        __import__(package_name)
    except ImportError:
        print(f"{PACKAGE_EMOJI} {RED}{package_name} not found. Installing {pip_name}...{RESET}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])
        print(f"{PACKAGE_EMOJI} {CHECK_EMOJI} {GREEN}{package_name} installed successfully!{RESET}")
    else:
        print(f" {CHECK_EMOJI} {YELLOW}{package_name} is already installed.{RESET}")

# List of required Python packages
required_packages = {
    "pandas": "pandas",
    "numpy": "numpy",
    "sklearn": "scikit-learn",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "scipy": "scipy",
    "skfuzzy": "scikit-fuzzy",
    "rpy2": "rpy2",
    "pickle": "pickle5",  # pickle5 is the package name for pip
    "os": None,  # 'os' is part of the Python standard library, no pip installation required
}

# Install each Python package if not present
for package_name, pip_name in required_packages.items():
    if pip_name is not None:
        install_package(package_name, pip_name)

# Import the installed Python packages
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering, DBSCAN, MeanShift
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, ClusterMixin
import skfuzzy as fuzz
import os
import pickle
import rpy2.robjects as robjects

# Function to install R packages using rpy2
def install_r_package(r_package_name):
    try:
        # Import the necessary rpy2 functions
        from rpy2.robjects.packages import importr
        utils = importr('utils')
        
        # Check if the R package is installed
        utils.chooseCRANmirror(ind=1)
        installed_packages = utils.installed_packages()
        
        if r_package_name not in installed_packages:
            print(f"{PACKAGE_EMOJI} {RED}{r_package_name} not found. Installing {r_package_name}...{RESET}")
            utils.install_packages(r_package_name)
            print(f"{PACKAGE_EMOJI} {CHECK_EMOJI} {GREEN}{r_package_name} installed successfully!{RESET}")
        else:
            print(f" {CHECK_EMOJI} {YELLOW}{r_package_name} is already installed.{RESET}")
    
    except Exception as e:
        print(f"{ERROR_EMOJI} {RED}Failed to install {r_package_name}: {e}{RESET}")

# List of required R packages
r_packages = [
    "SNFtool",
    "data.table",
    "dplyr",
    "stringr",
    "factoextra",
    "cluster",
    "pryr",
    "bigmemory",
    "data.table",
    "survival",
    "survminer",
    "ggfortify",
    "ggplot2",
    "ggtext",
    "cowplot"
]

# Install each R package if not present
for r_package in r_packages:
    install_r_package(r_package)

# Final confirmation message
# print(f"{CHECK_EMOJI} {CYAN}All Python and R packages are installed and imported successfully!{RESET}")
