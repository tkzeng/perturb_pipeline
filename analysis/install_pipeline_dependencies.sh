#!/bin/bash
# Installation script for Perturb-seq Pipeline Dependencies
# This script creates a conda environment with all required dependencies

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ENV_NAME="perturb_pipeline"

echo "=================================================="
echo "Perturb-seq Pipeline Installation Script"
echo "=================================================="
echo ""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --env-name NAME    Name for conda environment (default: perturb_pipeline)"
            echo "  --help            Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "Creating conda environment: $ENV_NAME"
echo ""

# Check if mamba or conda is available
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba for faster installation"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo "Using conda (consider installing mamba for faster installation)"
else
    echo "Error: Neither mamba nor conda found. Please install Anaconda/Miniconda/Mambaforge first."
    exit 1
fi

# Check if environment already exists
#if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
#    echo "Warning: Environment '$ENV_NAME' already exists."
#    read -p "Do you want to remove and recreate it? (y/N): " -n 1 -r
#    echo
#    if [[ $REPLY =~ ^[Yy]$ ]]; then
#        echo "Removing existing environment..."
#        $CONDA_CMD env remove -n $ENV_NAME -y
#    else
#        echo "Exiting without changes."
#        exit 0
#    fi
#fi

# Create environment with Python 3.12
echo "Creating base environment with Python 3.12..."
#$CONDA_CMD create -n $ENV_NAME python=3.12 -y

# Activate environment
echo "Activating environment..."
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# Add conda channels
echo "Adding conda channels..."
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict

# Install kallisto-bustools and related tools
echo ""
echo "Installing kallisto-bustools and alignment tools..."
#$CONDA_CMD install -y \
#    kallisto \
#    bustools \
#    samtools \
#    pigz \
#    seqtk

# Install custom kb-python with remove_collisions function for kite indexing
echo ""
echo "Installing custom kb-python with collision detection..."
echo "This version includes the remove_collisions function required for kite indexing"
#pip install git+https://github.com/tkzeng/kb_python.git

# Install R and R packages
echo ""
echo "Installing R base and dependencies..."
$CONDA_CMD install -y \
    r-base=4.4 \
    r-argparse \
    r-matrix || { echo "Failed to install R packages"; exit 1; }

# Note: We use a standalone barcodeRanks implementation that doesn't require full DropletUtils
echo ""
echo "Using minimal barcodeRanks implementation (no DropletUtils needed)"
echo "This uses just the knee/inflection detection algorithm without full DropletUtils dependencies"

# Install Python data science packages
echo ""
echo "Installing Python data science packages..."
$CONDA_CMD install -y \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    scanpy \
    anndata \
    pysam \
    pyyaml \
    openpyxl \
    igraph \
    leidenalg

# Install additional Python packages
echo ""
echo "Installing additional Python packages..."
$CONDA_CMD install -y \
    streamlit \
    gffutils

# Install Snakemake
echo ""
echo "Installing Snakemake..."
$CONDA_CMD install -y snakemake

# Install Snakemake Slurm executor plugin for cluster execution
echo ""
echo "Installing Snakemake Slurm executor plugin..."
pip install snakemake-executor-plugin-slurm

# Create environment.yml for reproducibility
echo ""
echo "Exporting environment to environment.yml..."
$CONDA_CMD env export --no-builds > "${SCRIPT_DIR}/environment.yml"

echo ""
echo "=================================================="
echo "Installation completed!"
echo "=================================================="
echo ""
echo "To activate the environment, run:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To run the pipeline:"
echo "  snakemake -s preprocessing.smk --configfile config.yaml"
echo "  snakemake -s downstream.smk --configfile config.yaml"
echo ""
