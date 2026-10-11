#!/bin/bash

# ==============================================================================
# Setup: Metagenomics 16S Pipeline Environment
# ==============================================================================
# This script sets up a Conda/Mamba environment and downloads the Kraken2
# database for 16S rRNA classification.
# ==============================================================================

set -euo pipefail

# Initialize conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

ENV_NAME="metagenomics_16s"
DB_DIR="./kraken2_db"

echo "Creating Conda environment: $ENV_NAME"
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found, using Conda."
    conda create -n "$ENV_NAME" -y
else
    mamba create -n "$ENV_NAME" -y
fi

# Activate environment (or use 'conda run -n')
conda activate "$ENV_NAME"

echo "Installing pipeline dependencies into $ENV_NAME..."
if ! command -v mamba &> /dev/null; then
    conda install -n "$ENV_NAME" -c bioconda -c conda-forge kraken2 bracken krona vsearch pandas scipy matplotlib seaborn scikit-learn -y
else
    mamba install -n "$ENV_NAME" -c bioconda -c conda-forge kraken2 bracken krona vsearch pandas scipy matplotlib seaborn scikit-learn -y
fi

# Optional PICRUSt2 functional profiling environment
echo "Checking PICRUSt2 environment..."
if ! conda info --envs | grep -q "qiime2-metagenome"; then
    echo "[+] Creating dedicated picrust2 environment..."
    conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.0.3_b -y || true
fi

# Setup Database
mkdir -p "$DB_DIR"
echo "Downloading and building Kraken2 16S database (greengenes)..."
# Download and build greengenes standard 16S library
kraken2-build --special greengenes --db "$DB_DIR"

echo "Building Bracken database..."
# Required step for Bracken to work later. Assumes 150bp read length as default.
bracken-build -d "$DB_DIR" -t 8 -l 150

echo "Setup complete. Database is in $DB_DIR"
