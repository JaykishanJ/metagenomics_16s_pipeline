#!/usr/bin/env bash
set -eo pipefail

# ==============================================================================
# 04_diversity_analysis.sh
# Performs Alpha Diversity, Beta Diversity PCoA ordination, compositional stacked
# barplots, and clustered heatmaps. Exports TSVs, 300 DPI PNGs, and vector PDFs.
# ==============================================================================

eval "$(/home/genome/miniconda3/bin/conda shell.bash hook 2>/dev/null || conda shell.bash hook 2>/dev/null)"
conda activate metagenomics_16s

MATRIX_FILE="results/taxonomy/species_abundance_matrix.txt"
METADATA_FILE=""
DIVERSITY_DIR="results/diversity"
FIGURES_DIR="results/figures"

usage() {
    echo "Usage: $0 [-i matrix_file] [-m metadata_tsv] [-d diversity_dir] [-f figures_dir]"
    echo "  -i  Path to species abundance matrix (Default: results/taxonomy/species_abundance_matrix.txt)"
    echo "  -m  Path to optional sample metadata TSV (e.g. SampleID, Group)"
    echo "  -d  Directory for diversity metric tables (Default: results/diversity)"
    echo "  -f  Directory for publication figures (Default: results/figures)"
    exit 1
}

while getopts "i:m:d:f:h" opt; do
    case "$opt" in
        i) MATRIX_FILE="$OPTARG" ;;
        m) METADATA_FILE="$OPTARG" ;;
        d) DIVERSITY_DIR="$OPTARG" ;;
        f) FIGURES_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ ! -f "$MATRIX_FILE" ]; then
    echo "[-] Error: Abundance matrix $MATRIX_FILE not found."
    exit 1
fi

CMD=(python3 "$SCRIPT_DIR/scripts/diversity_analysis.py" -i "$MATRIX_FILE" -d "$DIVERSITY_DIR" -f "$FIGURES_DIR")

if [ -n "$METADATA_FILE" ] && [ -f "$METADATA_FILE" ]; then
    CMD+=(-m "$METADATA_FILE")
fi

"${CMD[@]}"
