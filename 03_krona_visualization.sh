#!/usr/bin/env bash
set -eo pipefail

# ==============================================================================
# 03_krona_visualization.sh
# Generate interactive Krona radial pie charts from Kraken2 reports.
# Produces both per-sample HTML charts and a unified multi-sample comparison chart.
# ==============================================================================

# Activate Conda environment
eval "$(/home/genome/miniconda3/bin/conda shell.bash hook 2>/dev/null || conda shell.bash hook 2>/dev/null)"
conda activate metagenomics_16s

INPUT_DIR="results"
OUTPUT_DIR="results/krona"

usage() {
    echo "Usage: $0 [-i input_dir] [-o output_dir]"
    echo "  -i  Directory containing Kraken2 .report files (Default: results)"
    echo "  -o  Output directory for Krona charts (Default: results/krona)"
    exit 1
}

while getopts "i:o:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUTPUT_DIR"

echo "=========================================================="
echo " Starting Krona Interactive Taxonomic Visualization"
echo " Input Directory  : $INPUT_DIR"
echo " Output Directory : $OUTPUT_DIR"
echo "=========================================================="

REPORT_FILES=("$INPUT_DIR"/*.report)

# Filter out _bracken_species.report so we use the full Kraken hierarchical tree
RAW_REPORTS=()
for rep in "${REPORT_FILES[@]}"; do
    if [[ -f "$rep" && ! "$rep" =~ _bracken_species\.report$ ]]; then
        RAW_REPORTS+=("$rep")
    fi
done

if [ ${#RAW_REPORTS[@]} -eq 0 ]; then
    echo "[-] Error: No valid .report files found in $INPUT_DIR"
    exit 1
fi

echo "[+] Found ${#RAW_REPORTS[@]} Kraken2 reports to process."

KRONA_INPUTS=()
for rep in "${RAW_REPORTS[@]}"; do
    sample_name=$(basename "$rep" .report)
    sample_txt="$OUTPUT_DIR/${sample_name}.krona.txt"
    sample_html="$OUTPUT_DIR/${sample_name}.html"
    
    echo "[+] Converting $sample_name to Krona text format..."
    python3 "$SCRIPT_DIR/scripts/kraken2krona.py" -i "$rep" -o "$sample_txt"
    
    echo "[+] Generating individual Krona HTML for $sample_name..."
    ktImportText "$sample_txt" -o "$sample_html"
    
    # Store with sample label for combined chart: path,SampleName
    KRONA_INPUTS+=("$sample_txt,$sample_name")
done

echo "[+] Generating unified multi-sample comparison Krona chart..."
ktImportText "${KRONA_INPUTS[@]}" -o "$OUTPUT_DIR/all_samples_krona.html"

echo "=========================================================="
echo " Krona visualization completed successfully!"
echo " Individual sample charts: $OUTPUT_DIR/<sample>.html"
echo " Multi-sample interactive chart: $OUTPUT_DIR/all_samples_krona.html"
echo "=========================================================="
