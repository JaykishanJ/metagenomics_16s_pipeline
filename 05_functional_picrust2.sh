#!/usr/bin/env bash
set -eo pipefail

# ==============================================================================
# 05_functional_picrust2.sh
# End-to-end PICRUSt2 functional annotation pipeline:
# 1. Clusters amplicons into representative OTUs (rep_seqs.fasta) & otu_table.tsv
# 2. Runs PICRUSt2 pipeline (EC numbers, KEGG KO, and MetaCyc pathways)
# 3. Generates publication-ready functional heatmaps and barplots (300 DPI PNG & PDF)
# ==============================================================================

INPUT_DIR="data/fastq"
OUTPUT_DIR="results/functional"
FIGURES_DIR="results/figures"
THREADS=8

usage() {
    echo "Usage: $0 [-i fastq_dir] [-o output_dir] [-f figures_dir] [-t threads]"
    echo "  -i  Directory containing raw FastQ files (Default: data/fastq)"
    echo "  -o  Directory for functional output tables (Default: results/functional)"
    echo "  -f  Directory for publication figures (Default: results/figures)"
    echo "  -t  Number of threads (Default: 8)"
    exit 1
}

while getopts "i:o:f:t:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        f) FIGURES_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FIGURES_DIR"

eval "$(/home/genome/miniconda3/bin/conda shell.bash hook 2>/dev/null || conda shell.bash hook 2>/dev/null)"

echo "=========================================================="
echo " Starting PICRUSt2 Metagenomic Functional Profiling"
echo " FastQ Directory : $INPUT_DIR"
echo " Functional Dir  : $OUTPUT_DIR"
echo " Figures Dir     : $FIGURES_DIR"
echo " Threads         : $THREADS"
echo "=========================================================="

# 1. Generate representative sequences and abundance table
conda activate metagenomics_16s
echo "[+] Step 1: Generating representative sequences and OTU table..."
python3 "$SCRIPT_DIR/scripts/prepare_picrust2.py" \
    -i "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    -t "$THREADS"

REP_SEQS="$OUTPUT_DIR/rep_seqs.fasta"
OTU_TABLE="$OUTPUT_DIR/otu_table.tsv"
PICRUST_OUT="$OUTPUT_DIR/picrust2_out"

if [ ! -f "$REP_SEQS" ] || [ ! -f "$OTU_TABLE" ]; then
    echo "[-] Error: Failed to generate PICRUSt2 inputs."
    exit 1
fi

# 2. Run PICRUSt2 Pipeline
echo "[+] Step 2: Running PICRUSt2 prediction pipeline (EC, KO, MetaCyc)..."
conda activate qiime2-metagenome-2024.10

rm -rf "$PICRUST_OUT"
picrust2_pipeline.py \
    -s "$REP_SEQS" \
    -i "$OTU_TABLE" \
    -o "$PICRUST_OUT" \
    --threads "$THREADS" \
    --verbose

PATHWAY_FILE="$PICRUST_OUT/pathways_out/path_abun_unstrat.tsv.gz"
if [ ! -f "$PATHWAY_FILE" ]; then
    PATHWAY_FILE="$PICRUST_OUT/pathways_out/path_abun_unstrat.tsv"
fi

if [ ! -f "$PATHWAY_FILE" ]; then
    echo "[-] Error: PICRUSt2 pathway output not found."
    exit 1
fi

# 3. Plot functional pathways
echo "[+] Step 3: Generating publication-grade functional figures..."
conda activate metagenomics_16s
python3 "$SCRIPT_DIR/scripts/plot_functional_pathways.py" \
    -i "$PATHWAY_FILE" \
    -f "$FIGURES_DIR" \
    -o "$OUTPUT_DIR" \
    --top-n 20

echo "=========================================================="
echo " PICRUSt2 functional annotation completed successfully!"
echo " Enzyme predictions    : $PICRUST_OUT/EC_metagenome_out/"
echo " KEGG KO predictions   : $PICRUST_OUT/KO_metagenome_out/"
echo " MetaCyc pathways      : $PICRUST_OUT/pathways_out/"
echo " Publication figures   : $FIGURES_DIR/functional_pathways_*"
echo "=========================================================="
