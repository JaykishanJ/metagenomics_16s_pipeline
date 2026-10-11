#!/usr/bin/env bash
set -eo pipefail

# ==============================================================================
# run_pipeline.sh - Master Turnkey Metagenomics 16S Pipeline
# Complete publication-ready pipeline executing end-to-end:
#   Stage 1: Environment & Database Verification
#   Stage 2: Kraken2 Taxonomic Profiling & Bracken Bayesian Abundance
#   Stage 3: Merged Species Abundance Matrix Construction
#   Stage 4: Interactive Krona Multi-Level Radial Charts (HTML)
#   Stage 5: Publication-Grade Diversity Analysis (Alpha & Beta PCoA, Heatmaps, Stacked Bars)
#   Stage 6: PICRUSt2 Functional Annotation (EC numbers, KEGG KO, MetaCyc Pathways)
# ==============================================================================

# Default settings
INPUT_DIR="data/fastq"
OUTPUT_DIR="results"
DB_DIR="kraken2_db"
METADATA_FILE=""
THREADS=8
JOBS=4
READ_LEN=150
CONFIDENCE=0.0
SKIP_FUNCTIONAL=false
SKIP_KRONA=false
SKIP_DIVERSITY=false

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required / Main Options:
  -i  Input directory containing paired-end FastQ files (Default: data/fastq)
  -o  Main output directory (Default: results)
  -d  Kraken2 database directory (Default: kraken2_db)
  -m  Optional sample metadata file (TSV with 'SampleID' and 'Group')
  -t  Total threads to allocate (Default: 8)
  -j  Number of parallel classification jobs (Default: 4)
  -r  Read length for Bracken estimation (Default: 150)
  -c  Kraken2 classification confidence threshold (Default: 0.0)

Modular Flags:
  --skip-functional   Skip PICRUSt2 functional pathway prediction
  --skip-krona        Skip Krona HTML visualization
  --skip-diversity    Skip Alpha/Beta diversity calculations and figures
  -h, --help          Display this help message and exit

Examples:
  # Standard turnkey execution:
  $0 -i /path/to/fastqs -o results -t 16

  # Execution with experimental group metadata:
  $0 -i data/fastq -m metadata.tsv -o results -t 8
EOF
    exit 1
}

# Parse options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT_DIR="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -d) DB_DIR="$2"; shift 2 ;;
        -m) METADATA_FILE="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -j) JOBS="$2"; shift 2 ;;
        -r) READ_LEN="$2"; shift 2 ;;
        -c) CONFIDENCE="$2"; shift 2 ;;
        --skip-functional) SKIP_FUNCTIONAL=true; shift ;;
        --skip-krona) SKIP_KRONA=true; shift ;;
        --skip-diversity) SKIP_DIVERSITY=true; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUTPUT_DIR"

echo "================================================================================"
echo "          TURNKEY PUBLICATION-READY METAGENOMICS 16S PIPELINE                   "
echo "================================================================================"
echo " Input FastQ Dir  : $INPUT_DIR"
echo " Output Dir       : $OUTPUT_DIR"
echo " Reference DB     : $DB_DIR"
echo " Metadata File    : ${METADATA_FILE:-'None provided (cohort analysis)'}"
echo " Threads / Jobs   : $THREADS threads | $JOBS parallel jobs"
echo " Bracken Read Len : $READ_LEN bp"
echo "================================================================================"

START_TIME=$(date +%s)

# ------------------------------------------------------------------------------
# Stage 1: Setup & Database Indexing Check
# ------------------------------------------------------------------------------
echo ""
echo ">>> [STAGE 1/6] Checking Environment & Reference Databases..."
if [ ! -f "$DB_DIR/hash.k2d" ] || [ ! -f "$DB_DIR/database${READ_LEN}mers.kmer_distrib" ]; then
    echo "[+] Reference DB or Bracken distribution missing. Running setup..."
    bash "$SCRIPT_DIR/01_setup.sh"
else
    echo "[+] Reference database verified in $DB_DIR."
fi

# ------------------------------------------------------------------------------
# Stage 2: Kraken2 Classification & Bracken Bayesian Abundance
# ------------------------------------------------------------------------------
echo ""
echo ">>> [STAGE 2/6] Running Kraken2 Classification & Bracken Bayesian Re-estimation..."
bash "$SCRIPT_DIR/02_classify_pipeline.sh" \
    -i "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    -d "$DB_DIR" \
    -t "$THREADS" \
    -j "$JOBS" \
    -r "$READ_LEN" \
    -c "$CONFIDENCE"

# ------------------------------------------------------------------------------
# Stage 3: Downstream Abundance Matrix Construction
# ------------------------------------------------------------------------------
echo ""
echo ">>> [STAGE 3/6] Aggregating Species Abundance Matrix..."
MATRIX_FILE="$OUTPUT_DIR/taxonomy/species_abundance_matrix.txt"
bash "$SCRIPT_DIR/03_downstream_analysis.sh" \
    -i "$OUTPUT_DIR" \
    -o "$MATRIX_FILE"

# ------------------------------------------------------------------------------
# Stage 4: Krona Interactive Radial Visualizations
# ------------------------------------------------------------------------------
if [ "$SKIP_KRONA" = false ]; then
    echo ""
    echo ">>> [STAGE 4/6] Generating Interactive Krona Radial Charts..."
    bash "$SCRIPT_DIR/03_krona_visualization.sh" \
        -i "$OUTPUT_DIR" \
        -o "$OUTPUT_DIR/krona"
else
    echo ""
    echo ">>> [STAGE 4/6] Skipping Krona visualization (--skip-krona enabled)."
fi

# ------------------------------------------------------------------------------
# Stage 5: Comprehensive Downstream Diversity & Statistical Analysis
# ------------------------------------------------------------------------------
if [ "$SKIP_DIVERSITY" = false ]; then
    echo ""
    echo ">>> [STAGE 5/6] Generating Publication Diversity Figures, QC Curves & Statistics..."
    
    # 5a. Sequencing Depth QC & Rarefaction Curves
    echo "[+] Running Sequencing Depth QC & In-Silico Rarefaction Analysis..."
    python3 "$SCRIPT_DIR/scripts/qc_depth_analysis.py" \
        -r "$OUTPUT_DIR" \
        -m "$MATRIX_FILE" \
        -o "$OUTPUT_DIR/qc" \
        -f "$OUTPUT_DIR/figures"
        
    # 5b. Multi-Rank Taxonomic Composition (Phylum, Genus & Group-Mean)
    echo "[+] Generating Multi-Rank Stacked Barplots & Group-Mean Comparisons..."
    MULTI_ARGS=(-r "$OUTPUT_DIR" -t "$OUTPUT_DIR/taxonomy" -f "$OUTPUT_DIR/figures")
    if [ -n "$METADATA_FILE" ] && [ -f "$METADATA_FILE" ]; then
        MULTI_ARGS+=(-m "$METADATA_FILE")
    fi
    python3 "$SCRIPT_DIR/scripts/plot_multilevel_taxa.py" "${MULTI_ARGS[@]}"

    # 5c. Alpha & Beta Diversity (PCoA Ordination & Species Heatmap)
    echo "[+] Computing Alpha & Beta Diversity (Bray-Curtis PCoA)..."
    DIV_ARGS=(-i "$MATRIX_FILE" -d "$OUTPUT_DIR/diversity" -f "$OUTPUT_DIR/figures")
    if [ -n "$METADATA_FILE" ] && [ -f "$METADATA_FILE" ]; then
        DIV_ARGS+=(-m "$METADATA_FILE")
    fi
    bash "$SCRIPT_DIR/04_diversity_analysis.sh" "${DIV_ARGS[@]}"

    # 5d. Differential Abundance Testing & PERMANOVA
    if [ -n "$METADATA_FILE" ] && [ -f "$METADATA_FILE" ]; then
        echo "[+] Performing Differential Abundance (Volcano, Forest) & PERMANOVA..."
        python3 "$SCRIPT_DIR/scripts/differential_abundance.py" \
            -i "$MATRIX_FILE" \
            -m "$METADATA_FILE" \
            -o "$OUTPUT_DIR/differential" \
            -f "$OUTPUT_DIR/figures"
    fi

    # 5e. Microbial Co-occurrence & Association Networks
    echo "[+] Inferring Microbial Co-occurrence & Correlation Networks..."
    python3 "$SCRIPT_DIR/scripts/network_analysis.py" \
        -i "$MATRIX_FILE" \
        -o "$OUTPUT_DIR/diversity" \
        -f "$OUTPUT_DIR/figures"
else
    echo ""
    echo ">>> [STAGE 5/6] Skipping diversity analysis (--skip-diversity enabled)."
fi

# ------------------------------------------------------------------------------
# Stage 6: Functional & Metabolic Pathway Inference
# ------------------------------------------------------------------------------
if [ "$SKIP_FUNCTIONAL" = false ]; then
    echo ""
    echo ">>> [STAGE 6/6] Inferring Microbial Metabolic & Functional Pathways..."
    python3 "$SCRIPT_DIR/scripts/functional_profiling.py" \
        -i "$MATRIX_FILE" \
        -o "$OUTPUT_DIR/functional" \
        -f "$OUTPUT_DIR/figures"
        
    echo "[+] Running PICRUSt2 amplicon sequence preparation..."
    bash "$SCRIPT_DIR/05_functional_picrust2.sh" \
        -i "$INPUT_DIR" \
        -o "$OUTPUT_DIR/functional" \
        -f "$OUTPUT_DIR/figures" \
        -t "$THREADS" || true
else
    echo ""
    echo ">>> [STAGE 6/6] Skipping functional profiling (--skip-functional enabled)."
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "================================================================================"
echo "          PIPELINE EXECUTION COMPLETED SUCCESSFULLY IN ${ELAPSED}s              "
echo "================================================================================"
echo " Summary of Key Deliverables:"
echo "  - Species Abundance Matrix: $MATRIX_FILE"
echo "  - Interactive Krona Chart : $OUTPUT_DIR/krona/all_samples_krona.html"
echo "  - Alpha Diversity Metrics : $OUTPUT_DIR/diversity/alpha_diversity_metrics.tsv"
echo "  - Beta Distance Matrix    : $OUTPUT_DIR/diversity/beta_bray_curtis_distance.tsv"
echo "  - High-Res Figures (300DPI): $OUTPUT_DIR/figures/ (*.png and *.pdf)"
if [ "$SKIP_FUNCTIONAL" = false ]; then
echo "  - Functional Pathways     : $OUTPUT_DIR/functional/top_functional_pathways.tsv"
fi
echo "================================================================================"
