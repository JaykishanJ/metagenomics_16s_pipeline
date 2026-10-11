#!/bin/bash

# ==============================================================================
# Metagenomics 16S Pipeline: Kraken2 + Bracken (Production Ready)
# ==============================================================================
# This script performs species-level classification on paired-end FASTQ files.
# It uses GNU Parallel to maximize CPU utilization across samples.
# ==============================================================================

set -euo pipefail

# Initialize conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate metagenomics_16s 2>/dev/null || true
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
    conda activate metagenomics_16s 2>/dev/null || true
fi

# Default values
THREADS_PER_SAMPLE=4
MAX_JOBS=2
CONFIDENCE="0.0"
READ_LEN=150
DB_PATH=""
INPUT_DIR=""
OUTPUT_DIR=""

# Parse arguments
while getopts "i:o:d:t:j:c:r:" opt; do
  case $opt in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    d) DB_PATH="$OPTARG" ;;
    t) THREADS_PER_SAMPLE="$OPTARG" ;;
    j) MAX_JOBS="$OPTARG" ;;
    c) CONFIDENCE="$OPTARG" ;;
    r) READ_LEN="$OPTARG" ;;
    *) echo "Usage: $0 -i <input_dir> -o <output_dir> -d <db_path> -t <threads> -j <jobs> -c <confidence> -r <read_len>"; exit 1 ;;
  esac
done

# Validation
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$DB_PATH" ]]; then
    echo "Error: Missing required arguments."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Function to process a single sample
process_sample() {
    local r1_file="$1"
    local output_dir="$2"
    local db_path="$3"
    local threads="$4"
    local confidence="$5"
    local read_len="$6"

    local sample=""
    local r2_file=""

    # Flexible sample naming handling
    if [[ "$r1_file" == *_R1.fastq.gz ]]; then
        sample=$(basename "$r1_file" _R1.fastq.gz)
        r2_file="${r1_file/_R1.fastq.gz/_R2.fastq.gz}"
    elif [[ "$r1_file" == *_1.fastq.gz ]]; then
        sample=$(basename "$r1_file" _1.fastq.gz)
        r2_file="${r1_file/_1.fastq.gz/_2.fastq.gz}"
    else
        echo "Error: Unexpected filename format: $r1_file"
        return 1
    fi

    echo "Processing sample: $sample"

    # Kraken2
    kraken2 --db "$db_path" \
            --threads "$threads" \
            --confidence "$confidence" \
            --paired "$r1_file" "$r2_file" \
            --report "$output_dir/${sample}.report" \
            --output "$output_dir/${sample}.kraken"

    # Bracken
    bracken -d "$db_path" \
            -i "$output_dir/${sample}.report" \
            -o "$output_dir/${sample}.bracken" \
            -r "$read_len" \
            -l S \
            -t "$threads"
}

export -f process_sample

# Collect all R1 files (supporting both *_R1.fastq.gz and *_1.fastq.gz)
shopt -s nullglob
r1_files=("$INPUT_DIR"/*_R1.fastq.gz "$INPUT_DIR"/*_1.fastq.gz)
shopt -u nullglob

if [ ${#r1_files[@]} -eq 0 ]; then
    echo "Error: No matching FASTQ files (*_R1.fastq.gz or *_1.fastq.gz) found in $INPUT_DIR"
    exit 1
fi

echo "Found ${#r1_files[@]} sample(s) to process."

# Use GNU Parallel for parallel sample processing if available
if command -v parallel &> /dev/null; then
    echo "Using GNU Parallel to process samples ($MAX_JOBS concurrent jobs)..."
    parallel -j "$MAX_JOBS" process_sample {} "$OUTPUT_DIR" "$DB_PATH" "$THREADS_PER_SAMPLE" "$CONFIDENCE" "$READ_LEN" ::: "${r1_files[@]}"
else
    echo "GNU Parallel not found. Falling back to sequential processing."
    for r1 in "${r1_files[@]}"; do
        process_sample "$r1" "$OUTPUT_DIR" "$DB_PATH" "$THREADS_PER_SAMPLE" "$CONFIDENCE" "$READ_LEN"
    done
fi

echo "Pipeline completed successfully."
