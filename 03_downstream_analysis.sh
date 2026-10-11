#!/bin/bash

# ==============================================================================
# Metagenomics 16S Pipeline: Downstream Analysis
# ==============================================================================
# This script processes the Bracken species-level output files to generate
# a summary table of taxonomic abundances for all samples.
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
INPUT_DIR=""
OUTPUT_FILE=""

# Parse arguments
while getopts "i:o:" opt; do
  case $opt in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_FILE="$OPTARG" ;;
    *) echo "Usage: $0 -i <input_dir> -o <output_file>"; exit 1 ;;
  esac
done

# Validation
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_FILE" ]]; then
    echo "Error: Missing required arguments."
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"

echo "Combining Bracken outputs from $INPUT_DIR..."

# Create a temporary directory for processing
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Combine files using a simple python script for robust merging
cat << 'EOF' > "$TEMP_DIR/merge.py"
import pandas as pd
import glob
import os
import sys

input_dir = sys.argv[1]
output_file = sys.argv[2]

all_data = []

# Iterate over all .bracken files
for file_path in glob.glob(os.path.join(input_dir, "*.bracken")):
    sample_name = os.path.basename(file_path).replace(".bracken", "")

    # Read bracken output (tab-separated)
    df = pd.read_csv(file_path, sep='\t')

    # Select relevant columns (name, taxonomy_id, fraction_total_reads)
    df = df[['name', 'fraction_total_reads']]
    df = df.rename(columns={'fraction_total_reads': sample_name})

    all_data.append(df)

# Merge all dataframes on 'name'
if all_data:
    final_df = all_data[0]
    for df in all_data[1:]:
        final_df = pd.merge(final_df, df, on='name', how='outer')

    final_df.fillna(0, inplace=True)
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"Summary table saved to {output_file}")
else:
    print("No .bracken files found.")
EOF

# Run the merge script
python3 "$TEMP_DIR/merge.py" "$INPUT_DIR" "$OUTPUT_FILE"

echo "Downstream analysis completed."
