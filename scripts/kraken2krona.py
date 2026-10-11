#!/usr/bin/env python3
"""
kraken2krona.py - Convert Kraken2 / Bracken reports to Krona hierarchical text format.
Parses the indentation tree in Kraken report format (col 3: direct read count, col 6+: taxon name)
and emits: <count>\t<Domain>\t<Phylum>\t...\t<Species>
"""

import sys
import os
import argparse

def convert_report(report_file, output_file):
    lineage = {}
    total_written = 0
    with open(report_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            try:
                direct_count = int(parts[2])
            except ValueError:
                continue
            
            raw_name = parts[5]
            depth = (len(raw_name) - len(raw_name.lstrip(" "))) // 2
            name = raw_name.strip()
            
            if name.lower() in ("unclassified", "root"):
                continue
            
            lineage[depth] = name
            # Remove any deeper levels from previous branches
            for k in list(lineage.keys()):
                if k > depth:
                    del lineage[k]
                    
            if direct_count > 0:
                current_lineage = [lineage[i] for i in sorted(lineage.keys())]
                outfile.write(f"{direct_count}\t" + "\t".join(current_lineage) + "\n")
                total_written += 1
                
    return total_written

def main():
    parser = argparse.ArgumentParser(description="Convert Kraken2 reports to Krona text format")
    parser.add_argument("-i", "--input", required=True, help="Input Kraken2 report file (.report)")
    parser.add_argument("-o", "--output", required=True, help="Output Krona text file (.krona.txt)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"Error: input file {args.input} does not exist.")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    count = convert_report(args.input, args.output)
    print(f"Converted {args.input} -> {args.output} ({count} taxonomic lineages)")

if __name__ == "__main__":
    main()
