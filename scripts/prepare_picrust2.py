#!/usr/bin/env python3
"""
prepare_picrust2.py - Prepare representative sequences and abundance table for PICRUSt2.
Uses vsearch to:
  1. Merge paired-end amplicon reads per sample.
  2. Pool and dereplicate reads across all samples.
  3. Cluster into 97% representative OTUs / centroids (rep_seqs.fasta).
  4. Quantify OTU counts across each sample to produce otu_table.tsv (#OTU ID header).
"""

import os
import sys
import glob
import subprocess
import argparse
import pandas as pd
from collections import defaultdict

def run_cmd(cmd, desc=""):
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"[-] Error during {desc}:\n{res.stderr}", file=sys.stderr)
        sys.exit(res.returncode)
    return res.stdout

def find_vsearch():
    paths = [
        "vsearch",
        "/home/genome/miniconda3/envs/qiime2-metagenome-2024.10/bin/vsearch",
        "/home/genome/miniconda3/envs/metagenomics_16s/bin/vsearch"
    ]
    for p in paths:
        res = subprocess.run(f"which {p} || ls {p} 2>/dev/null", shell=True, stdout=subprocess.PIPE, text=True)
        out = res.stdout.strip()
        if out and os.path.exists(out):
            return out
    return "vsearch"

def main():
    parser = argparse.ArgumentParser(description="Generate PICRUSt2 inputs from FastQ files")
    parser.add_argument("-i", "--input-dir", required=True, help="Directory containing FastQ files")
    parser.add_argument("-o", "--output-dir", default="results/functional", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()

    vsearch = find_vsearch()
    os.makedirs(args.output_dir, exist_ok=True)
    scratch_dir = os.path.join(args.output_dir, "scratch")
    os.makedirs(scratch_dir, exist_ok=True)

    print("==========================================================")
    print(" Preparing PICRUSt2 Amplicon Sequences & Abundance Table")
    print(f" FastQ Directory : {args.input_dir}")
    print(f" Output Directory: {args.output_dir}")
    print(f" VSEARCH Binary  : {vsearch}")
    print("==========================================================")

    # Detect pairs
    r1_files = sorted(glob.glob(os.path.join(args.input_dir, "*_1.fastq*")) + glob.glob(os.path.join(args.input_dir, "*_R1.fastq*")))
    if not r1_files:
        sys.exit(f"[-] Error: No R1 FastQ files found in {args.input_dir}")

    samples = []
    merged_files = []

    for r1 in r1_files:
        if "_1.fastq" in r1:
            r2 = r1.replace("_1.fastq", "_2.fastq")
            sname = os.path.basename(r1).split("_1.fastq")[0]
        else:
            r2 = r1.replace("_R1.fastq", "_R2.fastq")
            sname = os.path.basename(r1).split("_R1.fastq")[0]

        if not os.path.isfile(r2):
            print(f"[!] Warning: Mate {r2} not found for {r1}, skipping...")
            continue

        samples.append(sname)
        merged_fa = os.path.join(scratch_dir, f"{sname}_merged.fa")
        merged_files.append((sname, merged_fa))

        if not os.path.isfile(merged_fa) or os.path.getsize(merged_fa) == 0:
            print(f"[+] Merging paired reads for sample {sname}...")
            cmd = f"{vsearch} --fastq_mergepairs {r1} --reverse {r2} --fastaout {merged_fa} --threads {args.threads} --quiet"
            run_cmd(cmd, f"merging {sname}")

    # Pool merged reads
    pooled_fa = os.path.join(scratch_dir, "all_pooled.fa")
    print(f"[+] Pooling merged reads across {len(samples)} samples...")
    with open(pooled_fa, "w") as out_f:
        for sname, mfa in merged_files:
            if os.path.isfile(mfa):
                with open(mfa, "r") as in_f:
                    for line in in_f:
                        out_f.write(line)

    # Dereplicate
    derep_fa = os.path.join(scratch_dir, "derep.fa")
    print("[+] Dereplicating pooled sequences...")
    cmd_derep = f"{vsearch} --derep_fulllength {pooled_fa} --output {derep_fa} --sizeout --minuniquesize 2 --threads {args.threads} --quiet"
    run_cmd(cmd_derep, "dereplicating")

    # Cluster at 97% to generate representative sequences
    raw_rep_fa = os.path.join(scratch_dir, "raw_centroids.fa")
    print("[+] Clustering into 97% representative OTUs...")
    cmd_cluster = f"{vsearch} --cluster_size {derep_fa} --id 0.97 --centroids {raw_rep_fa} --threads {args.threads} --quiet"
    run_cmd(cmd_cluster, "clustering")

    # Clean FASTA headers to OTU_1, OTU_2...
    rep_seqs_file = os.path.join(args.output_dir, "rep_seqs.fasta")
    otu_id_map = {}
    otu_counter = 1
    with open(raw_rep_fa, "r") as in_f, open(rep_seqs_file, "w") as out_f:
        for line in in_f:
            if line.startswith(">"):
                old_id = line.strip().lstrip(">").split(";")[0]
                new_id = f"OTU_{otu_counter}"
                otu_id_map[old_id] = new_id
                otu_counter += 1
                out_f.write(f">{new_id}\n")
            else:
                out_f.write(line)

    print(f"[+] Generated {otu_counter - 1} representative sequences in {rep_seqs_file}")

    # Map each sample back to rep_seqs.fasta to create abundance table
    print("[+] Quantifying OTU abundance across samples...")
    otu_counts = defaultdict(lambda: defaultdict(int))
    all_otus = [f"OTU_{i}" for i in range(1, otu_counter)]

    for sname, mfa in merged_files:
        tab_out = os.path.join(scratch_dir, f"{sname}_map.txt")
        cmd_map = f"{vsearch} --usearch_global {mfa} --db {rep_seqs_file} --id 0.97 --userout {tab_out} --userfields query+target --threads {args.threads} --quiet"
        run_cmd(cmd_map, f"mapping {sname}")

        if os.path.isfile(tab_out):
            with open(tab_out, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        target_otu = parts[1]
                        otu_counts[target_otu][sname] += 1

    # Build DataFrame
    matrix_data = []
    for otu in all_otus:
        row = {'#OTU ID': otu}
        for s in samples:
            row[s] = otu_counts[otu][s]
        matrix_data.append(row)

    df_otu = pd.DataFrame(matrix_data)
    otu_table_file = os.path.join(args.output_dir, "otu_table.tsv")
    df_otu.to_csv(otu_table_file, sep="\t", index=False)
    print(f"[+] PICRUSt2 abundance table generated: {otu_table_file}")
    print(f"[+] Matrix dimensions: {len(df_otu)} OTUs across {len(samples)} samples.")

if __name__ == "__main__":
    main()
