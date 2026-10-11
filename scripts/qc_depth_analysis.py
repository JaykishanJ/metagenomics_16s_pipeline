#!/usr/bin/env python3
"""
qc_depth_analysis.py - Sequencing Depth, QC Metrics & Rarefaction Curves.
Parses Kraken reports to extract:
  - Total read counts per sample
  - Classified vs unclassified reads and percentages
Generates:
  - Two-panel sequencing depth and classification efficiency figure (PNG 300 DPI & vector PDF)
  - Rarefaction curves (Observed Richness & Shannon index vs Sequencing Depth)
  - Summary table: sequencing_depth_qc_summary.tsv
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

def parse_kraken_report_qc(report_path):
    total_unclassified = 0
    total_classified = 0
    with open(report_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            rank = parts[3]
            clade_count = int(parts[1])
            name = parts[5].strip()
            if rank == "U" or name == "unclassified":
                total_unclassified = clade_count
            elif rank == "R" or name == "root":
                total_classified = clade_count
    total_reads = total_unclassified + total_classified
    pct_classified = (total_classified / total_reads * 100) if total_reads > 0 else 0
    return {
        'Total_Reads': total_reads,
        'Classified_Reads': total_classified,
        'Unclassified_Reads': total_unclassified,
        'Pct_Classified': round(pct_classified, 2)
    }

def compute_rarefaction(matrix_df, steps=10, reps=10):
    """
    Subsamples read abundance profiles at multiple depth fractions to generate
    rarefaction curves (Observed Richness and Shannon Diversity).
    """
    depth_fractions = np.linspace(0.1, 1.0, steps)
    samples = matrix_df.columns
    rarefaction_records = []

    for s in samples:
        # Approximate integer counts (assuming 50,000 depth scale for proportions)
        profile = matrix_df[s].values
        # non-zero taxa
        profile = profile[profile > 0]
        if len(profile) == 0:
            continue
        p_rel = profile / np.sum(profile)
        # Scaled pseudo-counts
        total_counts = 50000
        taxa_counts = np.round(p_rel * total_counts).astype(int)
        taxa_counts = taxa_counts[taxa_counts > 0]

        for frac in depth_fractions:
            subsample_size = int(total_counts * frac)
            obs_list = []
            shannon_list = []
            for _ in range(reps):
                # Multinomial sampling
                sampled = np.random.multinomial(subsample_size, p_rel)
                sampled_nz = sampled[sampled > 0]
                obs_list.append(len(sampled_nz))
                p_sub = sampled_nz / subsample_size
                shannon_list.append(-np.sum(p_sub * np.log(p_sub)))
                
            rarefaction_records.append({
                'SampleID': s,
                'Sequencing_Depth': subsample_size,
                'Fraction': frac,
                'Mean_Observed': np.mean(obs_list),
                'Mean_Shannon': np.mean(shannon_list)
            })

    return pd.DataFrame(rarefaction_records)

def plot_qc_sequencing_depth(qc_df, out_prefix):
    """Generates a two-panel publication plot: reads per sample and % classified."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # Panel A: Total Reads per Sample
    palette = sns.color_palette("mako", len(qc_df))
    bars = ax1.bar(qc_df['SampleID'], qc_df['Total_Reads'], color=palette, edgecolor='#333333', linewidth=0.6)
    median_depth = qc_df['Total_Reads'].median()
    ax1.axhline(median_depth, color='red', linestyle='--', linewidth=1.2, label=f"Median Depth ({median_depth:,.0f})")
    ax1.set_ylabel("Total Paired Reads", fontweight='bold', fontsize=11)
    ax1.set_title("A. Sequencing Depth per Sample", fontweight='bold', fontsize=13, pad=12)
    ax1.set_xticklabels(qc_df['SampleID'], rotation=45, ha='right', fontsize=8.5)
    ax1.legend(frameon=True, fontsize=9)
    ax1.grid(axis='y', linestyle=':', alpha=0.5)
    sns.despine(ax=ax1)

    # Panel B: Classified vs Unclassified Reads Proportion
    x = np.arange(len(qc_df))
    width = 0.65
    ax2.bar(x, qc_df['Pct_Classified'], width, label='Classified Reads', color='#2ca02c', edgecolor='white', linewidth=0.5)
    ax2.bar(x, 100.0 - qc_df['Pct_Classified'], width, bottom=qc_df['Pct_Classified'], label='Unclassified Reads', color='#d62728', edgecolor='white', linewidth=0.5)
    ax2.set_ylabel("Percentage of Total Reads (%)", fontweight='bold', fontsize=11)
    ax2.set_title("B. Taxonomic Classification Efficiency", fontweight='bold', fontsize=13, pad=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(qc_df['SampleID'], rotation=45, ha='right', fontsize=8.5)
    ax2.set_ylim(0, 100)
    ax2.legend(loc='lower left', frameon=True, fontsize=9)
    ax2.grid(axis='y', linestyle=':', alpha=0.5)
    sns.despine(ax=ax2)

    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] QC Sequencing depth plot saved to {out_prefix}.png and .pdf")

def plot_rarefaction_curves(rare_df, out_prefix):
    """Generates two-panel rarefaction curves for richness and Shannon diversity."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)

    # Panel A: Observed Richness vs Depth
    sns.lineplot(data=rare_df, x='Sequencing_Depth', y='Mean_Observed', hue='SampleID',
                 legend=False, ax=ax1, alpha=0.8, linewidth=1.5, palette='tab20')
    ax1.set_xlabel("Sequencing Depth (Subsampled Reads)", fontweight='bold', fontsize=10)
    ax1.set_ylabel("Observed Taxa Richness (S)", fontweight='bold', fontsize=10)
    ax1.set_title("A. Rarefaction Curve - Observed Richness", fontweight='bold', fontsize=12, pad=10)
    ax1.grid(True, linestyle=':', alpha=0.5)
    sns.despine(ax=ax1)

    # Panel B: Shannon Index vs Depth
    sns.lineplot(data=rare_df, x='Sequencing_Depth', y='Mean_Shannon', hue='SampleID',
                 legend=False, ax=ax2, alpha=0.8, linewidth=1.5, palette='tab20')
    ax2.set_xlabel("Sequencing Depth (Subsampled Reads)", fontweight='bold', fontsize=10)
    ax2.set_ylabel("Shannon Diversity Index (H')", fontweight='bold', fontsize=10)
    ax2.set_title("B. Rarefaction Curve - Shannon Diversity", fontweight='bold', fontsize=12, pad=10)
    ax2.grid(True, linestyle=':', alpha=0.5)
    sns.despine(ax=ax2)

    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Rarefaction curves saved to {out_prefix}.png and .pdf")

def main():
    parser = argparse.ArgumentParser(description="QC Sequencing Depth & Rarefaction Curves")
    parser.add_argument("-r", "--reports-dir", required=True, help="Directory containing .report files")
    parser.add_argument("-m", "--matrix", default=None, help="Species abundance matrix for rarefaction")
    parser.add_argument("-o", "--qc-out", default="results/qc", help="Directory for QC tables")
    parser.add_argument("-f", "--figures-out", default="results/figures", help="Directory for figures")
    args = parser.parse_args()

    os.makedirs(args.qc_out, exist_ok=True)
    os.makedirs(args.figures_out, exist_ok=True)

    # Parse reports
    report_files = sorted(glob.glob(os.path.join(args.reports_dir, "*.report")))
    clean_reports = [r for r in report_files if not r.endswith("_bracken_species.report")]

    if not clean_reports:
        sys.exit(f"[-] Error: No valid .report files found in {args.reports_dir}")

    records = []
    for r in clean_reports:
        sname = os.path.basename(r).split(".report")[0]
        data = parse_kraken_report_qc(r)
        data['SampleID'] = sname
        records.append(data)

    qc_df = pd.DataFrame(records).sort_values('SampleID')
    summary_file = os.path.join(args.qc_out, "sequencing_depth_qc_summary.tsv")
    qc_df.to_csv(summary_file, sep="\t", index=False)
    print(f"[+] QC summary table saved to {summary_file}")

    # Plot 1: Sequencing depth and classification
    depth_plot_prefix = os.path.join(args.figures_out, "qc_sequencing_depth_and_classification")
    plot_qc_sequencing_depth(qc_df, depth_plot_prefix)

    # Plot 2: Rarefaction curves
    if args.matrix and os.path.isfile(args.matrix):
        print("[+] Calculating in-silico rarefaction profiles...")
        matrix_df = pd.read_csv(args.matrix, sep="\t", index_col=0)
        rare_df = compute_rarefaction(matrix_df, steps=10, reps=5)
        rare_file = os.path.join(args.qc_out, "rarefaction_metrics.tsv")
        rare_df.to_csv(rare_file, sep="\t", index=False)
        plot_rarefaction_curves(rare_df, os.path.join(args.figures_out, "qc_rarefaction_curves"))

if __name__ == "__main__":
    main()
