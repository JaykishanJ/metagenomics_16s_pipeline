#!/usr/bin/env python3
"""
plot_functional_pathways.py - Plot predicted MetaCyc functional pathways from PICRUSt2.
Extracts top 20 dominant pathways across samples, creates publication-ready
clustered heatmap and relative abundance stacked barplot (PNG 300 DPI & vector PDF).
"""

import os
import sys
import gzip
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

def load_pathway_table(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            df = pd.read_csv(f, sep="\t", index_col=0)
    else:
        df = pd.read_csv(path, sep="\t", index_col=0)
    
    # If description column exists
    desc_col = None
    for c in ['description', 'Description']:
        if c in df.columns:
            desc_col = df[c]
            df = df.drop(columns=[c])
            break
            
    return df, desc_col

def main():
    parser = argparse.ArgumentParser(description="Plot PICRUSt2 functional pathway predictions")
    parser.add_argument("-i", "--pathways", required=True, help="Path to path_abun_unstrat.tsv or .tsv.gz")
    parser.add_argument("-f", "--figures-out", default="results/figures", help="Directory for figures")
    parser.add_argument("-o", "--tables-out", default="results/functional", help="Directory for tables")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top pathways to display")
    args = parser.parse_args()

    os.makedirs(args.figures_out, exist_ok=True)
    os.makedirs(args.tables_out, exist_ok=True)

    if not os.path.isfile(args.pathways):
        sys.exit(f"[-] Error: Pathway file {args.pathways} not found.")

    df, desc = load_pathway_table(args.pathways)
    print(f"[+] Loaded functional pathway table: {df.shape[0]} pathways across {df.shape[1]} samples.")

    # Calculate mean abundance to find top pathways
    mean_abun = df.mean(axis=1).sort_values(ascending=False)
    top_pws = mean_abun.head(args.top_n).index

    sub_df = df.loc[top_pws].copy()

    # Rename indices with description if available
    if desc is not None:
        new_names = []
        for pw in sub_df.index:
            d_str = str(desc.get(pw, ""))
            if d_str and d_str != "nan":
                # Truncate if too long
                d_clean = (d_str[:45] + '...') if len(d_str) > 45 else d_str
                new_names.append(f"{pw}: {d_clean}")
            else:
                new_names.append(pw)
        sub_df.index = new_names

    # Save summary table
    summary_file = os.path.join(args.tables_out, "top_functional_pathways.tsv")
    sub_df.to_csv(summary_file, sep="\t")
    print(f"[+] Top pathways summary table saved to {summary_file}")

    # 1. Clustered Heatmap
    log_df = np.log10(sub_df + 1.0)
    g = sns.clustermap(log_df,
                       figsize=(13, 9.5),
                       cmap="YlGnBu",
                       linewidths=0.5,
                       linecolor='white',
                       dendrogram_ratio=(0.15, 0.15),
                       cbar_kws={'label': 'log10(Pathway Abundance + 1)'})

    g.fig.suptitle(f"Top {args.top_n} Predicted MetaCyc Metabolic Pathways (PICRUSt2)",
                   fontweight='bold', fontsize=13, y=1.02)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8)

    heatmap_prefix = os.path.join(args.figures_out, "functional_pathways_heatmap")
    g.savefig(f"{heatmap_prefix}.png", dpi=300, bbox_inches='tight')
    g.savefig(f"{heatmap_prefix}.pdf", bbox_inches='tight')
    plt.close(g.fig)
    print(f"[+] Functional heatmap saved to {heatmap_prefix}.png and .pdf")

    # 2. Horizontal Barplot of Average Abundance
    fig, ax = plt.subplots(figsize=(10, 7.5), dpi=300)
    mean_top = sub_df.mean(axis=1).sort_values()
    colors = sns.color_palette("viridis", len(mean_top))

    y_pos = np.arange(len(mean_top))
    ax.barh(y_pos, mean_top.values, color=colors, edgecolor='#333333', linewidth=0.5, height=0.65)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(mean_top.index, fontsize=8.5)
    ax.set_xlabel("Mean Predicted Abundance Across Samples", fontweight='bold', fontsize=10)
    ax.set_title(f"Top {args.top_n} Abundant Microbial Metabolic Pathways", fontweight='bold', fontsize=12, pad=12)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    sns.despine(ax=ax)

    plt.tight_layout()
    bar_prefix = os.path.join(args.figures_out, "functional_pathways_barplot")
    fig.savefig(f"{bar_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{bar_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Functional barplot saved to {bar_prefix}.png and .pdf")

if __name__ == "__main__":
    main()
