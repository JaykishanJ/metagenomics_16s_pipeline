#!/usr/bin/env python3
"""
plot_multilevel_taxa.py - Multi-Rank Taxonomic Composition & Group-Mean Comparisons.
Extracts:
  - Phylum-level abundance profiles
  - Genus-level abundance profiles
Generates:
  - Phylum-level stacked bar plot across samples (PNG 300 DPI & vector PDF)
  - Genus-level stacked bar plot across samples (PNG 300 DPI & vector PDF)
  - Group-mean composition comparative bar plot (Cohort A vs Cohort B)
  - Phylum and Genus abundance matrices
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

def parse_kraken_report_ranks(report_path):
    phyla = {}
    genera = {}
    with open(report_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            rank = parts[3]
            clade_count = int(parts[1])
            name = parts[5].strip()
            if rank == "P" and name.lower() not in ("unclassified", "root"):
                phyla[name] = clade_count
            elif rank == "G" and name.lower() not in ("unclassified", "root"):
                genera[name] = clade_count
    return phyla, genera

def plot_stacked_bars(df, title, ylabel, out_prefix, top_n=10, palette_name="tab20"):
    # Convert to relative abundance %
    rel_df = df.div(df.sum(axis=0), axis=1) * 100
    mean_abun = rel_df.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_abun.head(top_n).index.tolist()

    top_df = rel_df.loc[top_taxa].copy()
    other_df = 100.0 - top_df.sum(axis=0)
    top_df.loc['Other (< 1%)'] = other_df

    fig, ax = plt.subplots(figsize=(12, 6.5), dpi=300)
    colors = sns.color_palette(palette_name, n_colors=top_n + 1)

    bottom = np.zeros(len(top_df.columns))
    for i, taxon in enumerate(top_df.index):
        vals = top_df.loc[taxon].values
        ax.bar(top_df.columns, vals, bottom=bottom, label=taxon, color=colors[i], width=0.68, edgecolor='white', linewidth=0.5)
        bottom += vals

    ax.set_ylabel(ylabel, fontweight='bold', fontsize=11)
    ax.set_xlabel("Samples", fontweight='bold', fontsize=11)
    ax.set_title(title, fontweight='bold', fontsize=13, pad=12)
    ax.set_ylim(0, 100.5)
    plt.xticks(rotation=45, ha='right', fontsize=8.5)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, fontsize=8.5, title="Taxa", title_fontsize=9.5)
    sns.despine(ax=ax)

    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Stacked barplot saved to {out_prefix}.png and .pdf")

def plot_group_mean(df, meta_df, group_col, out_prefix, top_n=8):
    if meta_df is None or group_col not in meta_df.columns:
        return
    common = [s for s in df.columns if s in meta_df.index]
    sub_df = df[common].div(df[common].sum(axis=0), axis=1) * 100
    sub_meta = meta_df.loc[common]

    groups = sub_meta[group_col].unique()
    if len(groups) < 2:
        return

    # Calculate mean per group
    group_means = {}
    for g in groups:
        s_g = sub_meta[sub_meta[group_col] == g].index.tolist()
        group_means[g] = sub_df[s_g].mean(axis=1)

    mean_df = pd.DataFrame(group_means)
    top_taxa = mean_df.mean(axis=1).sort_values(ascending=False).head(top_n).index.tolist()

    plot_df = mean_df.loc[top_taxa].copy()
    plot_df.loc['Other (< 1%)'] = 100.0 - plot_df.sum(axis=0)

    fig, ax = plt.subplots(figsize=(7, 6.5), dpi=300)
    colors = sns.color_palette("tab10", n_colors=top_n + 1)

    bottom = np.zeros(len(plot_df.columns))
    for i, taxon in enumerate(plot_df.index):
        vals = plot_df.loc[taxon].values
        ax.bar(plot_df.columns, vals, bottom=bottom, label=taxon, color=colors[i], width=0.45, edgecolor='white', linewidth=0.6)
        bottom += vals

    ax.set_ylabel("Mean Relative Abundance (%)", fontweight='bold', fontsize=11)
    ax.set_xlabel("Experimental Cohort", fontweight='bold', fontsize=11)
    ax.set_title("Group-Mean Comparative Composition (Phylum)", fontweight='bold', fontsize=13, pad=12)
    ax.set_ylim(0, 100.5)
    ax.legend(bbox_to_anchor=(1.03, 1), loc='upper left', frameon=True, fontsize=8.5, title="Phylum", title_fontsize=9.5)
    sns.despine(ax=ax)

    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Group-mean composition plot saved to {out_prefix}.png and .pdf")

def main():
    parser = argparse.ArgumentParser(description="Multi-Rank Composition Analysis")
    parser.add_argument("-r", "--reports-dir", required=True, help="Directory containing .report files")
    parser.add_argument("-m", "--metadata", default=None, help="Sample metadata TSV")
    parser.add_argument("-t", "--taxa-out", default="results/taxonomy", help="Directory for taxonomy matrices")
    parser.add_argument("-f", "--figures-out", default="results/figures", help="Directory for figures")
    args = parser.parse_args()

    os.makedirs(args.taxa_out, exist_ok=True)
    os.makedirs(args.figures_out, exist_ok=True)

    report_files = sorted(glob.glob(os.path.join(args.reports_dir, "*.report")))
    clean_reports = [r for r in report_files if not r.endswith("_bracken_species.report")]

    all_phyla = set()
    all_genera = set()
    sample_phyla = {}
    sample_genera = {}

    for r in clean_reports:
        sname = os.path.basename(r).split(".report")[0]
        p_dict, g_dict = parse_kraken_report_ranks(r)
        sample_phyla[sname] = p_dict
        sample_genera[sname] = g_dict
        all_phyla.update(p_dict.keys())
        all_genera.update(g_dict.keys())

    samples = sorted(sample_phyla.keys())

    # Build Phylum DataFrame
    phyla_mat = {s: [sample_phyla[s].get(p, 0) for p in sorted(all_phyla)] for s in samples}
    df_phyla = pd.DataFrame(phyla_mat, index=sorted(all_phyla))
    phyla_file = os.path.join(args.taxa_out, "phylum_abundance_matrix.tsv")
    df_phyla.to_csv(phyla_file, sep="\t")

    # Build Genus DataFrame
    genus_mat = {s: [sample_genera[s].get(g, 0) for g in sorted(all_genera)] for s in samples}
    df_genus = pd.DataFrame(genus_mat, index=sorted(all_genera))
    genus_file = os.path.join(args.taxa_out, "genus_abundance_matrix.tsv")
    df_genus.to_csv(genus_file, sep="\t")

    print(f"[+] Phylum abundance matrix saved to {phyla_file} ({df_phyla.shape[0]} phyla)")
    print(f"[+] Genus abundance matrix saved to {genus_file} ({df_genus.shape[0]} genera)")

    # Plot Phylum Composition
    plot_stacked_bars(df_phyla, "Microbial Taxonomic Composition (Phylum Level)",
                      "Relative Abundance (%)", os.path.join(args.figures_out, "taxa_composition_phylum"),
                      top_n=8, palette_name="Set3")

    # Plot Genus Composition
    plot_stacked_bars(df_genus, "Microbial Taxonomic Composition (Top 15 Genera)",
                      "Relative Abundance (%)", os.path.join(args.figures_out, "taxa_composition_genus"),
                      top_n=15, palette_name="tab20")

    # Plot Group-Mean Comparison if metadata available
    if args.metadata and os.path.isfile(args.metadata):
        meta_df = pd.read_csv(args.metadata, sep="\t")
        if 'SampleID' in meta_df.columns:
            meta_df = meta_df.set_index('SampleID')
        plot_group_mean(df_phyla, meta_df, 'Group', os.path.join(args.figures_out, "group_mean_composition_phylum"))

if __name__ == "__main__":
    main()
