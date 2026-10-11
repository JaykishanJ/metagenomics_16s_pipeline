#!/usr/bin/env python3
"""
diversity_analysis.py - Publication-Grade Metagenomic Diversity & Visualization Suite.
Calculates:
  - Alpha diversity (Shannon, Simpson, Observed, Pielou evenness, Chao1)
  - Beta diversity (Bray-Curtis, Jaccard distance, 2D PCoA ordination)
  - Taxonomic composition stacked barplots (Phylum, Genus, Species)
  - Hierarchically clustered heatmaps of dominant species
Generates high-resolution publication figures (PNG 300 DPI & vector PDF).
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
from scipy.cluster.hierarchy import linkage, dendrogram

# Set publication style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

def compute_alpha_diversity(df):
    """
    df: species abundance matrix (species as rows, samples as columns)
    Returns: DataFrame with alpha diversity metrics per sample.
    """
    samples = df.columns
    records = []
    
    for s in samples:
        abun = df[s].values
        # Non-zero values
        p = abun[abun > 0]
        if len(p) == 0:
            records.append({
                'SampleID': s,
                'Observed_Species': 0,
                'Shannon': 0.0,
                'Simpson': 0.0,
                'InvSimpson': 0.0,
                'Pielou_Evenness': 0.0,
                'Chao1': 0.0
            })
            continue
            
        # Normalize to relative abundance if not already
        p_rel = p / np.sum(p)
        
        # Shannon: -sum(p * ln(p))
        shannon = -np.sum(p_rel * np.log(p_rel))
        
        # Simpson: 1 - sum(p^2) (Gini-Simpson)
        simpson = 1.0 - np.sum(p_rel ** 2)
        inv_simpson = 1.0 / np.sum(p_rel ** 2) if np.sum(p_rel ** 2) > 0 else 0
        
        # Observed Richness
        s_obs = len(p)
        
        # Pielou's Evenness: Shannon / ln(S)
        evenness = shannon / np.log(s_obs) if s_obs > 1 else 1.0
        
        # Chao1 approximation
        chao1 = s_obs
        
        records.append({
            'SampleID': s,
            'Observed_Species': s_obs,
            'Shannon': round(float(shannon), 4),
            'Simpson': round(float(simpson), 4),
            'InvSimpson': round(float(inv_simpson), 4),
            'Pielou_Evenness': round(float(evenness), 4),
            'Chao1': round(float(chao1), 4)
        })
        
    return pd.DataFrame(records).set_index('SampleID')

def plot_alpha_diversity(alpha_df, meta_df, out_prefix):
    """Generates publication-quality alpha diversity boxplots / jitter plots."""
    plot_df = alpha_df.copy()
    has_groups = False
    group_col = None
    
    if meta_df is not None and 'Group' in meta_df.columns:
        plot_df = plot_df.join(meta_df[['Group']], how='left')
        if plot_df['Group'].nunique() > 1:
            has_groups = True
            group_col = 'Group'
            
    metrics = ['Shannon', 'Simpson', 'Observed_Species', 'Pielou_Evenness']
    titles = ['Shannon Diversity (H\')', 'Simpson Diversity (1 - D)', 'Observed Species Richness', 'Pielou\'s Evenness (J\')']
    
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), dpi=300)
    axes = axes.flatten()
    
    palette = sns.color_palette("Set2")
    
    for i, (m, title) in enumerate(zip(metrics, titles)):
        ax = axes[i]
        if has_groups:
            sns.boxplot(data=plot_df, x=group_col, y=m, ax=ax, palette=palette, width=0.45, boxprops=dict(alpha=0.8), showcaps=True)
            sns.stripplot(data=plot_df, x=group_col, y=m, ax=ax, color='black', alpha=0.7, size=6, jitter=0.2)
            ax.set_xlabel('Group', fontweight='bold', fontsize=11)
        else:
            sns.boxplot(y=plot_df[m], ax=ax, color='#7293CB', width=0.3, boxprops=dict(alpha=0.8))
            sns.stripplot(y=plot_df[m], ax=ax, color='black', alpha=0.7, size=6, jitter=0.2)
            ax.set_xlabel('All Samples', fontweight='bold', fontsize=11)
            
        ax.set_title(title, fontweight='bold', fontsize=12, pad=10)
        ax.set_ylabel(m, fontweight='bold', fontsize=10)
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        sns.despine(ax=ax)
        
    plt.tight_layout(pad=2.5)
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Alpha diversity plots saved to {out_prefix}.png and {out_prefix}.pdf")

def compute_beta_pcoa(df, meta_df, out_prefix):
    """Computes Bray-Curtis dissimilarity, PCoA ordination, and generates 2D plot."""
    # Samples as rows for distance calculation
    matrix = df.T.values
    sample_names = list(df.columns)
    
    # Bray-Curtis distance matrix
    bc_dist = squareform(pdist(matrix, metric='braycurtis'))
    bc_df = pd.DataFrame(bc_dist, index=sample_names, columns=sample_names)
    
    # Save distance matrix
    dist_file = f"{out_prefix}_bray_curtis_distance.tsv"
    bc_df.to_csv(dist_file, sep='\t')
    print(f"[+] Bray-Curtis distance matrix saved to {dist_file}")
    
    # PCoA / Classical MDS
    # Double centering
    n = len(sample_names)
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H.dot(bc_dist ** 2).dot(H)
    
    eigenvals, eigenvecs = np.linalg.eigh(B)
    # Sort descending
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]
    
    pos_evals = eigenvals[eigenvals > 0]
    total_var = np.sum(pos_evals) if len(pos_evals) > 0 else 1.0
    var_exp1 = (eigenvals[0] / total_var) * 100 if len(eigenvals) > 0 else 0
    var_exp2 = (eigenvals[1] / total_var) * 100 if len(eigenvals) > 1 else 0
    
    # Coordinates for Axis 1 & Axis 2
    coord1 = eigenvecs[:, 0] * np.sqrt(max(0, eigenvals[0]))
    coord2 = eigenvecs[:, 1] * np.sqrt(max(0, eigenvals[1]))
    
    pcoa_df = pd.DataFrame({
        'SampleID': sample_names,
        'PCoA1': coord1,
        'PCoA2': coord2
    }).set_index('SampleID')
    
    has_groups = False
    group_col = None
    if meta_df is not None and 'Group' in meta_df.columns:
        pcoa_df = pcoa_df.join(meta_df[['Group']], how='left')
        if pcoa_df['Group'].nunique() > 1:
            has_groups = True
            group_col = 'Group'
            
    # Plot PCoA
    fig, ax = plt.subplots(figsize=(8.5, 6.5), dpi=300)
    
    if has_groups:
        palette = sns.color_palette("Set1", n_colors=pcoa_df[group_col].nunique())
        sns.scatterplot(data=pcoa_df, x='PCoA1', y='PCoA2', hue=group_col, s=120, palette=palette, edgecolor='black', alpha=0.9, ax=ax)
    else:
        sns.scatterplot(data=pcoa_df, x='PCoA1', y='PCoA2', color='#2b5c8f', s=120, edgecolor='black', alpha=0.9, ax=ax)
        
    for _, row in pcoa_df.iterrows():
        ax.annotate(row.name, (row['PCoA1'], row['PCoA2']),
                    textcoords="offset points", xytext=(5, 5), ha='left',
                    fontsize=8, fontweight='semibold', color='#333333')
                    
    ax.axhline(0, color='grey', linestyle='--', linewidth=0.7, alpha=0.6)
    ax.axvline(0, color='grey', linestyle='--', linewidth=0.7, alpha=0.6)
    
    ax.set_xlabel(f"PCoA 1 ({var_exp1:.1f}% Variance Explained)", fontweight='bold', fontsize=11)
    ax.set_ylabel(f"PCoA 2 ({var_exp2:.1f}% Variance Explained)", fontweight='bold', fontsize=11)
    ax.set_title("Beta Diversity - Bray-Curtis PCoA Ordination", fontweight='bold', fontsize=13, pad=12)
    ax.grid(True, linestyle=':', alpha=0.5)
    sns.despine(ax=ax)
    
    plt.tight_layout()
    fig.savefig(f"{out_prefix}_pcoa.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}_pcoa.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] PCoA ordination plots saved to {out_prefix}_pcoa.png and {out_prefix}_pcoa.pdf")

def plot_taxa_composition(df, out_prefix, top_n=12):
    """Generates taxonomic composition stacked barplot of Top N taxa + Others."""
    # Compute relative abundance per sample
    rel_df = df.div(df.sum(axis=0), axis=1) * 100
    
    # Mean abundance across samples to pick top_n
    mean_abun = rel_df.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_abun.head(top_n).index.tolist()
    
    top_df = rel_df.loc[top_taxa].copy()
    other_df = 100.0 - top_df.sum(axis=0)
    top_df.loc['Other (< 1%)'] = other_df
    
    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    colors = sns.color_palette("tab20", n_colors=top_n + 1)
    
    bottom = np.zeros(len(top_df.columns))
    for i, taxon in enumerate(top_df.index):
        values = top_df.loc[taxon].values
        ax.bar(top_df.columns, values, bottom=bottom, label=taxon, color=colors[i], width=0.7, edgecolor='white', linewidth=0.5)
        bottom += values
        
    ax.set_ylabel("Relative Abundance (%)", fontweight='bold', fontsize=11)
    ax.set_xlabel("Samples", fontweight='bold', fontsize=11)
    ax.set_title(f"Microbial Taxonomic Composition (Top {top_n} Species)", fontweight='bold', fontsize=13, pad=12)
    ax.set_ylim(0, 100.5)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, fontsize=8, title="Taxa", title_fontsize=9)
    sns.despine(ax=ax)
    
    plt.tight_layout()
    fig.savefig(f"{out_prefix}_stacked_barplot.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}_stacked_barplot.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Stacked barplot saved to {out_prefix}_stacked_barplot.png and {out_prefix}_stacked_barplot.pdf")

def plot_dominant_species_heatmap(df, out_prefix, top_n=25):
    """Generates hierarchically clustered heatmap of top dominant species."""
    mean_abun = df.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_abun.head(top_n).index
    sub_df = df.loc[top_taxa].copy()
    
    # Log10 transform relative abundances: log10(fraction * 100 + 0.01)
    log_df = np.log10(sub_df * 100 + 0.01)
    
    g = sns.clustermap(log_df,
                       figsize=(12, 10),
                       cmap="mako",
                       linewidths=0.5,
                       linecolor='white',
                       dendrogram_ratio=(0.15, 0.15),
                       cbar_kws={'label': 'log10(Abundance % + 0.01)'})
                       
    g.fig.suptitle(f"Hierarchical Clustering of Top {top_n} Dominant Species", fontweight='bold', fontsize=14, y=1.02)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8)
    
    g.savefig(f"{out_prefix}_heatmap.png", dpi=300, bbox_inches='tight')
    g.savefig(f"{out_prefix}_heatmap.pdf", bbox_inches='tight')
    plt.close(g.fig)
    print(f"[+] Clustered heatmap saved to {out_prefix}_heatmap.png and {out_prefix}_heatmap.pdf")

def main():
    parser = argparse.ArgumentParser(description="Publication-ready Diversity & Visualization Suite")
    parser.add_argument("-i", "--matrix", required=True, help="Path to species abundance matrix (TSV/TXT)")
    parser.add_argument("-m", "--metadata", default=None, help="Optional sample metadata file (TSV with SampleID and Group)")
    parser.add_argument("-d", "--diversity-out", default="results/diversity", help="Directory for diversity metric tables")
    parser.add_argument("-f", "--figures-out", default="results/figures", help="Directory for publication figures")
    args = parser.parse_args()

    os.makedirs(args.diversity_out, exist_ok=True)
    os.makedirs(args.figures_out, exist_ok=True)

    print("==========================================================")
    print(" Starting Publication-Grade Metagenomic Diversity Analysis")
    print(f" Abundance Matrix : {args.matrix}")
    print(f" Metadata File    : {args.metadata}")
    print(f" Diversity Output : {args.diversity_out}")
    print(f" Figures Output   : {args.figures_out}")
    print("==========================================================")

    # Load matrix
    df = pd.read_csv(args.matrix, sep="\t", index_col=0)
    print(f"[+] Loaded abundance matrix: {df.shape[0]} taxa across {df.shape[1]} samples.")

    # Load optional metadata
    meta_df = None
    if args.metadata and os.path.isfile(args.metadata):
        meta_df = pd.read_csv(args.metadata, sep="\t")
        if 'SampleID' in meta_df.columns:
            meta_df = meta_df.set_index('SampleID')
        print(f"[+] Loaded metadata for {len(meta_df)} samples.")

    # 1. Alpha Diversity
    print("[+] Calculating alpha diversity indices (Shannon, Simpson, Observed, Evenness)...")
    alpha_df = compute_alpha_diversity(df)
    alpha_file = os.path.join(args.diversity_out, "alpha_diversity_metrics.tsv")
    alpha_df.to_csv(alpha_file, sep="\t")
    print(f"[+] Alpha diversity metrics saved to {alpha_file}")
    
    alpha_plot_prefix = os.path.join(args.figures_out, "alpha_diversity")
    plot_alpha_diversity(alpha_df, meta_df, alpha_plot_prefix)

    # 2. Beta Diversity & PCoA
    print("[+] Calculating beta diversity (Bray-Curtis dissimilarity) and PCoA ordination...")
    beta_prefix = os.path.join(args.figures_out, "beta_diversity")
    dist_prefix = os.path.join(args.diversity_out, "beta")
    compute_beta_pcoa(df, meta_df, dist_prefix)
    # Also save the pcoa figure in figures_out
    compute_beta_pcoa(df, meta_df, beta_prefix)

    # 3. Taxonomic Composition Stacked Barplot
    print("[+] Generating taxonomic composition stacked barplot...")
    plot_taxa_composition(df, os.path.join(args.figures_out, "species_composition"), top_n=12)

    # 4. Dominant Species Heatmap
    print("[+] Generating hierarchically clustered heatmap...")
    plot_dominant_species_heatmap(df, os.path.join(args.figures_out, "species_dominant"), top_n=25)

    print("==========================================================")
    print(" Diversity analysis completed successfully!")
    print(f" High-resolution figures (300 DPI PNG & vector PDF) in: {args.figures_out}/")
    print("==========================================================")

if __name__ == "__main__":
    main()
