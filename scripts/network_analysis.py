#!/usr/bin/env python3
"""
network_analysis.py - Microbial Co-occurrence Network & Correlation Suite.
Calculates:
  1. Pairwise Spearman rank correlations between dominant taxa
  2. Co-occurrence correlation heatmap (hierarchically clustered)
  3. Association network plot (nodes = taxa, edges = significant positive/negative correlations)
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
from scipy import stats

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

def compute_spearman_network(df, top_n=20, r_thresh=0.45, p_thresh=0.05):
    # Select top N taxa by mean abundance
    mean_abun = df.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_abun.head(top_n).index.tolist()
    sub_df = df.loc[top_taxa]
    
    n = len(top_taxa)
    corr_matrix = np.zeros((n, n))
    pval_matrix = np.zeros((n, n))
    
    edge_list = []
    
    for i in range(n):
        for j in range(n):
            if i == j:
                corr_matrix[i, j] = 1.0
                pval_matrix[i, j] = 0.0
            else:
                rho, p = stats.spearmanr(sub_df.iloc[i].values, sub_df.iloc[j].values)
                corr_matrix[i, j] = rho
                pval_matrix[i, j] = p
                
                if i < j and abs(rho) >= r_thresh and p < p_thresh:
                    edge_list.append({
                        'Taxon1': top_taxa[i],
                        'Taxon2': top_taxa[j],
                        'Spearman_Rho': round(float(rho), 4),
                        'P_Value': float(p),
                        'Type': 'Positive' if rho > 0 else 'Negative'
                    })
                    
    corr_df = pd.DataFrame(corr_matrix, index=top_taxa, columns=top_taxa)
    pval_df = pd.DataFrame(pval_matrix, index=top_taxa, columns=top_taxa)
    edges_df = pd.DataFrame(edge_list)
    
    return corr_df, pval_df, edges_df, top_taxa, mean_abun.loc[top_taxa]

def plot_correlation_heatmap(corr_df, out_prefix):
    g = sns.clustermap(corr_df,
                       figsize=(11, 9.5),
                       cmap="vlag",
                       vmin=-1, vmax=1,
                       linewidths=0.5,
                       linecolor='white',
                       dendrogram_ratio=(0.15, 0.15),
                       cbar_kws={'label': "Spearman's Rank Correlation (rho)"})
                       
    g.fig.suptitle("Microbial Co-occurrence Correlation Heatmap", fontweight='bold', fontsize=13, y=1.02)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8.5)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8)
    
    g.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    g.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(g.fig)
    print(f"[+] Correlation heatmap saved to {out_prefix}.png and .pdf")

def plot_network(edges_df, taxa, abundances, out_prefix):
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    
    n = len(taxa)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    radius = 1.0
    pos = {taxa[i]: (radius * np.cos(angles[i]), radius * np.sin(angles[i])) for i in range(n)}
    
    # Draw edges
    for _, edge in edges_df.iterrows():
        p1 = pos[edge['Taxon1']]
        p2 = pos[edge['Taxon2']]
        color = '#2ca02c' if edge['Type'] == 'Positive' else '#d62728'
        linestyle = '-' if edge['Type'] == 'Positive' else '--'
        lw = 1.0 + abs(edge['Spearman_Rho']) * 2.5
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linestyle=linestyle, linewidth=lw, alpha=0.6, zorder=1)
        
    # Draw nodes
    abun_norm = (abundances / abundances.max()) * 400 + 150
    for taxon in taxa:
        p = pos[taxon]
        s = abun_norm.get(taxon, 200)
        ax.scatter(p[0], p[1], s=s, color='#1f77b4', edgecolors='#333333', linewidth=1.2, zorder=2)
        
        # Label with offset
        ha = 'left' if p[0] >= 0 else 'right'
        offset = 0.12 * np.sign(p[0]) if p[0] != 0 else 0.05
        ax.text(p[0] * 1.15, p[1] * 1.15, taxon, fontsize=8, fontweight='semibold',
                ha=ha, va='center', color='#222222', zorder=3)
                
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.axis('off')
    
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='#2ca02c', lw=2, label='Positive Co-occurrence (rho >= 0.45)'),
        Line2D([0], [0], color='#d62728', lw=2, linestyle='--', label='Co-exclusion / Competition (rho <= -0.45)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markeredgecolor='#333333', markersize=10, label='Taxon (Node size ~ Abundance)')
    ]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=True, fontsize=8.5)
    
    ax.set_title("Microbial Co-occurrence & Association Network", fontweight='bold', fontsize=14, pad=15)
    
    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Co-occurrence network saved to {out_prefix}.png and .pdf")

def main():
    parser = argparse.ArgumentParser(description="Microbial Co-occurrence Network Analysis")
    parser.add_argument("-i", "--matrix", required=True, help="Species abundance matrix")
    parser.add_argument("-o", "--net-out", default="results/diversity", help="Directory for tables")
    parser.add_argument("-f", "--figures-out", default="results/figures", help="Directory for figures")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top taxa to include")
    args = parser.parse_args()

    os.makedirs(args.net_out, exist_ok=True)
    os.makedirs(args.figures_out, exist_ok=True)

    df = pd.read_csv(args.matrix, sep="\t", index_col=0)

    print("==========================================================")
    print(" Microbial Co-occurrence Network & Association Suite")
    print(f" Abundance Matrix : {args.matrix}")
    print(f" Top Taxa Evaluated: {args.top_n}")
    print("==========================================================")

    corr_df, pval_df, edges_df, taxa, abundances = compute_spearman_network(df, top_n=args.top_n)

    corr_file = os.path.join(args.net_out, "taxa_spearman_correlations.tsv")
    edges_file = os.path.join(args.net_out, "network_edges.tsv")
    corr_df.to_csv(corr_file, sep="\t")
    edges_df.to_csv(edges_file, sep="\t", index=False)
    print(f"[+] Correlation matrix saved to {corr_file}")
    print(f"[+] Significant network edges saved to {edges_file} ({len(edges_df)} edges)")

    plot_correlation_heatmap(corr_df, os.path.join(args.figures_out, "taxa_correlation_heatmap"))
    plot_network(edges_df, taxa, abundances, os.path.join(args.figures_out, "taxa_cooccurrence_network"))

if __name__ == "__main__":
    main()
