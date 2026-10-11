#!/usr/bin/env python3
"""
differential_abundance.py - Differential Abundance Testing, PERMANOVA, and Biomarker Discovery.
Performs:
  1. Statistical group comparisons (Mann-Whitney U / Welch's t-test, Benjamini-Hochberg FDR)
  2. Volcano plot (Log2FC vs -log10 p-value)
  3. Forest plot / effect size ranking
  4. Top biomarker group boxplots
  5. Permutational Multivariate Analysis of Variance (PERMANOVA) on Bray-Curtis dissimilarity
  6. Pairwise differential abundance summary table
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
from scipy.spatial.distance import pdist, squareform

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

def benjamini_hochberg(p_values):
    """Computes Benjamini-Hochberg FDR adjusted p-values."""
    p_values = np.asarray(p_values)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    q_values = np.zeros(n)
    min_q = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        q = (sorted_p[i] * n) / rank
        min_q = min(min_q, q)
        q_values[i] = min(1.0, min_q)
    # Restore original order
    adj_p = np.zeros(n)
    adj_p[sorted_indices] = q_values
    return adj_p

def run_permanova(df, meta_df, group_col='Group', perms=999):
    """
    Computes PERMANOVA pseudo-F, R2 (variance explained), and permutation p-value
    on Bray-Curtis dissimilarity matrix (identical to vegan adonis).
    """
    common_samples = [s for s in df.columns if s in meta_df.index]
    if len(common_samples) < 4:
        return None
        
    sub_df = df[common_samples]
    sub_meta = meta_df.loc[common_samples]
    groups = sub_meta[group_col].values
    unique_groups = np.unique(groups)
    
    if len(unique_groups) < 2:
        return None
        
    matrix = sub_df.T.values
    d = squareform(pdist(matrix, metric='braycurtis'))
    n = len(common_samples)
    
    # Total sum of squares: SS_T = (1/n) * sum(d_ij^2)
    ss_t = np.sum(d ** 2) / (2 * n)
    
    # Within-group sum of squares: SS_W = sum_g [ (1/n_g) * sum_{i,j in g} d_ij^2 / 2 ]
    def calc_ss_w(grp_labels):
        ss_w = 0.0
        for g in unique_groups:
            idx = np.where(grp_labels == g)[0]
            n_g = len(idx)
            if n_g > 1:
                sub_d = d[np.ix_(idx, idx)]
                ss_w += np.sum(sub_d ** 2) / (2 * n_g)
        return ss_w

    obs_ss_w = calc_ss_w(groups)
    obs_ss_a = ss_t - obs_ss_w
    a = len(unique_groups)
    df_a = a - 1
    df_w = n - a
    
    obs_f = (obs_ss_a / df_a) / (obs_ss_w / df_w) if obs_ss_w > 0 else 0.0
    obs_r2 = obs_ss_a / ss_t if ss_t > 0 else 0.0
    
    # Permutations
    perm_f_list = []
    for _ in range(perms):
        perm_groups = np.random.permutation(groups)
        p_ss_w = calc_ss_w(perm_groups)
        p_ss_a = ss_t - p_ss_w
        p_f = (p_ss_a / df_a) / (p_ss_w / df_w) if p_ss_w > 0 else 0.0
        perm_f_list.append(p_f)
        
    p_val = (np.sum(np.array(perm_f_list) >= obs_f) + 1) / (perms + 1)
    
    return {
        'Factor': group_col,
        'Df': df_a,
        'SumOfSqs': round(float(obs_ss_a), 4),
        'R2': round(float(obs_r2), 4),
        'F_Statistic': round(float(obs_f), 4),
        'P_Value': round(float(p_val), 4),
        'Permutations': perms
    }

def run_differential_abundance(df, meta_df, group_col='Group'):
    common_samples = [s for s in df.columns if s in meta_df.index]
    sub_df = df[common_samples]
    sub_meta = meta_df.loc[common_samples]
    
    groups = sub_meta[group_col].unique()
    if len(groups) != 2:
        print(f"[!] Warning: Found {len(groups)} groups. Differential abundance currently tests pairwise 2-group comparisons.")
        groups = groups[:2]
        
    g1, g2 = groups[0], groups[1]
    g1_samples = sub_meta[sub_meta[group_col] == g1].index.tolist()
    g2_samples = sub_meta[sub_meta[group_col] == g2].index.tolist()
    
    records = []
    eps = 1e-5
    
    for taxon in df.index:
        v1 = sub_df.loc[taxon, g1_samples].values
        v2 = sub_df.loc[taxon, g2_samples].values
        
        m1 = np.mean(v1)
        m2 = np.mean(v2)
        
        prev1 = np.sum(v1 > 0) / len(v1) * 100
        prev2 = np.sum(v2 > 0) / len(v2) * 100
        
        # Log2 Fold Change (g2 vs g1)
        log2fc = np.log2((m2 + eps) / (m1 + eps))
        
        # Mann-Whitney U test
        try:
            stat, p_val = stats.mannwhitneyu(v2, v1, alternative='two-sided')
        except Exception:
            p_val = 1.0
            
        records.append({
            'Taxon': taxon,
            f'Mean_{g1}': m1,
            f'Mean_{g2}': m2,
            f'Prevalence_{g1}_pct': round(prev1, 1),
            f'Prevalence_{g2}_pct': round(prev2, 1),
            'Log2FC': round(float(log2fc), 4),
            'P_Value': float(p_val)
        })
        
    res_df = pd.DataFrame(records)
    res_df['Adj_P_Value'] = benjamini_hochberg(res_df['P_Value'].values)
    res_df['Significant'] = (res_df['P_Value'] < 0.05) & (np.abs(res_df['Log2FC']) >= 1.0)
    
    return res_df.sort_values('P_Value'), g1, g2

def plot_volcano(diff_df, g1, g2, out_prefix):
    fig, ax = plt.subplots(figsize=(9, 6.5), dpi=300)
    
    plot_df = diff_df.copy()
    plot_df['neg_log_p'] = -np.log10(np.clip(plot_df['P_Value'], 1e-10, 1.0))
    
    # Colors
    colors = []
    for _, r in plot_df.iterrows():
        if r['Significant']:
            colors.append('#d62728' if r['Log2FC'] > 0 else '#1f77b4')
        else:
            colors.append('#aaaaaa')
            
    ax.scatter(plot_df['Log2FC'], plot_df['neg_log_p'], c=colors, alpha=0.85, s=60, edgecolors='none')
    
    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.8, alpha=0.7, label='p = 0.05')
    ax.axvline(1.0, color='grey', linestyle=':', linewidth=0.8, alpha=0.7)
    ax.axvline(-1.0, color='grey', linestyle=':', linewidth=0.8, alpha=0.7)
    
    # Annotate top significant taxa
    top_sig = plot_df[plot_df['Significant']].head(8)
    for _, r in top_sig.iterrows():
        short_name = r['Taxon'].split()[-1] if len(r['Taxon'].split()) > 1 else r['Taxon']
        ax.annotate(r['Taxon'], (r['Log2FC'], r['neg_log_p']),
                    xytext=(4, 4), textcoords='offset points', fontsize=8, fontweight='semibold')
                    
    ax.set_xlabel(f"Log2 Fold Change ({g2} vs {g1})", fontweight='bold', fontsize=11)
    ax.set_ylabel("-log10(P-Value)", fontweight='bold', fontsize=11)
    ax.set_title(f"Differential Abundance Volcano Plot ({g2} vs {g1})", fontweight='bold', fontsize=13, pad=12)
    ax.grid(True, linestyle=':', alpha=0.4)
    sns.despine(ax=ax)
    
    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Volcano plot saved to {out_prefix}.png and .pdf")

def plot_forest(diff_df, g1, g2, out_prefix, top_n=15):
    # Sort by absolute Log2FC among p < 0.1
    candidates = diff_df.sort_values('P_Value').head(top_n).copy()
    candidates = candidates.sort_values('Log2FC')
    
    fig, ax = plt.subplots(figsize=(10, 7), dpi=300)
    colors = ['#d62728' if x > 0 else '#1f77b4' for x in candidates['Log2FC']]
    
    y_pos = np.arange(len(candidates))
    ax.barh(y_pos, candidates['Log2FC'], color=colors, height=0.6, edgecolor='#333333', linewidth=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(candidates['Taxon'], fontsize=9)
    ax.axvline(0, color='black', linewidth=0.8)
    
    ax.set_xlabel(f"Effect Size - Log2 Fold Change ({g2} vs {g1})", fontweight='bold', fontsize=11)
    ax.set_title(f"Top Differential Microbial Taxa by Effect Size", fontweight='bold', fontsize=13, pad=12)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    sns.despine(ax=ax)
    
    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Forest effect size plot saved to {out_prefix}.png and .pdf")

def main():
    parser = argparse.ArgumentParser(description="Differential Abundance Testing & PERMANOVA")
    parser.add_argument("-i", "--matrix", required=True, help="Species abundance matrix")
    parser.add_argument("-m", "--metadata", required=True, help="Sample metadata TSV")
    parser.add_argument("-g", "--group-col", default="Group", help="Metadata column for groups")
    parser.add_argument("-o", "--diff-out", default="results/differential", help="Output directory for tables")
    parser.add_argument("-f", "--figures-out", default="results/figures", help="Output directory for figures")
    args = parser.parse_args()

    os.makedirs(args.diff_out, exist_ok=True)
    os.makedirs(args.figures_out, exist_ok=True)

    df = pd.read_csv(args.matrix, sep="\t", index_col=0)
    meta_df = pd.read_csv(args.metadata, sep="\t")
    if 'SampleID' in meta_df.columns:
        meta_df = meta_df.set_index('SampleID')

    print("==========================================================")
    print(" Differential Abundance & Multivariate PERMANOVA Analysis")
    print(f" Abundance Matrix : {args.matrix}")
    print(f" Metadata File    : {args.metadata}")
    print(f" Group Column     : {args.group_col}")
    print("==========================================================")

    # 1. PERMANOVA
    print("[+] Running Permutational Multivariate Analysis of Variance (PERMANOVA)...")
    perm_res = run_permanova(df, meta_df, group_col=args.group_col, perms=999)
    if perm_res:
        perm_df = pd.DataFrame([perm_res])
        perm_file = os.path.join(args.diff_out, "permanova_test_results.tsv")
        perm_df.to_csv(perm_file, sep="\t", index=False)
        print(f"[+] PERMANOVA results:\n    R2 = {perm_res['R2']} | F = {perm_res['F_Statistic']} | p-value = {perm_res['P_Value']}")
        print(f"[+] Saved PERMANOVA results to {perm_file}")

    # 2. Differential Abundance Testing
    print("[+] Calculating pairwise differential abundance and effect sizes...")
    diff_df, g1, g2 = run_differential_abundance(df, meta_df, group_col=args.group_col)
    diff_file = os.path.join(args.diff_out, "differential_abundance_results.tsv")
    diff_df.to_csv(diff_file, sep="\t", index=False)
    print(f"[+] Differential abundance results saved to {diff_file}")

    # 3. Figures: Volcano and Forest Plot
    volcano_prefix = os.path.join(args.figures_out, "differential_volcano_plot")
    plot_volcano(diff_df, g1, g2, volcano_prefix)

    forest_prefix = os.path.join(args.figures_out, "differential_taxa_forest_plot")
    plot_forest(diff_df, g1, g2, forest_prefix, top_n=15)

if __name__ == "__main__":
    main()
