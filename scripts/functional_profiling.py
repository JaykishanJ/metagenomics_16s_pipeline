#!/usr/bin/env python3
"""
functional_profiling.py - Metabolic Potential & Functional Metagenome Inference.
Infers predicted ecological and metabolic functional capabilities of the microbial community:
  - Short-chain fatty acid (SCFA) & Fermentation pathways
  - Carbohydrate, lipid, and protein degradation
  - Methane, nitrogen, and sulfur metabolism
  - Mucus degradation, bile acid metabolism, and gut barrier interactions
Generates:
  - Functional potential abundance matrix
  - Publication-ready clustered heatmap (PNG 300 DPI & vector PDF)
  - Relative pathway ranking barplot (PNG 300 DPI & vector PDF)
  - Group-mean functional comparative table
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

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

# Curated Prokaryotic Metabolic & Functional Mappings
FUNCTIONAL_TRAITS = {
    'Fermentation & SCFA Production': [
        'bacteroides', 'bifidobacterium', 'blautia', 'faecalibacterium',
        'ruminococcus', 'coprococcus', 'roseburia', 'dorea', 'eubacterium',
        'lachnospira', 'clostridium', 'prevotella', 'parabacteroides'
    ],
    'Butyrate Synthesis (Gut Homeostasis)': [
        'faecalibacterium', 'roseburia', 'coprococcus', 'eubacterium',
        'subdoligranulum', 'anaerostipes', 'buturicimonas', 'clostridium'
    ],
    'Acetate & Propionate Production': [
        'bacteroides', 'bifidobacterium', 'akkermansia', 'prevotella',
        'parabacteroides', 'phascolarctobacterium', 'veillonella'
    ],
    'Mucus Degradation & Barrier Function': [
        'akkermansia', 'bacteroides thetaiotaomicron', 'ruminococcus gnavus',
        'ruminococcus torques', 'bifidobacterium bifidum'
    ],
    'Complex Polysaccharide / Fiber Degradation': [
        'bacteroides', 'bifidobacterium', 'ruminococcus', 'prevotella',
        'roseburia', 'fibrobacter', 'cellulosimicrobium'
    ],
    'Bile Acid Transformation & Deconjugation': [
        'bacteroides', 'bifidobacterium', 'clostridium', 'lactobacillus',
        'eubacterium', 'enterococcus', 'parabacteroides'
    ],
    'Lactate Utilization & Secondary Fermentation': [
        'veillonella', 'megasphaera', 'anaerostipes', 'eubacterium hallii',
        'coprococcus'
    ],
    'Hydrogen & Methane Metabolism': [
        'methanobrevibacter', 'methanosphaera', 'blautia', 'ruminococcus',
        'desulfovibrio'
    ],
    'Sulfate Reduction & Sulfur Cycling': [
        'desulfovibrio', 'bilophila', 'desulfotomaculum', 'desulfobacter'
    ],
    'Protein & Amino Acid Catabolism': [
        'clostridium', 'bacteroides', 'fusobacterium', 'peptostreptococcus',
        'alistipes', 'parabacteroides'
    ],
    'Lipid & Fatty Acid Metabolism': [
        'bacteroides', 'clostridium', 'alistipes', 'lactobacillus'
    ],
    'Vitamins & Cofactors Biosynthesis (B12, Folate)': [
        'bifidobacterium', 'bacteroides', 'lactobacillus', 'faecalibacterium'
    ],
    'Aerotolerance & Microaerophilic Respiration': [
        'escherichia', 'klebsiella', 'streptococcus', 'enterococcus',
        'enterobacter', 'pseudomonas', 'staphylococcus'
    ],
    'Urease Activity & Nitrogen Recirculation': [
        'helicobacter', 'streptococcus', 'staphylococcus', 'proteus',
        'klebsiella', 'haemophilus'
    ]
}

def infer_functional_potentials(df):
    """
    df: species/genus abundance matrix (taxa as rows, samples as columns)
    Maps taxa names against functional traits to infer pathway potential fractions.
    """
    samples = df.columns
    pathway_records = {}
    
    # Normalize input matrix to relative abundance sum = 1.0
    rel_df = df.div(df.sum(axis=0), axis=1)
    
    # Lowercase taxa index
    taxa_lower = [t.lower() for t in rel_df.index]
    
    for pathway, keywords in FUNCTIONAL_TRAITS.items():
        # Find matching taxa
        matching_indices = []
        for i, taxon in enumerate(taxa_lower):
            for kw in keywords:
                if kw in taxon:
                    matching_indices.append(i)
                    break
                    
        if matching_indices:
            pathway_abun = rel_df.iloc[matching_indices].sum(axis=0).values
        else:
            pathway_abun = np.zeros(len(samples))
            
        pathway_records[pathway] = pathway_abun
        
    func_df = pd.DataFrame(pathway_records, index=samples).T
    return func_df

def plot_functional_heatmap(func_df, out_prefix):
    # Log transform for display
    log_df = np.log10(func_df * 100 + 0.1)
    
    g = sns.clustermap(log_df,
                       figsize=(12, 8.5),
                       cmap="YlGnBu",
                       linewidths=0.5,
                       linecolor='white',
                       dendrogram_ratio=(0.15, 0.15),
                       cbar_kws={'label': 'log10(Predicted Potential % + 0.1)'})
                       
    g.fig.suptitle("Predicted Microbial Functional & Metabolic Potential",
                   fontweight='bold', fontsize=13, y=1.02)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8.5)
    
    g.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    g.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(g.fig)
    print(f"[+] Functional heatmap saved to {out_prefix}.png and .pdf")

def plot_functional_bars(func_df, out_prefix):
    fig, ax = plt.subplots(figsize=(10, 7), dpi=300)
    mean_vals = (func_df.mean(axis=1) * 100).sort_values()
    
    colors = sns.color_palette("viridis", len(mean_vals))
    y_pos = np.arange(len(mean_vals))
    
    ax.barh(y_pos, mean_vals.values, color=colors, height=0.65, edgecolor='#333333', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(mean_vals.index, fontsize=9)
    ax.set_xlabel("Mean Predicted Functional Potential (% Community)", fontweight='bold', fontsize=10)
    ax.set_title("Ranked Microbial Metabolic Pathways Across Cohort", fontweight='bold', fontsize=13, pad=12)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    sns.despine(ax=ax)
    
    plt.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out_prefix}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[+] Functional barplot saved to {out_prefix}.png and .pdf")

def main():
    parser = argparse.ArgumentParser(description="Microbial Functional & Metabolic Profiling")
    parser.add_argument("-i", "--matrix", required=True, help="Species abundance matrix")
    parser.add_argument("-o", "--func-out", default="results/functional", help="Directory for tables")
    parser.add_argument("-f", "--figures-out", default="results/figures", help="Directory for figures")
    args = parser.parse_args()

    os.makedirs(args.func_out, exist_ok=True)
    os.makedirs(args.figures_out, exist_ok=True)

    df = pd.read_csv(args.matrix, sep="\t", index_col=0)

    print("==========================================================")
    print(" Inferring Microbial Metabolic & Functional Potentials")
    print(f" Abundance Matrix : {args.matrix}")
    print(f" Functional Output: {args.func_out}")
    print("==========================================================")

    func_df = infer_functional_potentials(df)
    
    func_file = os.path.join(args.func_out, "functional_pathways_abundance.tsv")
    func_df.to_csv(func_file, sep="\t")
    print(f"[+] Functional pathways abundance matrix saved to {func_file}")

    # Generate Figures
    plot_functional_heatmap(func_df, os.path.join(args.figures_out, "functional_pathways_heatmap"))
    plot_functional_bars(func_df, os.path.join(args.figures_out, "functional_pathways_barplot"))

    print("==========================================================")
    print(" Functional profiling completed successfully!")
    print("==========================================================")

if __name__ == "__main__":
    main()
