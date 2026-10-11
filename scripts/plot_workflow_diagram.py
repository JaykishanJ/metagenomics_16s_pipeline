#!/usr/bin/env python3
"""
plot_workflow_diagram.py
Generates a publication-grade vector (PDF) and 300 DPI raster (PNG) workflow diagram
for the Turnkey 16S rRNA Metagenomics Pipeline.
Suitable for import into Adobe Illustrator, Inkscape, or direct publication.
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

def create_workflow_figure(output_png, output_pdf):
    # Set publication figure canvas (16:10 aspect ratio, widescreen publication layout)
    fig, ax = plt.subplots(figsize=(18, 11), dpi=300)
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 11)
    ax.axis('off')

    # Background canvas
    fig.patch.set_facecolor('#F8FAFC')
    ax.set_facecolor('#F8FAFC')

    # Color Palette (Nature/Cell inspired, colorblind accessible)
    PALETTE = {
        'input_bg': '#EFF6FF',      'input_border': '#3B82F6',    'input_title': '#1E40AF',
        'db_bg': '#F1F5F9',         'db_border': '#64748B',       'db_title': '#334155',
        'prof_bg': '#ECFDF5',       'prof_border': '#10B981',     'prof_title': '#065F46',
        'krona_bg': '#F0FDFA',      'krona_border': '#14B8A6',    'krona_title': '#0F766E',
        'div_bg': '#FAF5FF',        'div_border': '#A855F7',      'div_title': '#6B21A8',
        'diff_bg': '#FFF1F2',       'diff_border': '#F43F5E',     'diff_title': '#9F1239',
        'func_bg': '#FFFBEB',       'func_border': '#F59E0B',     'func_title': '#92400E',
        'card_bg': '#FFFFFF',       'card_border': '#CBD5E1',     'text_dark': '#0F172A',
        'text_muted': '#475569',    'arrow': '#475569'
    }

    # Helper: Draw Subgraph Container
    def draw_subgraph(x, y, w, h, bg_col, border_col, title, subtitle=""):
        bg = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15,rounding_size=0.2",
                            facecolor=bg_col, edgecolor=border_col, linewidth=1.5, zorder=1)
        ax.add_patch(bg)
        # Subgraph Header (well above cards)
        ax.text(x + 0.25, y + h - 0.35, title, fontsize=11.5, fontweight='bold', color=border_col, zorder=3)
        if subtitle:
            ax.text(x + 0.25, y + h - 0.62, subtitle, fontsize=8.2, fontstyle='italic', color=PALETTE['text_muted'], zorder=3)

    # Helper: Draw Process Card
    def draw_card(x, y, w, h, title, tool_text, details, border_col=PALETTE['card_border'], header_col=PALETTE['text_dark']):
        card = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1,rounding_size=0.15",
                              facecolor=PALETTE['card_bg'], edgecolor=border_col, linewidth=1.1, zorder=2)
        ax.add_patch(card)
        # Title & Tool badge
        ax.text(x + 0.20, y + h - 0.28, title, fontsize=9.2, fontweight='bold', color=header_col, zorder=4)
        if tool_text:
            ax.text(x + 0.20, y + h - 0.52, tool_text, fontsize=7.8, fontweight='semibold', color=border_col, zorder=4)
        # Detail lines
        start_y = y + h - 0.78
        for i, line in enumerate(details):
            ax.text(x + 0.20, start_y - (i * 0.25), line, fontsize=7.3, color=PALETTE['text_muted'], zorder=4)

    # Helper: Draw Flow Arrow
    def draw_arrow(x1, y1, x2, y2, label="", dashed=False, color=PALETTE['arrow'], connectionstyle="arc3,rad=0", rad_label_offset=(0, 0.08)):
        ls = '--' if dashed else '-'
        arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                connectionstyle=connectionstyle,
                                arrowstyle="-|>",
                                linestyle=ls,
                                linewidth=1.4,
                                mutation_scale=12,
                                color=color,
                                zorder=5)
        ax.add_patch(arrow)
        if label:
            mx = (x1 + x2)/2 + rad_label_offset[0]
            my = (y1 + y2)/2 + rad_label_offset[1]
            ax.text(mx, my, label, fontsize=7.2, ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.18', facecolor='#FFFFFF', edgecolor=color, linewidth=0.6, alpha=0.95),
                    color=color, fontweight='semibold', zorder=6)

    # =========================================================================
    # Header Banner
    # =========================================================================
    ax.text(0.6, 10.45, "Turnkey 16S rRNA Metagenomics Analysis Pipeline", fontsize=17, fontweight='bold', color='#0F172A')
    ax.text(0.6, 10.12, "End-to-End Workflow: Quality Control  ->  Taxonomic Profiling  ->  Diversity Statistics  ->  Functional Annotation",
            fontsize=9.8, color='#475569')

    # =========================================================================
    # 1. Inputs Subgraph (Left top: x=0.6, w=3.4)
    # =========================================================================
    draw_subgraph(0.6, 5.5, 3.4, 4.3, PALETTE['input_bg'], PALETTE['input_border'], "1. Raw Sequencing Inputs", "Cohort: BioProject PRJNA1088906")
    draw_card(0.8, 7.6, 3.0, 1.4, "Paired-End FASTQ Reads", "Illumina NovaSeq / MiSeq",
              ["* File: data/fastq/*_1.fastq.gz, *_2.fastq.gz", "* 17 samples (72k-184k reads/sample)", "* 16S rRNA V3-V4 hypervariable region"],
              border_col=PALETTE['input_border'], header_col=PALETTE['input_title'])
    draw_card(0.8, 5.7, 3.0, 1.6, "Clinical Metadata", "Patient Cohort Stratification",
              ["* File: metadata.tsv", "* Cohorts: Healthy (n=4) vs. UC (n=13)", "* Variables: SampleID, Group, SRA Accession"],
              border_col=PALETTE['input_border'], header_col=PALETTE['input_title'])

    # =========================================================================
    # Reference DB & Indexing (Left bottom: x=0.6, w=3.4)
    # =========================================================================
    draw_subgraph(0.6, 0.6, 3.4, 4.6, PALETTE['db_bg'], PALETTE['db_border'], "Reference DB & Indexing", "01_setup.sh Automated Build")
    draw_card(0.8, 2.75, 3.0, 1.5, "Greengenes 16S Database", "kraken2-build --special greengenes",
              ["* Reference: GG 16S rRNA sequences", "* Lineage: NCBI + Greengenes taxonomy", "* Files: hash.k2d, opts.k2d, taxo.k2d"],
              border_col=PALETTE['db_border'], header_col=PALETTE['db_title'])
    draw_card(0.8, 0.8, 3.0, 1.65, "Bracken Read-Length Index", "bracken-build -l 150 -t 8",
              ["* Parameter: Read length = 150 bp", "* Precomputed k-mer distribution", "* File: database150mers.kmer_distrib"],
              border_col=PALETTE['db_border'], header_col=PALETTE['db_title'])

    # =========================================================================
    # 2. Profiling & Abundance Subgraph (Column 2: x=4.4, w=4.0)
    # =========================================================================
    draw_subgraph(4.4, 0.6, 4.0, 9.2, PALETTE['prof_bg'], PALETTE['prof_border'], "2. Taxonomic Profiling & Matrix", "02_classify_pipeline.sh & 03_downstream_analysis.sh")
    draw_card(4.6, 7.3, 3.6, 1.5, "Kraken2 Exact k-mer Matching", "kraken2 --confidence 0.0 --paired",
              ["* Algorithm: Lowest Common Ancestor (LCA)", "* Multi-threaded execution (GNU Parallel)", "* Output: *.kraken & *.report (98.4-99.9% classified)"],
              border_col=PALETTE['prof_border'], header_col=PALETTE['prof_title'])
    draw_card(4.6, 5.25, 3.6, 1.65, "Bracken Bayesian Redistribution", "bracken -r 150 -l S -t 4",
              ["* Algorithm: Bayesian probability re-estimation", "* Resolves ambiguous genus/family reads to species", "* Output: *.bracken abundance files"],
              border_col=PALETTE['prof_border'], header_col=PALETTE['prof_title'])
    draw_card(4.6, 3.15, 3.6, 1.7, "Merged Species Abundance Matrix", "03_downstream_analysis.sh (pandas outer join)",
              ["* File: results/taxonomy/species_abundance_matrix.txt", "* Dimensions: 290 species x 17 clinical samples", "* Relative abundance fractions (sum = 1.0)"],
              border_col=PALETTE['prof_border'], header_col=PALETTE['prof_title'])
    draw_card(4.6, 0.8, 3.6, 1.95, "Sequencing QC & Rarefaction", "scripts/qc_depth_analysis.py",
              ["* Reads/sample depth boxplot (median: 117k reads)", "* Classification efficiency distribution (>98.4%)", "* In-silico rarefaction curves (Shannon & Richness)", "* Confirms sampling saturation plateau"],
              border_col=PALETTE['db_border'], header_col=PALETTE['db_title'])

    # =========================================================================
    # 3. Interactive Krona Visualization (Column 3 Top: x=8.8, w=4.1)
    # =========================================================================
    draw_subgraph(8.8, 7.5, 4.1, 2.3, PALETTE['krona_bg'], PALETTE['krona_border'], "3. Interactive Exploration", "03_krona_visualization.sh")
    draw_card(9.0, 7.65, 3.7, 1.25, "Radial Multi-Rank Krona Charts", "ktImportText (Hierarchical Pie)",
              ["* Interactive HTML: results/krona/all_samples_krona.html", "* Full multi-level hierarchy: Domain -> Phylum -> Species"],
              border_col=PALETTE['krona_border'], header_col=PALETTE['krona_title'])

    # =========================================================================
    # 4. Diversity & Statistical Suite (Column 3 Bottom: x=8.8, w=4.1)
    # =========================================================================
    draw_subgraph(8.8, 0.6, 4.1, 6.7, PALETTE['div_bg'], PALETTE['div_border'], "4. Diversity & Statistical Suite", "04_diversity_analysis.sh & statistical scripts")
    draw_card(9.0, 5.0, 3.7, 1.7, "Multi-Level Taxonomic Profiles", "scripts/plot_multilevel_taxa.py",
              ["* Stacked barplots: Phylum & Genus compositions", "* Group-mean comparative composition (Healthy vs UC)", "* Top-25 dominant species hierarchical clustered heatmap"],
              border_col=PALETTE['div_border'], header_col=PALETTE['div_title'])
    draw_card(9.0, 3.0, 3.7, 1.7, "Alpha & Beta Diversity Suite", "scripts/diversity_analysis.py",
              ["* Alpha diversity: Shannon, Simpson, Observed, Evenness", "* Beta diversity: Bray-Curtis dissimilarity distance matrix", "* PCoA 2D Ordination (Axis 1: 28.4%, Axis 2: 19.8%)"],
              border_col=PALETTE['div_border'], header_col=PALETTE['div_title'])
    draw_card(9.0, 0.8, 3.7, 1.9, "Hypothesis Testing & Networks", "differential_abundance.py & network_analysis.py",
              ["* PERMANOVA: R2 = 0.156, F = 2.772, p = 0.016*", "* Biomarker Volcano & Forest plots (BH-FDR < 0.1)", "* Microbial co-occurrence network (Spearman |r| >= 0.6)"],
              border_col=PALETTE['diff_border'], header_col=PALETTE['diff_title'])

    # =========================================================================
    # 5. Functional Annotation Subgraph (Column 4: x=13.3, w=4.2)
    # =========================================================================
    draw_subgraph(13.3, 0.6, 4.2, 9.2, PALETTE['func_bg'], PALETTE['func_border'], "5. Functional Annotation", "05_functional_picrust2.sh & profiling.py")
    draw_card(13.5, 7.0, 3.8, 1.85, "OTU Clustering & Centroids", "prepare_picrust2.py (VSEARCH)",
              ["* Dereplication & sequence clustering at 97% identity", "* Centroid representative fasta: rep_seqs.fasta", "* Sample OTU abundance table: otu_table.tsv", "* Formats amplicons for phylogenetic placement"],
              border_col=PALETTE['func_border'], header_col=PALETTE['func_title'])
    draw_card(13.5, 4.2, 3.8, 2.5, "PICRUSt2 Metabolic Inference", "picrust2_pipeline.py (MetaCyc/KEGG)",
              ["* Phylogenetic hidden state prediction (hsp.py)", "* Enzyme Commission (EC) numbers predicted", "* KEGG Orthology (KO) functional gene families", "* MetaCyc pathway abundance: path_abun_unstrat.tsv", "* Deliverables: Top-20 pathway barplots & heatmaps"],
              border_col=PALETTE['func_border'], header_col=PALETTE['func_title'])
    draw_card(13.5, 0.8, 3.8, 2.9, "Curated Microbial Pathways", "scripts/functional_profiling.py",
              ["* Short-Chain Fatty Acid (SCFA) metabolism:", "  - Acetate synthesis (EC 2.3.1.8 / ackA)", "  - Propionate production (acrylate pathway)", "  - Butyrate synthesis (buk / but pathway)", "* Mucus degradation (Akkermansia muciniphila signature)", "* Lactate fermentation & primary carbon cycling", "* Deliverable: Curated functional pathway comparative heatmap"],
              border_col=PALETTE['func_border'], header_col=PALETTE['func_title'])

    # =========================================================================
    # Connecting Arrows & Clean Data Flows
    # =========================================================================
    # 1. Primary Sequencing Input Flow: FASTQ -> Kraken2
    draw_arrow(3.8, 8.3, 4.6, 8.3, label="FASTQ")

    # 2. Reference DB Indexing -> Profiling (short clean arrows)
    draw_arrow(3.8, 3.7, 4.6, 7.6, dashed=True, label="K2 DB", connectionstyle="arc3,rad=-0.12")
    draw_arrow(3.8, 1.65, 4.6, 5.8, dashed=True, label="150mer Dist", connectionstyle="arc3,rad=-0.10")

    # 3. Column 2 Vertical Cascade (Kraken2 -> Bracken -> Matrix -> QC)
    draw_arrow(6.4, 7.3, 6.4, 6.9, label="*.report")
    draw_arrow(6.4, 5.25, 6.4, 4.85, label="*.bracken")
    draw_arrow(6.4, 3.15, 6.4, 2.75, label="Matrix QC")

    # 4. Profiling -> Krona Interactive Exploration
    draw_arrow(8.2, 8.0, 9.0, 8.25, label="Report Tree", connectionstyle="arc3,rad=0.08")

    # 5. Species Abundance Matrix -> Diversity & Statistical Suite (across empty gap x=8.2 to 9.0)
    draw_arrow(8.2, 4.15, 9.0, 5.85, label="Taxa Matrix", connectionstyle="arc3,rad=0.05")
    draw_arrow(8.2, 3.95, 9.0, 3.85, label="Abundance TSV")
    draw_arrow(8.2, 3.75, 9.0, 1.85, label="Species Counts", connectionstyle="arc3,rad=-0.05")

    # 6. Clinical Metadata -> Diversity & Statistics (Routed via inter-column gap x=4.2 and bottom corridor y=0.25)
    meta_path = patches.Path([(3.8, 6.4), (4.2, 6.4), (4.2, 0.25), (10.85, 0.25), (10.85, 0.8)],
                             [patches.Path.MOVETO, patches.Path.LINETO, patches.Path.LINETO, patches.Path.LINETO, patches.Path.LINETO])
    meta_patch = patches.PathPatch(meta_path, edgecolor='#9333EA', facecolor='none',
                                   linestyle='--', linewidth=1.4, zorder=5)
    ax.add_patch(meta_patch)
    meta_head = FancyArrowPatch((10.85, 0.65), (10.85, 0.8), arrowstyle="-|>", color='#9333EA', mutation_scale=12, zorder=6)
    ax.add_patch(meta_head)
    ax.text(7.5, 0.25, "Clinical Metadata (Healthy vs. UC Cohorts)", fontsize=7.2, fontweight='bold',
            color='#9333EA', ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.18', facecolor='#FFFFFF', edgecolor='#9333EA', linewidth=0.6), zorder=7)

    # 7. Raw FastQ -> PICRUSt2 VSEARCH OTU Clustering (Routed via inter-column gap x=4.2 and top corridor y=10.0)
    fastq_path = patches.Path([(3.8, 8.5), (4.2, 8.5), (4.2, 10.0), (15.4, 10.0), (15.4, 8.85)],
                              [patches.Path.MOVETO, patches.Path.LINETO, patches.Path.LINETO, patches.Path.LINETO, patches.Path.LINETO])
    fastq_patch = patches.PathPatch(fastq_path, edgecolor='#D97706', facecolor='none',
                                    linestyle='-', linewidth=1.4, zorder=5)
    ax.add_patch(fastq_patch)
    fastq_head = FancyArrowPatch((15.4, 8.95), (15.4, 8.85), arrowstyle="-|>", color='#D97706', mutation_scale=12, zorder=6)
    ax.add_patch(fastq_head)
    ax.text(9.8, 10.0, "Raw FastQ Amplicon Sequences (VSEARCH OTU Clustering)", fontsize=7.2, fontweight='bold',
            color='#D97706', ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.18', facecolor='#FFFFFF', edgecolor='#D97706', linewidth=0.6), zorder=7)

    # 8. Within Functional Annotation Cascade
    draw_arrow(15.4, 7.0, 15.4, 6.7, label="rep_seqs & otu_table", color='#D97706')
    draw_arrow(15.4, 4.2, 15.4, 3.7, label="MetaCyc & Pathways", color='#D97706')

    # 9. Species Abundance Matrix -> Curated Microbial Pathways (clean arc across bottom)
    draw_arrow(8.2, 3.4, 13.5, 1.8, dashed=True, label="Species Abundance Matrix",
               color='#D97706', connectionstyle="arc3,rad=-0.08", rad_label_offset=(0, 0.12))



    # Save to disk
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.savefig(output_pdf, format='pdf', bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    print(f"[+] Saved workflow overview to:\n  - {output_png}\n  - {output_pdf}")

if __name__ == '__main__':
    png_path = "results/figures/workflow_pipeline_overview.png"
    pdf_path = "results/figures/workflow_pipeline_overview.pdf"
    create_workflow_figure(png_path, pdf_path)
