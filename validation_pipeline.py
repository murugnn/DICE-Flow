#!/usr/bin/env python3
"""
DiCE-Duan Validation & Benchmarking Pipeline
=============================================
Publication-quality validation of disease gene predictions.

Usage:
    python validation_pipeline.py --all
    python validation_pipeline.py --analysis 1 2 3
    python validation_pipeline.py --all --skip-pubmed
"""

import argparse
import os
import sys
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import hypergeom, mannwhitneyu, ttest_ind

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ─── Style Configuration ────────────────────────────────────────────────────
STYLE = {
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (8, 6),
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
}
plt.rcParams.update(STYLE)
sns.set_theme(style='whitegrid', palette='deep')
PALETTE_MAIN = ['#2563EB', '#DC2626', '#16A34A', '#F59E0B', '#7C3AED']


# ═══════════════════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════

def find_project_root():
    """Find the DICE project root."""
    script_dir = Path(__file__).resolve().parent
    for d in [script_dir, script_dir.parent]:
        if (d / 'data' / 'results').exists():
            return d
    return script_dir


def load_dice_genes(root):
    """Load ranked DiCE genes."""
    path = root / 'data' / 'results' / 'dice_genes_top500.txt'
    df = pd.read_csv(path, sep='\t')
    log.info(f"Loaded {len(df)} DiCE genes from {path.name}")
    return df


def load_centrality(root, condition):
    """Load centrality file (normal or tumor)."""
    path = root / 'data' / 'results' / f'centrality_{condition}.txt'
    df = pd.read_csv(path, sep='\t')
    log.info(f"Loaded {len(df)} genes from centrality_{condition}.txt")
    return df


def load_knockout(root):
    """Load knockout analysis results."""
    path = root / 'data' / 'results' / 'knockout_analysis.txt'
    df = pd.read_csv(path, sep='\t')
    log.info(f"Loaded {len(df)} knockout results")
    return df


def load_expression(root):
    """Load normal and tumor expression matrices."""
    normal_path = root / 'data' / 'raw' / 'expression' / 'expression_normal.txt'
    tumor_path = root / 'data' / 'raw' / 'expression' / 'expression_tumor.txt'

    normal = pd.read_csv(normal_path, sep='\t', index_col=0)
    tumor = pd.read_csv(tumor_path, sep='\t', index_col=0)
    log.info(f"Loaded expression: normal={normal.shape}, tumor={tumor.shape}")
    return normal, tumor


def load_known_genes(root):
    """Load known disease gene lists."""
    known_dir = root / 'data' / 'known_genes'
    result = {}
    for name, fname in [('COSMIC', 'cosmic_census.txt'), ('DisGeNET', 'disgenet_prostate.txt')]:
        path = known_dir / fname
        if path.exists():
            genes = set(l.strip() for l in path.read_text().splitlines() if l.strip() and not l.startswith('#'))
            result[name] = genes
            log.info(f"Loaded {len(genes)} known genes from {fname}")
        else:
            log.warning(f"Known gene file not found: {path}")
    return result


def load_network_edges(root):
    """Load merged network edges."""
    path = root / 'data' / 'processed' / 'networks' / 'merged_network.txt'
    if not path.exists():
        return None
    edges = pd.read_csv(path, sep='\t', header=None, names=['gene1', 'gene2'])
    log.info(f"Loaded {len(edges)} network edges")
    return edges


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 1: PATHWAY ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════

def analysis_pathway_enrichment(dice_df, root, out_fig, out_tbl):
    """Enrichment analysis of top DiCE genes split by Upregulated vs Downregulated."""
    log.info("=" * 60)
    log.info("Analysis 1: Pathway Enrichment (Split Up/Down)")
    log.info("=" * 60)

    try:
        import gseapy
    except ImportError:
        log.error("gseapy not installed. Skipping pathway enrichment.")
        return None

    # Load DEA log2FC data to split genes
    dea_table_path = root / 'data' / 'results' / 'dea_results_table.txt'
    if not dea_table_path.exists():
        log.warning(f"DEA table not found at {dea_table_path}. Cannot perform split enrichment.")
        return None
        
    dea_df = pd.read_csv(dea_table_path, sep='\t')
    log2fc_map = dict(zip(dea_df['gene'], dea_df['log2FC']))
    
    # 1. Split top 150 DiCE genes into Up and Down
    top_genes = dice_df['gene'].head(150).tolist()
    
    up_genes = [g for g in top_genes if log2fc_map.get(g, 0) > 0]
    down_genes = [g for g in top_genes if log2fc_map.get(g, 0) < 0]
    
    log.info(f"Split {len(top_genes)} top genes into {len(up_genes)} Upregulated and {len(down_genes)} Downregulated.")

    gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021']
    
    # Process both lists
    def run_enrichr(genes, label):
        if len(genes) < 5:
            log.warning(f"Too few {label} genes for enrichr.")
            return None
        all_results = []
        for lib in gene_sets:
            try:
                enr = gseapy.enrichr(gene_list=genes, gene_sets=lib, organism='human', outdir=None, cutoff=0.05)
                res = enr.results
                res['Database'] = lib.split('_')[0]
                res['Regulation'] = label
                all_results.append(res)
                log.info(f"  {label} {lib}: {len(res[res['Adjusted P-value'] < 0.05])} significant terms")
            except Exception as e:
                log.warning(f"  Enrichr failed for {label} {lib}: {e}")
                
        if not all_results: return None
        return pd.concat(all_results, ignore_index=True)

    up_results = run_enrichr(up_genes, 'Upregulated')
    down_results = run_enrichr(down_genes, 'Downregulated')
    
    if up_results is None and down_results is None:
        return None
        
    combined = pd.concat([df for df in [up_results, down_results] if df is not None], ignore_index=True)
    sig = combined[combined['Adjusted P-value'] < 0.05].sort_values(['Regulation', 'Adjusted P-value'])
    sig.to_csv(out_tbl / '01_pathway_enrichment_split.csv', index=False)
    log.info(f"  Saved {len(sig)} significant terms to table")

    # ── Divergent Bar plot ──
    try:
        plot_df = []
        if up_results is not None:
            up_sig = up_results[up_results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value').head(10)
            plot_df.append(up_sig)
        if down_results is not None:
            down_sig = down_results[down_results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value').head(10)
            plot_df.append(down_sig)
            
        if not plot_df:
            return combined
            
        plot_df = pd.concat(plot_df, ignore_index=True)
        # Flip values for downregulated for divergent plot
        plot_df['neg_log10_padj'] = -np.log10(plot_df['Adjusted P-value'].clip(lower=1e-50))
        plot_df.loc[plot_df['Regulation'] == 'Downregulated', 'neg_log10_padj'] *= -1
        
        # Sort so Up are on top, Down on bottom
        plot_df = plot_df.sort_values('neg_log10_padj')
        plot_df['Term_short'] = plot_df['Term'].str[:50]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        colors = ['#EF4444' if r == 'Upregulated' else '#3B82F6' for r in plot_df['Regulation']]
        
        bars = ax.barh(range(len(plot_df)), plot_df['neg_log10_padj'], color=colors)
        
        # Add labels dynamically inside/outside bars
        for i, (val, term) in enumerate(zip(plot_df['neg_log10_padj'], plot_df['Term_short'])):
            if val > 0:
                ax.text(0.1, i, term, va='center', ha='left', fontsize=9, color='white', fontweight='bold')
            else:
                ax.text(-0.1, i, term, va='center', ha='right', fontsize=9, color='white', fontweight='bold')
                
        ax.set_yticks([]) # Hide traditional y-ticks
        
        # Format axes
        ax.axvline(0, color='black', linewidth=1)
        ax.set_xlabel('Significance ($-\\log_{10}$(Adjusted p-value))')
        # Labels for axes
        ticks = ax.get_xticks()
        ax.set_xticklabels([str(abs(int(t))) for t in ticks])
        
        ax.set_title('Pathway Enrichment (Upregulated DiCE vs Downregulated DiCE)')
        
        from matplotlib.patches import Patch
        legend_elems = [Patch(color='#EF4444', label='Upregulated (Log2FC > 0)'), 
                        Patch(color='#3B82F6', label='Downregulated (Log2FC < 0)')]
        ax.legend(handles=legend_elems, loc='upper left')
        
        plt.tight_layout()
        fig.savefig(out_fig / '01_pathway_enrichment_divergent_barplot.png')
        plt.close(fig)
        log.info("  Saved divergent barplot")
    except Exception as e:
        log.warning(f"  Failed generating divergent plot: {e}")

    return combined


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 2: HYPERGEOMETRIC SIGNIFICANCE TEST
# ═══════════════════════════════════════════════════════════════════════════

def analysis_hypergeometric(dice_df, known_genes, universe_size, out_fig, out_tbl):
    """Hypergeometric test for overlap with known disease genes."""
    log.info("=" * 60)
    log.info("Analysis 2: Hypergeometric Significance Test")
    log.info("=" * 60)

    predicted = set(dice_df['gene'].head(100).tolist())
    results = []

    for db_name, db_genes in known_genes.items():
        overlap = predicted & db_genes
        k = len(overlap)  # successes in sample
        M = universe_size  # population size
        n = len(db_genes & set(dice_df['gene']))  # total successes in population
        N = len(predicted)  # sample size

        # P(X >= k) — survival function
        pval = hypergeom.sf(k - 1, M, n, N)
        expected = N * n / M

        fold_enrichment = (k / expected) if expected > 0 else float('inf')

        results.append({
            'Database': db_name,
            'Known_genes_in_universe': n,
            'Predicted_top100': N,
            'Overlap': k,
            'Expected_overlap': round(expected, 2),
            'Fold_enrichment': round(fold_enrichment, 2),
            'P_value': pval,
            'Significant': pval < 0.05,
            'Overlapping_genes': ', '.join(sorted(overlap))
        })
        log.info(f"  {db_name}: overlap={k}, expected={expected:.2f}, fold={fold_enrichment:.2f}, p={pval:.2e}")

    res_df = pd.DataFrame(results)
    res_df.to_csv(out_tbl / '02_hypergeometric_results.csv', index=False)

    # ── Plot ──
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for i, (_, row) in enumerate(res_df.iterrows()):
        ax = axes[i]
        vals = [row['Expected_overlap'], row['Overlap']]
        colors_bar = ['#94A3B8', PALETTE_MAIN[0] if row['Significant'] else '#EF4444']
        bars = ax.bar(['Expected', 'Observed'], vals, color=colors_bar, width=0.5, edgecolor='black', linewidth=0.8)
        ax.set_title(f"{row['Database']}\n(p = {row['P_value']:.2e})", fontweight='bold')
        ax.set_ylabel('Gene Overlap Count')

        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    f'{val:.1f}' if isinstance(val, float) else str(val),
                    ha='center', va='bottom', fontweight='bold')

        ax.set_ylim(0, max(vals) * 1.3 + 1)

    fig.suptitle('Overlap with Known Disease Genes (Top 100 DiCE Genes)', fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(out_fig / '02_hypergeometric_overlap.png')
    plt.close(fig)
    log.info("  Saved overlap plot")
    return res_df


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 3: PRECISION@K
# ═══════════════════════════════════════════════════════════════════════════

def analysis_precision_at_k(dice_df, known_genes, out_fig, out_tbl):
    """Precision@K for K = 10, 20, 50, 100."""
    log.info("=" * 60)
    log.info("Analysis 3: Precision@K Evaluation")
    log.info("=" * 60)

    K_values = [10, 20, 50, 100, 200, 500]
    all_known = set()
    for genes in known_genes.values():
        all_known |= genes

    records = []
    per_db = {db: [] for db in known_genes}
    per_db['Combined'] = []

    for k in K_values:
        top_k = set(dice_df['gene'].head(k))
        # Combined
        prec = len(top_k & all_known) / k
        per_db['Combined'].append(prec)
        row = {'K': k, 'Combined_precision': round(prec, 4)}

        for db_name, db_genes in known_genes.items():
            p = len(top_k & db_genes) / k
            per_db[db_name].append(p)
            row[f'{db_name}_precision'] = round(p, 4)

        records.append(row)
        log.info(f"  K={k}: combined precision = {prec:.4f}")

    res_df = pd.DataFrame(records)
    res_df.to_csv(out_tbl / '03_precision_at_k.csv', index=False)

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(8, 5.5))
    markers = ['o', 's', 'D']
    colors_pk = [PALETTE_MAIN[0], PALETTE_MAIN[1], PALETTE_MAIN[2]]

    for i, (name, vals) in enumerate(per_db.items()):
        ax.plot(K_values, vals, marker=markers[i % len(markers)], color=colors_pk[i % len(colors_pk)],
                label=name, linewidth=2, markersize=7)

    ax.set_xlabel('K (Number of Top Genes)')
    ax.set_ylabel('Precision@K')
    ax.set_title('Precision@K — DiCE-Duan Gene Predictions')
    ax.legend(frameon=True)
    ax.set_xticks(K_values)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_fig / '03_precision_at_k.png')
    plt.close(fig)
    log.info("  Saved Precision@K plot")
    return res_df


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 4: NETWORK TOPOLOGY VALIDATION
# ═══════════════════════════════════════════════════════════════════════════

def analysis_network_topology(dice_df, cent_normal, cent_tumor, out_fig, out_tbl):
    """Centrality shift analysis: scatter + boxplots."""
    log.info("=" * 60)
    log.info("Analysis 4: Network Topology Validation")
    log.info("=" * 60)

    # Merge centralities
    norm = cent_normal.rename(columns={'gene': 'gene'}).copy()
    tum = cent_tumor.rename(columns={'gene': 'gene'}).copy()

    merged = norm.merge(tum, on='gene', suffixes=('_normal', '_tumor'))

    # Compute delta centrality
    for metric in ['betweenness_normalized', 'eigenvector_normalized']:
        col_n = f'{metric}_normal'
        col_t = f'{metric}_tumor'
        if col_n in merged.columns and col_t in merged.columns:
            merged[f'delta_{metric}'] = merged[col_t] - merged[col_n]

    # If the expected columns don't exist, try alternative naming
    if 'delta_betweenness_normalized' not in merged.columns:
        # Try to find betweenness columns
        bet_cols_n = [c for c in merged.columns if 'betweenness' in c.lower() and 'normal' in c.lower()]
        bet_cols_t = [c for c in merged.columns if 'betweenness' in c.lower() and 'tumor' in c.lower()]
        if bet_cols_n and bet_cols_t:
            merged['delta_betweenness'] = merged[bet_cols_t[0]] - merged[bet_cols_n[0]]
        eig_cols_n = [c for c in merged.columns if 'eigenvector' in c.lower() and 'normal' in c.lower()]
        eig_cols_t = [c for c in merged.columns if 'eigenvector' in c.lower() and 'tumor' in c.lower()]
        if eig_cols_n and eig_cols_t:
            merged['delta_eigenvector'] = merged[eig_cols_t[0]] - merged[eig_cols_n[0]]

    # Mark candidates
    top_genes = set(dice_df['gene'].head(100))
    merged['is_candidate'] = merged['gene'].isin(top_genes)

    # Find delta columns
    delta_cols = [c for c in merged.columns if c.startswith('delta_')]
    if len(delta_cols) < 2:
        log.warning("  Could not find enough delta centrality columns.")
        # Use whatever delta columns exist or create placeholders
        if len(delta_cols) == 0:
            log.error("  No delta columns found. Skipping topology analysis.")
            return None

    merged.to_csv(out_tbl / '04_centrality_statistics.csv', index=False)

    # ── Scatter plot ──
    if len(delta_cols) >= 2:
        fig, ax = plt.subplots(figsize=(8, 7))
        bg = merged[~merged['is_candidate']]
        fg = merged[merged['is_candidate']]

        ax.scatter(bg[delta_cols[0]], bg[delta_cols[1]], alpha=0.15, s=8, c='#94A3B8', label='Background')
        ax.scatter(fg[delta_cols[0]], fg[delta_cols[1]], alpha=0.8, s=25, c='#DC2626',
                   edgecolors='black', linewidths=0.3, label='DiCE Top 100', zorder=5)

        ax.set_xlabel(f'Δ{delta_cols[0].replace("delta_", "")}')
        ax.set_ylabel(f'Δ{delta_cols[1].replace("delta_", "")}')
        ax.set_title('Differential Centrality: Tumor − Normal')
        ax.axhline(0, color='black', linewidth=0.5, alpha=0.5)
        ax.axvline(0, color='black', linewidth=0.5, alpha=0.5)
        ax.legend(frameon=True, loc='upper left')
        plt.tight_layout()
        fig.savefig(out_fig / '04_centrality_shift_scatter.png')
        plt.close(fig)
        log.info("  Saved centrality scatter")

    # ── Boxplot ──
    fig, axes = plt.subplots(1, len(delta_cols), figsize=(5 * len(delta_cols), 5))
    if len(delta_cols) == 1:
        axes = [axes]

    for i, dc in enumerate(delta_cols):
        data_bg = merged.loc[~merged['is_candidate'], dc].dropna()
        data_fg = merged.loc[merged['is_candidate'], dc].dropna()

        stat, pval = mannwhitneyu(data_fg, data_bg, alternative='greater')

        box_data = [data_bg.values, data_fg.values]
        bp = axes[i].boxplot(box_data, labels=['Background', 'DiCE\nTop 100'],
                             patch_artist=True, widths=0.5)
        bp['boxes'][0].set_facecolor('#E2E8F0')
        bp['boxes'][1].set_facecolor('#FCA5A5')
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(2)

        axes[i].set_ylabel(f'Δ{dc.replace("delta_", "")}')
        axes[i].set_title(f'p = {pval:.2e}', fontsize=10)
        axes[i].axhline(0, color='grey', linestyle='--', alpha=0.5)

    fig.suptitle('Centrality Shift: DiCE Candidates vs Background', fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(out_fig / '04_centrality_boxplot.png')
    plt.close(fig)
    log.info("  Saved centrality boxplot")
    return merged


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 5: DIFFERENTIAL NETWORK REWIRING
# ═══════════════════════════════════════════════════════════════════════════

def analysis_rewiring(dice_df, cent_normal, cent_tumor, out_fig, out_tbl):
    """Volcano plot and histogram of network rewiring."""
    log.info("=" * 60)
    log.info("Analysis 5: Differential Network Rewiring")
    log.info("=" * 60)

    # Merge and compute changes
    norm = cent_normal.copy()
    tum = cent_tumor.copy()
    merged = norm.merge(tum, on='gene', suffixes=('_normal', '_tumor'))

    # Find betweenness columns
    bet_n = [c for c in merged.columns if 'betweenness' in c.lower() and 'normal' in c.lower()]
    bet_t = [c for c in merged.columns if 'betweenness' in c.lower() and 'tumor' in c.lower()]

    if not bet_n or not bet_t:
        log.warning("  Cannot find betweenness columns. Skipping rewiring analysis.")
        return None

    merged['delta_bet'] = merged[bet_t[0]] - merged[bet_n[0]]
    merged['abs_delta'] = merged['delta_bet'].abs()

    # Significance via z-score
    mean_d = merged['delta_bet'].mean()
    std_d = merged['delta_bet'].std()
    merged['z_score'] = (merged['delta_bet'] - mean_d) / std_d if std_d > 0 else 0
    merged['neg_log10_p'] = -np.log10(2 * stats.norm.sf(merged['z_score'].abs()).clip(min=1e-300))

    top_genes = set(dice_df['gene'].head(100))
    merged['is_candidate'] = merged['gene'].isin(top_genes)

    merged.to_csv(out_tbl / '05_rewiring_statistics.csv', index=False)

    # ── Volcano plot ──
    fig, ax = plt.subplots(figsize=(9, 7))
    bg = merged[~merged['is_candidate']]
    fg = merged[merged['is_candidate']]

    ax.scatter(bg['delta_bet'], bg['neg_log10_p'], alpha=0.2, s=8, c='#94A3B8', label='Background')
    ax.scatter(fg['delta_bet'], fg['neg_log10_p'], alpha=0.8, s=30, c='#DC2626',
               edgecolors='black', linewidths=0.3, label='DiCE Top 100', zorder=5)

    # Label top rewired genes
    top_rewired = fg.nlargest(5, 'neg_log10_p')
    for _, row in top_rewired.iterrows():
        ax.annotate(row['gene'], (row['delta_bet'], row['neg_log10_p']),
                    fontsize=8, fontweight='bold', color='#DC2626',
                    xytext=(5, 5), textcoords='offset points')

    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.7, label='p = 0.05')
    ax.set_xlabel('ΔBetweenness Centrality (Tumor − Normal)')
    ax.set_ylabel('$-\\log_{10}$(p-value)')
    ax.set_title('Network Rewiring Volcano Plot')
    ax.legend(frameon=True, loc='upper left')
    plt.tight_layout()
    fig.savefig(out_fig / '05_rewiring_volcano.png')
    plt.close(fig)
    log.info("  Saved volcano plot")

    # ── Histogram ──
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(bg['delta_bet'], bins=60, alpha=0.5, color='#94A3B8', label='Background', density=True)
    ax.hist(fg['delta_bet'], bins=20, alpha=0.7, color='#DC2626', label='DiCE Top 100', density=True)
    ax.set_yscale('log')
    ax.set_xlabel('ΔBetweenness Centrality')
    ax.set_ylabel('Density (log scale)')
    ax.set_title('Distribution of Network Centrality Change')
    ax.legend(frameon=True)
    plt.tight_layout()
    fig.savefig(out_fig / '05_rewiring_histogram.png')
    plt.close(fig)
    log.info("  Saved histogram")
    return merged


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 6: LITERATURE VALIDATION
# ═══════════════════════════════════════════════════════════════════════════

def analysis_literature(dice_df, out_fig, out_tbl, skip_pubmed=False):
    """Query PubMed for disease-related publications."""
    log.info("=" * 60)
    log.info("Analysis 6: Literature Validation")
    log.info("=" * 60)

    top_genes = dice_df['gene'].head(20).tolist()
    records = []

    if skip_pubmed:
        log.info("  PubMed queries skipped (--skip-pubmed flag)")
        # Generate mock data for demonstration
        np.random.seed(42)
        for gene in top_genes:
            count = int(np.random.exponential(50)) + 1
            records.append({'Gene': gene, 'Rank': top_genes.index(gene) + 1,
                            'PubMed_count': count, 'Source': 'mock'})
    else:
        try:
            from Bio import Entrez
            import time
            Entrez.email = 'dice.validation@example.com'

            for i, gene in enumerate(top_genes):
                query = f'"{gene}"[Title/Abstract] AND ("prostate cancer"[Title/Abstract] OR "prostate neoplasm"[MeSH])'
                try:
                    handle = Entrez.esearch(db='pubmed', term=query, retmax=0)
                    result = Entrez.read(handle)
                    handle.close()
                    count = int(result['Count'])
                    records.append({'Gene': gene, 'Rank': i + 1,
                                    'PubMed_count': count, 'Source': 'pubmed'})
                    log.info(f"  {gene}: {count} publications")
                    time.sleep(0.35)  # NCBI rate limit
                except Exception as e:
                    log.warning(f"  PubMed query failed for {gene}: {e}")
                    records.append({'Gene': gene, 'Rank': i + 1,
                                    'PubMed_count': 0, 'Source': 'error'})
        except ImportError:
            log.warning("  Biopython not installed. Using mock data.")
            np.random.seed(42)
            for i, gene in enumerate(top_genes):
                count = int(np.random.exponential(50)) + 1
                records.append({'Gene': gene, 'Rank': i + 1,
                                'PubMed_count': count, 'Source': 'mock'})

    res_df = pd.DataFrame(records)
    res_df.to_csv(out_tbl / '06_literature_counts.csv', index=False)

    # ── Bar chart ──
    fig, ax = plt.subplots(figsize=(10, 6))
    colors_lit = [PALETTE_MAIN[0] if c > 10 else '#94A3B8' for c in res_df['PubMed_count']]
    bars = ax.barh(range(len(res_df)), res_df['PubMed_count'], color=colors_lit, edgecolor='black', linewidth=0.3)
    ax.set_yticks(range(len(res_df)))
    ax.set_yticklabels(res_df['Gene'], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Publications (Prostate Cancer)')
    ax.set_title('Literature Validation — Top 20 DiCE Genes')

    source = res_df['Source'].iloc[0] if len(res_df) > 0 else 'unknown'
    if source == 'mock':
        ax.text(0.98, 0.02, 'Mock data (PubMed not queried)', transform=ax.transAxes,
                fontsize=8, ha='right', va='bottom', style='italic', color='grey')

    plt.tight_layout()
    fig.savefig(out_fig / '06_literature_validation.png')
    plt.close(fig)
    log.info("  Saved literature plot")
    return res_df


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS 7: VIRTUAL KNOCKOUT SIMULATION
# ═══════════════════════════════════════════════════════════════════════════

def analysis_knockout(knockout_df, out_fig, out_tbl):
    """Visualize virtual knockout disruption."""
    log.info("=" * 60)
    log.info("Analysis 7: Virtual Knockout Simulation")
    log.info("=" * 60)

    df = knockout_df.copy()

    # Identify the disruption metric column
    aspl_col = None
    for c in df.columns:
        if 'aspl_change' in c.lower() or 'aspl change' in c.lower():
            aspl_col = c
            break
    if aspl_col is None:
        for c in df.columns:
            if 'change' in c.lower() or 'vitality' in c.lower():
                aspl_col = c
                break
    if aspl_col is None:
        aspl_col = df.columns[-1]  # fallback to last column

    gene_col = df.columns[0]

    df_sorted = df.sort_values(aspl_col, ascending=False).head(30)
    df_sorted.to_csv(out_tbl / '07_knockout_rankings.csv', index=False)

    # ── Bar chart ──
    fig, ax = plt.subplots(figsize=(10, 8))
    values = df_sorted[aspl_col]
    norm_vals = (values - values.min()) / (values.max() - values.min() + 1e-10)
    cmap = plt.cm.RdYlBu_r
    colors_ko = [cmap(v) for v in norm_vals]

    bars = ax.barh(range(len(df_sorted)), df_sorted[aspl_col], color=colors_ko,
                   edgecolor='black', linewidth=0.3)
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted[gene_col], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(aspl_col.replace('_', ' ').title())
    ax.set_title('Virtual Knockout — Network Disruption Score')
    plt.tight_layout()
    fig.savefig(out_fig / '07_knockout_disruption.png')
    plt.close(fig)
    log.info("  Saved knockout plot")
    return df_sorted


# ═══════════════════════════════════════════════════════════════════════════
#  BENCHMARKING: DiCE vs DE vs Degree Centrality
# ═══════════════════════════════════════════════════════════════════════════

def analysis_benchmarking(dice_df, cent_tumor, expr_normal, expr_tumor,
                          known_genes, out_fig, out_tbl):
    """Compare DiCE-Duan against baseline methods."""
    log.info("=" * 60)
    log.info("Benchmarking: DiCE vs DE Ranking vs Degree Centrality")
    log.info("=" * 60)

    all_known = set()
    for genes in known_genes.values():
        all_known |= genes

    methods = {}

    # 1. DiCE ranking (already sorted)
    methods['DiCE-Duan'] = dice_df['gene'].tolist()

    # 2. Differential Expression ranking
    common_genes = list(set(expr_normal.index) & set(expr_tumor.index))
    if common_genes:
        log.info("  Computing differential expression ranking...")
        de_results = []
        for gene in common_genes:
            n_vals = expr_normal.loc[gene].values.astype(float)
            t_vals = expr_tumor.loc[gene].values.astype(float)
            # Filter zeros for log-space data
            n_vals = n_vals[n_vals > 0]
            t_vals = t_vals[t_vals > 0]
            if len(n_vals) >= 3 and len(t_vals) >= 3:
                stat, pval = ttest_ind(t_vals, n_vals, equal_var=False)
                log2fc = np.mean(t_vals) - np.mean(n_vals)  # already log-space
                de_results.append({'gene': gene, 'log2FC': log2fc, 'pvalue': pval,
                                   'abs_log2FC': abs(log2fc)})

        de_df = pd.DataFrame(de_results).sort_values('abs_log2FC', ascending=False)
        methods['DE Ranking'] = de_df['gene'].tolist()
        log.info(f"  DE ranking: {len(de_df)} genes")

    # 3. Degree centrality ranking (from tumor network centrality)
    # Use betweenness as a proxy for degree since we have that
    bet_col = None
    for c in cent_tumor.columns:
        if 'betweenness' in c.lower():
            bet_col = c
            break

    if bet_col:
        degree_ranked = cent_tumor.sort_values(bet_col, ascending=False)['gene'].tolist()
        methods['Degree Centrality'] = degree_ranked
        log.info(f"  Degree ranking: {len(degree_ranked)} genes")

    # Compute Precision@K for all methods
    K_values = [10, 20, 50, 100, 200]
    benchmark_records = []

    for method_name, gene_list in methods.items():
        for k in K_values:
            top_k = set(gene_list[:k])
            prec = len(top_k & all_known) / k
            benchmark_records.append({'Method': method_name, 'K': k, 'Precision': round(prec, 4)})

    bench_df = pd.DataFrame(benchmark_records)
    bench_df.to_csv(out_tbl / '08_benchmark_comparison.csv', index=False)

    # ── Precision@K Comparison Plot ──
    fig, ax = plt.subplots(figsize=(9, 6))
    method_colors = {'DiCE-Duan': PALETTE_MAIN[0], 'DE Ranking': PALETTE_MAIN[1],
                     'Degree Centrality': PALETTE_MAIN[3]}
    method_markers = {'DiCE-Duan': 'o', 'DE Ranking': 's', 'Degree Centrality': 'D'}

    for method_name in methods:
        subset = bench_df[bench_df['Method'] == method_name]
        ax.plot(subset['K'], subset['Precision'],
                marker=method_markers.get(method_name, '^'),
                color=method_colors.get(method_name, '#666'),
                label=method_name, linewidth=2.5, markersize=8)

    ax.set_xlabel('K (Number of Top Genes)')
    ax.set_ylabel('Precision@K (Known Disease Gene Overlap)')
    ax.set_title('Benchmarking: DiCE-Duan vs Baseline Methods')
    ax.legend(frameon=True, loc='best')
    ax.set_xticks(K_values)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_fig / '08_benchmark_comparison.png')
    plt.close(fig)
    log.info("  Saved benchmark plot")

    # ── Also create a grouped bar chart for a specific K ──
    fig, ax = plt.subplots(figsize=(8, 5))
    k_target = 100
    subset_k = bench_df[bench_df['K'] == k_target]
    if len(subset_k) > 0:
        x = range(len(subset_k))
        colors_bench = [method_colors.get(m, '#666') for m in subset_k['Method']]
        bars = ax.bar(x, subset_k['Precision'], color=colors_bench, width=0.5,
                      edgecolor='black', linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(subset_k['Method'], fontsize=11)
        ax.set_ylabel(f'Precision@{k_target}')
        ax.set_title(f'Method Comparison at K={k_target}')

        for bar, val in zip(bars, subset_k['Precision']):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                    f'{val:.3f}', ha='center', va='bottom', fontweight='bold')

        ax.set_ylim(0, max(subset_k['Precision']) * 1.25 + 0.01)
        plt.tight_layout()
        fig.savefig(out_fig / '08_benchmark_bar_k100.png')
        plt.close(fig)
        log.info("  Saved benchmark bar chart")

    return bench_df


# ═══════════════════════════════════════════════════════════════════════════
#  SUMMARY REPORT
# ═══════════════════════════════════════════════════════════════════════════

def write_summary(out_dir, results):
    """Write a text summary of all validation results."""
    summary_path = out_dir / 'validation_summary.txt'
    lines = [
        "=" * 70,
        "DiCE-Duan Validation Pipeline — Summary Report",
        "=" * 70,
        ""
    ]

    if 'hypergeometric' in results and results['hypergeometric'] is not None:
        lines.append("HYPERGEOMETRIC SIGNIFICANCE TEST")
        lines.append("-" * 40)
        for _, row in results['hypergeometric'].iterrows():
            lines.append(f"  {row['Database']}: overlap={row['Overlap']}, "
                         f"expected={row['Expected_overlap']}, "
                         f"fold={row['Fold_enrichment']}, p={row['P_value']:.2e}")
        lines.append("")

    if 'precision' in results and results['precision'] is not None:
        lines.append("PRECISION@K")
        lines.append("-" * 40)
        for _, row in results['precision'].iterrows():
            lines.append(f"  K={int(row['K'])}: combined={row['Combined_precision']:.4f}")
        lines.append("")

    if 'benchmark' in results and results['benchmark'] is not None:
        lines.append("BENCHMARKING")
        lines.append("-" * 40)
        for _, row in results['benchmark'].iterrows():
            lines.append(f"  {row['Method']} @ K={int(row['K'])}: precision={row['Precision']:.4f}")
        lines.append("")

    lines.extend([
        "=" * 70,
        "Output files generated in:",
        f"  Figures: {out_dir / 'figures'}",
        f"  Tables:  {out_dir / 'tables'}",
        "=" * 70
    ])

    summary_path.write_text('\n'.join(lines))
    log.info(f"Summary written to {summary_path}")


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='DiCE-Duan Validation Pipeline')
    parser.add_argument('--all', action='store_true', help='Run all analyses')
    parser.add_argument('--analysis', nargs='+', type=int, choices=[1,2,3,4,5,6,7,8,9,10],
                        help='Run specific analyses (1-8 + 9=fitness + 10=dual-dataset)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: data/results/validation/)')
    parser.add_argument('--skip-pubmed', action='store_true',
                        help='Skip PubMed queries (uses mock data for analysis 6)')
    args = parser.parse_args()

    if not args.all and not args.analysis:
        args.all = True  # default to all

    analyses = set(range(1, 11)) if args.all else set(args.analysis)

    # Setup paths
    root = find_project_root()
    out_dir = Path(args.output_dir) if args.output_dir else root / 'data' / 'results' / 'validation'
    out_fig = out_dir / 'figures'
    out_tbl = out_dir / 'tables'
    out_fig.mkdir(parents=True, exist_ok=True)
    out_tbl.mkdir(parents=True, exist_ok=True)

    log.info("DiCE-Duan Validation Pipeline")
    log.info(f"Project root: {root}")
    log.info(f"Output dir:   {out_dir}")
    log.info(f"Analyses:     {sorted(analyses)}")
    log.info("")

    # Load core data
    dice_df = load_dice_genes(root)
    universe_size = len(dice_df) if len(dice_df) > 500 else 8962  # approximate

    # Load full gene list for universe size
    full_genes_path = root / 'data' / 'results' / 'dice_genes.txt'
    if full_genes_path.exists():
        full_df = pd.read_csv(full_genes_path, sep='\t')
        universe_size = len(full_df)
        log.info(f"Gene universe size: {universe_size}")

    results = {}

    # Analysis 1: Pathway Enrichment
    if 1 in analyses:
        results['enrichment'] = analysis_pathway_enrichment(dice_df, root, out_fig, out_tbl)

    # Analysis 2: Hypergeometric Test
    if 2 in analyses:
        known = load_known_genes(root)
        if known:
            results['hypergeometric'] = analysis_hypergeometric(dice_df, known, universe_size, out_fig, out_tbl)

    # Analysis 3: Precision@K
    if 3 in analyses:
        known = load_known_genes(root)
        if known:
            results['precision'] = analysis_precision_at_k(dice_df, known, out_fig, out_tbl)

    # Analysis 4: Network Topology
    if 4 in analyses:
        cent_n = load_centrality(root, 'normal')
        cent_t = load_centrality(root, 'tumor')
        results['topology'] = analysis_network_topology(dice_df, cent_n, cent_t, out_fig, out_tbl)

    # Analysis 5: Differential Rewiring
    if 5 in analyses:
        cent_n = load_centrality(root, 'normal')
        cent_t = load_centrality(root, 'tumor')
        results['rewiring'] = analysis_rewiring(dice_df, cent_n, cent_t, out_fig, out_tbl)

    # Analysis 6: Literature Validation
    if 6 in analyses:
        results['literature'] = analysis_literature(dice_df, out_fig, out_tbl,
                                                    skip_pubmed=args.skip_pubmed)

    # Analysis 7: Virtual Knockout
    if 7 in analyses:
        knockout_df = load_knockout(root)
        results['knockout'] = analysis_knockout(knockout_df, out_fig, out_tbl)

    # Analysis 8: Benchmarking
    if 8 in analyses:
        known = load_known_genes(root)
        cent_t = load_centrality(root, 'tumor')
        expr_n, expr_t = load_expression(root)
        if known:
            results['benchmark'] = analysis_benchmarking(dice_df, cent_t, expr_n, expr_t,
                                                          known, out_fig, out_tbl)

    # Analysis 9: Cancer Fitness Gene Overlap
    if 9 in analyses:
        try:
            from fitness_analysis import run_fitness_gene_analysis
            results['fitness'] = run_fitness_gene_analysis(
                dice_df, root, out_fig, out_tbl,
                universe_size=universe_size
            )
        except ImportError:
            log.warning("Could not import fitness_analysis module. Skipping Analysis 9.")
        except Exception as e:
            log.error(f"Analysis 9 failed: {e}")

    # Analysis 10: Dual-Dataset Cross-Validation
    if 10 in analyses:
        try:
            from dual_dataset_analysis import run_dual_dataset_analysis
            # Look for second dataset results
            alt_paths = [
                root / 'data' / 'results_gse21032' / 'dice_genes_top500.txt',
                root / 'data' / 'results_alt' / 'dice_genes_top500.txt',
                root / 'data' / 'results_metastasis' / 'dice_genes_top500.txt',
            ]
            alt_df = None
            alt_label = 'GSE21032'
            for p in alt_paths:
                if p.exists():
                    alt_df = pd.read_csv(p, sep='\t')
                    alt_label = p.parent.name.replace('results_', '').upper()
                    log.info(f"Found second dataset: {p}")
                    break

            if alt_df is not None:
                results['dual_dataset'] = run_dual_dataset_analysis(
                    dice_df, alt_df, 'TCGA', alt_label,
                    root, out_fig, out_tbl,
                    universe_size=universe_size
                )
            else:
                log.info("Analysis 10: No second dataset found. Run the pipeline on a ")
                log.info("  second dataset and place results in data/results_gse21032/")
                log.info("  Then re-run: python validation_pipeline.py --analysis 10")
                results['dual_dataset'] = {'status': 'SKIPPED', 'reason': 'No second dataset'}
        except ImportError:
            log.warning("Could not import dual_dataset_analysis module. Skipping Analysis 10.")
        except Exception as e:
            log.error(f"Analysis 10 failed: {e}")

    # Summary
    write_summary(out_dir, results)

    log.info("")
    log.info("=" * 60)
    log.info("PIPELINE COMPLETE")
    log.info(f"Figures: {out_fig}")
    log.info(f"Tables:  {out_tbl}")
    log.info("=" * 60)


if __name__ == '__main__':
    main()
