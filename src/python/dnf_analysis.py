#!/usr/bin/env python3
"""
DNF Analysis & Visualization
==============================
Post-processes the C++ DNF pairwise knockout output to produce:
1. Ranked gene pair table (by DSS)
2. DSS heatmap of top genes
3. Scatter plot: ΔASPL_tumor vs ΔASPL_normal (colored by DSS)
4. SynLethDB cross-reference (if data available)
5. Summary statistics and publication-ready CSV

Usage:
    python dnf_analysis.py \
        --dnf-results data/results/dnf_pairwise_results.txt \
        --dice-genes data/results/dice_genes_top500.txt \
        --output-dir data/results/dnf
"""

import argparse
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
log = logging.getLogger(__name__)

# ─── Style ───────────────────────────────────────────────────────────────────
STYLE = {
    'figure.dpi': 300, 'savefig.dpi': 300, 'font.size': 11,
    'axes.titlesize': 13, 'axes.labelsize': 12,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 10, 'figure.figsize': (10, 7),
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.15,
}
plt.rcParams.update(STYLE)


# ═══════════════════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════

def load_dnf_results(filepath: Path) -> pd.DataFrame:
    """Load the C++ DNF pairwise output."""
    df = pd.read_csv(filepath, sep='\t')
    log.info(f"Loaded {len(df)} gene pairs from {filepath.name}")
    return df


def load_dice_genes(filepath: Path) -> pd.DataFrame:
    """Load DiCE gene rankings for context."""
    if filepath.exists():
        df = pd.read_csv(filepath, sep='\t')
        log.info(f"Loaded {len(df)} DiCE genes from {filepath.name}")
        return df
    log.warning(f"DiCE gene file not found: {filepath}")
    return pd.DataFrame()


def load_synlethdb(filepath: Optional[Path]) -> Optional[pd.DataFrame]:
    """Load SynLethDB data if available."""
    if filepath and filepath.exists():
        try:
            df = pd.read_csv(filepath, sep='\t')
            log.info(f"Loaded {len(df)} SL pairs from SynLethDB")
            return df
        except Exception as e:
            log.warning(f"Failed to load SynLethDB: {e}")
    return None


# ═══════════════════════════════════════════════════════════════════════════
#  ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

def compute_summary_stats(df: pd.DataFrame) -> dict:
    """Compute summary statistics for the DNF results."""
    summary = {
        'total_pairs': len(df),
        'mean_dss': df['dss'].mean(),
        'std_dss': df['dss'].std(),
        'max_dss': df['dss'].max(),
        'min_dss': df['dss'].min(),
        'pairs_positive_dss': (df['dss'] > 0).sum(),
        'pairs_negative_dss': (df['dss'] < 0).sum(),
        'mean_synergy_tumor': df['synergy_tumor'].mean(),
        'mean_synergy_normal': df['synergy_normal'].mean(),
        'pairs_synergistic_in_tumor': (df['synergy_tumor'] > 0).sum(),
        'pairs_diff_synergy_positive': (df['differential_synergy'] > 0).sum(),
    }

    # Statistical test: Is the mean DSS significantly different from 0?
    if len(df) > 2:
        t_stat, p_val = stats.ttest_1samp(df['dss'].dropna(), 0)
        summary['dss_ttest_statistic'] = t_stat
        summary['dss_ttest_pvalue'] = p_val

    return summary


def cross_reference_synlethdb(
    dnf_df: pd.DataFrame,
    sl_df: pd.DataFrame,
    top_k: int = 50
) -> dict:
    """
    Cross-reference top DNF pairs against SynLethDB known synthetic lethal pairs.
    Returns enrichment statistics.
    """
    # Identify gene pair columns in SynLethDB
    gene_col_a, gene_col_b = None, None
    for ca, cb in [('gene_a', 'gene_b'), ('GeneA', 'GeneB'),
                   ('gene1', 'gene2'), ('Gene1Symbol', 'Gene2Symbol')]:
        if ca in sl_df.columns and cb in sl_df.columns:
            gene_col_a, gene_col_b = ca, cb
            break

    if gene_col_a is None:
        log.warning(f"Cannot identify gene columns in SynLethDB. Columns: {sl_df.columns.tolist()}")
        return {'status': 'failed', 'reason': 'column_mismatch'}

    # Build set of known SL pairs (unordered)
    known_sl = set()
    for _, row in sl_df.iterrows():
        pair = frozenset([str(row[gene_col_a]).upper(), str(row[gene_col_b]).upper()])
        known_sl.add(pair)

    log.info(f"Known SL pairs in database: {len(known_sl)}")

    # Check top-K DNF pairs
    top_dnf = dnf_df.head(top_k)
    predicted_pairs = []
    hits = 0

    for _, row in top_dnf.iterrows():
        pair = frozenset([str(row['gene_a']).upper(), str(row['gene_b']).upper()])
        is_sl = pair in known_sl
        if is_sl:
            hits += 1
        predicted_pairs.append({
            'gene_a': row['gene_a'],
            'gene_b': row['gene_b'],
            'dss': row['dss'],
            'known_sl': is_sl,
        })

    # Hypergeometric test for enrichment
    # M = total possible pairs in universe (~20000 genes → ~200M pairs)
    # n = known SL pairs
    # N = our top-K predictions
    # k = hits
    from scipy.stats import hypergeom
    M = 20000 * 19999 // 2  # Gene universe pairs
    n = len(known_sl)
    N = top_k
    k = hits

    p_value = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
    expected = n * N / M if M > 0 else 0
    fold_enrichment = k / expected if expected > 0 else 0

    result = {
        'status': 'complete',
        'top_k': top_k,
        'hits': hits,
        'hit_rate': hits / top_k * 100 if top_k > 0 else 0,
        'known_sl_total': len(known_sl),
        'expected_hits': expected,
        'fold_enrichment': fold_enrichment,
        'pvalue': p_value,
        'hit_pairs': [p for p in predicted_pairs if p['known_sl']],
    }

    log.info(f"SynLethDB validation: {hits}/{top_k} hits "
             f"({result['hit_rate']:.1f}%), fold={fold_enrichment:.1f}x, P={p_value:.2e}")

    return result


# ═══════════════════════════════════════════════════════════════════════════
#  VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════

def plot_dss_heatmap(df: pd.DataFrame, output_path: Path, top_n: int = 30):
    """
    Generate a symmetric heatmap of DSS values for top genes.
    """
    # Get all unique genes from top pairs
    top_pairs = df.head(top_n * 3)  # Get enough pairs to fill the matrix
    genes = list(set(top_pairs['gene_a'].tolist() + top_pairs['gene_b'].tolist()))[:top_n]

    if len(genes) < 3:
        log.warning("Not enough genes for heatmap")
        return

    # Build matrix
    gene_idx = {g: i for i, g in enumerate(genes)}
    matrix = np.zeros((len(genes), len(genes)))
    matrix[:] = np.nan

    for _, row in df.iterrows():
        ga, gb = row['gene_a'], row['gene_b']
        if ga in gene_idx and gb in gene_idx:
            i, j = gene_idx[ga], gene_idx[gb]
            matrix[i][j] = row['dss']
            matrix[j][i] = row['dss']

    # Plot
    fig, ax = plt.subplots(figsize=(14, 12))
    mask = np.isnan(matrix)

    sns.heatmap(matrix, mask=mask, ax=ax,
                xticklabels=genes, yticklabels=genes,
                cmap='RdBu_r', center=0,
                linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Differential Synergy Score (DSS)',
                          'shrink': 0.8})

    ax.set_title('Differential Network Fragility: Pairwise DSS Heatmap\n'
                 'Positive (Red) = Tumor-selective vulnerability',
                 fontsize=14, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log.info(f"Saved DSS heatmap → {output_path}")


def plot_differential_scatter(df: pd.DataFrame, output_path: Path, top_n: int = 200):
    """
    Scatter plot: ΔASPL_tumor vs ΔASPL_normal, colored by DSS.
    Points above the diagonal = tumor-selective vulnerability.
    """
    plot_df = df.head(top_n).copy()

    fig, ax = plt.subplots(figsize=(10, 10))

    scatter = ax.scatter(
        plot_df['delta_aspl_normal'],
        plot_df['delta_aspl_tumor'],
        c=plot_df['dss'],
        cmap='RdBu_r',
        s=40,
        alpha=0.7,
        edgecolors='white',
        linewidth=0.5
    )

    plt.colorbar(scatter, ax=ax, label='DSS (Differential Synergy Score)', shrink=0.8)

    # Diagonal line (equal disruption)
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, 'k--', alpha=0.3, linewidth=1, label='Equal disruption')
    ax.fill_between(lims, lims, [lims[1], lims[1]], alpha=0.05, color='red',
                     label='Tumor-selective zone')

    # Label top 5 pairs
    top5 = plot_df.head(5)
    for _, row in top5.iterrows():
        label = f"{row['gene_a']}-{row['gene_b']}"
        ax.annotate(label,
                    (row['delta_aspl_normal'], row['delta_aspl_tumor']),
                    fontsize=7, fontweight='bold', color='darkred',
                    textcoords="offset points", xytext=(5, 5),
                    arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

    ax.set_xlabel('ΔASPL Normal (Network Disruption in Healthy Tissue)', fontsize=12)
    ax.set_ylabel('ΔASPL Tumor (Network Disruption in Cancer)', fontsize=12)
    ax.set_title('Differential Network Fragility: Tumor vs Normal Disruption\n'
                 'Points above diagonal = tumor-selective targets',
                 fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log.info(f"Saved differential scatter → {output_path}")


def plot_synergy_distribution(df: pd.DataFrame, output_path: Path):
    """
    Distribution of DSS scores with tumor-selective zone highlighted.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # DSS distribution
    ax = axes[0]
    ax.hist(df['dss'], bins=50, color='#7C3AED', alpha=0.7, edgecolor='white')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=df['dss'].mean(), color='red', linestyle='-', linewidth=1.5,
               label=f"Mean DSS = {df['dss'].mean():.4f}")

    pos_pct = (df['dss'] > 0).sum() / len(df) * 100
    ax.set_title(f'Distribution of DSS\n{pos_pct:.0f}% pairs have positive DSS (tumor-selective)',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('Differential Synergy Score (DSS)')
    ax.set_ylabel('Number of Pairs')
    ax.legend()

    # Differential synergy distribution
    ax = axes[1]
    ax.hist(df['differential_synergy'], bins=50, color='#DC2626', alpha=0.7, edgecolor='white')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)

    syn_pos = (df['differential_synergy'] > 0).sum() / len(df) * 100
    ax.set_title(f'Distribution of Differential Synergy\n{syn_pos:.0f}% pairs are synergistic in tumor',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('Differential Synergy (Tumor − Normal)')
    ax.set_ylabel('Number of Pairs')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log.info(f"Saved distribution plot → {output_path}")


def plot_top_pairs_bar(df: pd.DataFrame, output_path: Path, top_n: int = 20):
    """Horizontal bar chart of top DSS pairs."""
    top = df.head(top_n).copy()
    top['pair_label'] = top['gene_a'] + ' — ' + top['gene_b']
    top = top.iloc[::-1]  # Reverse for horizontal bar (top at top)

    fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.35)))

    colors = ['#DC2626' if d > 0 else '#2563EB' for d in top['dss']]
    ax.barh(range(len(top)), top['dss'], color=colors, edgecolor='white', linewidth=0.5, height=0.7)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['pair_label'], fontsize=9, fontfamily='monospace')
    ax.set_xlabel('Differential Synergy Score (DSS)', fontsize=12)
    ax.set_title(f'Top {top_n} Tumor-Selective Vulnerable Gene Pairs',
                 fontsize=13, fontweight='bold')
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log.info(f"Saved top pairs bar → {output_path}")


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='DNF Pairwise Knockout Analysis & Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--dnf-results', required=True,
                        help='C++ DNF pairwise results (TSV)')
    parser.add_argument('--dice-genes', default=None,
                        help='DiCE gene rankings for context')
    parser.add_argument('--synlethdb', default=None,
                        help='SynLethDB file for validation')
    parser.add_argument('--output-dir', default='data/results/dnf',
                        help='Output directory')
    parser.add_argument('--top-k', type=int, default=50,
                        help='Top K pairs for SynLethDB validation')

    args = parser.parse_args()
    out_dir = Path(args.output_dir)
    fig_dir = out_dir / 'figures'
    tbl_dir = out_dir / 'tables'
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load data
    log.info("=" * 60)
    log.info("DIFFERENTIAL NETWORK FRAGILITY — Post-Analysis")
    log.info("=" * 60)

    dnf_df = load_dnf_results(Path(args.dnf_results))

    dice_df = pd.DataFrame()
    if args.dice_genes:
        dice_df = load_dice_genes(Path(args.dice_genes))

    # 2. Compute summary statistics
    summary = compute_summary_stats(dnf_df)
    log.info(f"\n  Summary Statistics:")
    log.info(f"    Total pairs:          {summary['total_pairs']}")
    log.info(f"    Mean DSS:             {summary['mean_dss']:.6f}")
    log.info(f"    Max DSS:              {summary['max_dss']:.6f}")
    log.info(f"    Positive DSS pairs:   {summary['pairs_positive_dss']}")
    log.info(f"    Diff-synergistic:     {summary['pairs_diff_synergy_positive']}")

    if 'dss_ttest_pvalue' in summary:
        log.info(f"    DSS ≠ 0 (t-test):     P={summary['dss_ttest_pvalue']:.2e}")

    # Save summary
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(tbl_dir / 'dnf_summary.csv', index=False)

    # 3. Save ranked pairs
    dnf_df['pair_rank'] = range(1, len(dnf_df) + 1)
    dnf_df.to_csv(tbl_dir / 'dnf_ranked_pairs.csv', index=False, sep='\t')
    log.info(f"Saved ranked pairs → {tbl_dir / 'dnf_ranked_pairs.csv'}")

    # Save top 50
    dnf_df.head(50).to_csv(tbl_dir / 'dnf_top50_pairs.csv', index=False, sep='\t')

    # 4. Generate visualizations
    log.info("\nGenerating visualizations...")

    plot_dss_heatmap(dnf_df, fig_dir / 'dnf_heatmap.png')
    plot_differential_scatter(dnf_df, fig_dir / 'dnf_scatter.png')
    plot_synergy_distribution(dnf_df, fig_dir / 'dnf_distribution.png')
    plot_top_pairs_bar(dnf_df, fig_dir / 'dnf_top_pairs.png')

    # 5. SynLethDB cross-reference
    sl_result = None
    if args.synlethdb:
        sl_df = load_synlethdb(Path(args.synlethdb))
        if sl_df is not None:
            sl_result = cross_reference_synlethdb(dnf_df, sl_df, top_k=args.top_k)

            if sl_result['status'] == 'complete':
                # Save SL validation
                sl_summary = {k: v for k, v in sl_result.items() if k != 'hit_pairs'}
                pd.DataFrame([sl_summary]).to_csv(
                    tbl_dir / 'dnf_synlethdb_validation.csv', index=False)

                if sl_result['hit_pairs']:
                    pd.DataFrame(sl_result['hit_pairs']).to_csv(
                        tbl_dir / 'dnf_synlethdb_hits.csv', index=False)

    # 6. Print top results
    log.info("\n" + "=" * 70)
    log.info("TOP 10 TUMOR-SELECTIVE VULNERABLE PAIRS")
    log.info("=" * 70)
    top10 = dnf_df.head(10)
    for _, row in top10.iterrows():
        log.info(f"  {row['gene_a']:>10} — {row['gene_b']:<10}  "
                 f"DSS={row['dss']:+.4f}  "
                 f"Syn_T={row['synergy_tumor']:+.4f}  "
                 f"Syn_N={row['synergy_normal']:+.4f}")
    log.info("=" * 70)

    log.info(f"\nAll outputs saved to: {out_dir}")
    log.info("DNF post-analysis complete.")


if __name__ == '__main__':
    main()
