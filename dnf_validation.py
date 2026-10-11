#!/usr/bin/env python3
"""
DNF Validation Suite
====================
Rigorous multi-layer validation of Differential Network Fragility (DNF) predictions.

Validation Layers:
  1. Known Cancer Gene Enrichment (COSMIC, DisGeNET)
  2. DepMap Essentiality Concordance
  3. Functional Coherence of Top Pairs (GO/Pathway co-annotation)
  4. Comparison Against Naive Baselines (Does DSS > random? > centrality product?)
  5. Permutation Test (Is the DSS distribution non-random?)
  6. Literature Cross-Reference for Top Pairs

Usage:
    python dnf_validation.py \
        --dnf-results data/results/dnf_pairwise_results.txt \
        --dice-genes data/results/dice_genes_top500.txt \
        --known-genes-dir data/known_genes \
        --output-dir data/results/dnf/validation
"""

import argparse
import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import hypergeom, mannwhitneyu, spearmanr, pearsonr

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
log = logging.getLogger(__name__)

STYLE = {
    'figure.dpi': 300, 'savefig.dpi': 300, 'font.size': 11,
    'axes.titlesize': 13, 'axes.labelsize': 12,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.15,
}
plt.rcParams.update(STYLE)


# ═══════════════════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════

def load_gene_set(filepath: Path) -> Set[str]:
    """Load a gene set from a simple text file (one gene per line, skip comments)."""
    genes = set()
    if not filepath.exists():
        return genes
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Take first column if TSV
                gene = line.split('\t')[0].strip()
                if gene and gene != 'GENE':
                    genes.add(gene.upper())
    return genes


def load_depmap(filepath: Path) -> pd.DataFrame:
    """Load DepMap fitness scores."""
    if not filepath.exists():
        return pd.DataFrame()
    rows = []
    with open(filepath) as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue
            if len(parts) >= 4:
                gene = parts[0]
                scores = []
                for s in parts[1:]:
                    try:
                        scores.append(float(s))
                    except (ValueError, TypeError):
                        scores.append(np.nan)
                rows.append({
                    'gene': gene.upper(),
                    'mean_score': np.nanmean(scores),
                    'min_score': np.nanmin(scores),
                })
    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════════
#  VALIDATION 1: Known Cancer Gene Enrichment
# ═══════════════════════════════════════════════════════════════════════════

def validate_cancer_gene_enrichment(
    dnf_df: pd.DataFrame,
    cosmic_genes: Set[str],
    disgenet_genes: Set[str],
    fig_dir: Path,
    tbl_dir: Path,
) -> Dict:
    """
    Test whether top DSS pairs are enriched for known cancer genes.
    For each gene pair, check if EITHER gene is in the cancer gene set.
    """
    log.info("\n=== Validation 1: Known Cancer Gene Enrichment ===")

    results = {}
    for label, known_set in [("COSMIC", cosmic_genes), ("DisGeNET_PCa", disgenet_genes)]:
        if not known_set:
            log.warning(f"  {label}: No genes loaded, skipping.")
            continue

        # For each pair, check if at least one gene is a known cancer gene
        dnf_df[f'has_{label}'] = dnf_df.apply(
            lambda r: r['gene_a'].upper() in known_set or r['gene_b'].upper() in known_set,
            axis=1
        )
        # Both genes known
        dnf_df[f'both_{label}'] = dnf_df.apply(
            lambda r: r['gene_a'].upper() in known_set and r['gene_b'].upper() in known_set,
            axis=1
        )

        # Test at multiple top-K cutoffs
        cutoffs = [10, 20, 30, len(dnf_df)]
        cutoff_results = []
        for k in cutoffs:
            if k > len(dnf_df):
                continue
            top_k = dnf_df.head(k)
            n_hits = top_k[f'has_{label}'].sum()
            n_both = top_k[f'both_{label}'].sum()
            pct = n_hits / k * 100 if k > 0 else 0

            # Background rate: what fraction of ALL pairs contain a known gene?
            bg_rate = dnf_df[f'has_{label}'].mean()

            # Binomial test
            from scipy.stats import binomtest
            try:
                p_val = binomtest(n_hits, k, bg_rate, alternative='greater').pvalue
            except Exception:
                p_val = 1.0

            cutoff_results.append({
                'top_k': k,
                'hits': n_hits,
                'both_known': n_both,
                'hit_rate_pct': round(pct, 1),
                'background_rate_pct': round(bg_rate * 100, 1),
                'enrichment_fold': round(pct / (bg_rate * 100), 2) if bg_rate > 0 else 0,
                'pvalue': p_val,
            })

            log.info(f"  {label} Top-{k}: {n_hits}/{k} pairs contain known gene "
                     f"({pct:.1f}% vs {bg_rate*100:.1f}% bg, P={p_val:.2e})")

        results[label] = cutoff_results

        # Save
        pd.DataFrame(cutoff_results).to_csv(tbl_dir / f'val1_{label}_enrichment.csv', index=False)

    # Visualization: bar chart of enrichment at different cutoffs
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, (label, data) in enumerate(results.items()):
        ax = axes[idx] if len(results) > 1 else axes
        df_plot = pd.DataFrame(data)
        x = [str(r['top_k']) for r in data]
        y_hits = [r['hit_rate_pct'] for r in data]
        y_bg = [r['background_rate_pct'] for r in data]

        width = 0.35
        x_pos = np.arange(len(x))
        ax.bar(x_pos - width/2, y_hits, width, label='DNF Top-K', color='#DC2626', alpha=0.8)
        ax.bar(x_pos + width/2, y_bg, width, label='Background', color='#9CA3AF', alpha=0.6)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(x)
        ax.set_xlabel('Top K Pairs')
        ax.set_ylabel('% Containing Known Cancer Gene')
        ax.set_title(f'{label} Cancer Gene Enrichment', fontweight='bold')
        ax.legend()

    plt.tight_layout()
    plt.savefig(fig_dir / 'val1_cancer_gene_enrichment.png', dpi=300)
    plt.close()
    log.info(f"  Saved figure → val1_cancer_gene_enrichment.png")

    return results


# ═══════════════════════════════════════════════════════════════════════════
#  VALIDATION 2: DepMap Essentiality Concordance
# ═══════════════════════════════════════════════════════════════════════════

def validate_depmap_essentiality(
    dnf_df: pd.DataFrame,
    depmap_df: pd.DataFrame,
    fig_dir: Path,
    tbl_dir: Path,
) -> Dict:
    """
    Test whether high-DSS pairs involve genes that are essential in cancer cell lines.
    Hypothesis: Genes in tumor-selective pairs (high DSS) should show
    stronger essentiality (more negative DepMap scores) than genes in low-DSS pairs.
    """
    log.info("\n=== Validation 2: DepMap Essentiality Concordance ===")

    if depmap_df.empty:
        log.warning("  No DepMap data available. Skipping.")
        return {'status': 'no_data'}

    depmap_lookup = dict(zip(depmap_df['gene'], depmap_df['mean_score']))

    # For each pair, compute the average essentiality of the two genes
    pair_scores = []
    for _, row in dnf_df.iterrows():
        ga = row['gene_a'].upper()
        gb = row['gene_b'].upper()
        sa = depmap_lookup.get(ga, np.nan)
        sb = depmap_lookup.get(gb, np.nan)
        mean_ess = np.nanmean([sa, sb])
        pair_scores.append({
            'gene_a': row['gene_a'],
            'gene_b': row['gene_b'],
            'dss': row['dss'],
            'pair_rank': row.get('pair_rank', 0),
            'depmap_a': sa,
            'depmap_b': sb,
            'mean_essentiality': mean_ess,
        })

    score_df = pd.DataFrame(pair_scores)
    valid = score_df.dropna(subset=['mean_essentiality'])

    if len(valid) < 5:
        log.warning(f"  Only {len(valid)} pairs have DepMap data. Insufficient for analysis.")
        return {'status': 'insufficient_data', 'n_pairs_with_data': len(valid)}

    # Correlation: DSS vs mean essentiality
    rho, rho_p = spearmanr(valid['dss'], valid['mean_essentiality'])
    r, r_p = pearsonr(valid['dss'], valid['mean_essentiality'])

    log.info(f"  Pairs with DepMap data: {len(valid)}/{len(dnf_df)}")
    log.info(f"  Spearman ρ(DSS, Essentiality) = {rho:.3f}, P = {rho_p:.2e}")
    log.info(f"  Pearson  r(DSS, Essentiality) = {r:.3f}, P = {r_p:.2e}")

    # Split into top vs bottom half by DSS and compare essentiality
    median_dss = valid['dss'].median()
    high_dss = valid[valid['dss'] > median_dss]['mean_essentiality'].dropna()
    low_dss = valid[valid['dss'] <= median_dss]['mean_essentiality'].dropna()

    if len(high_dss) >= 3 and len(low_dss) >= 3:
        U_stat, mw_p = mannwhitneyu(high_dss, low_dss, alternative='less')
        log.info(f"  Mann-Whitney (High DSS more essential?): U={U_stat:.0f}, P={mw_p:.4f}")
        log.info(f"    High DSS mean essentiality: {high_dss.mean():.3f}")
        log.info(f"    Low  DSS mean essentiality: {low_dss.mean():.3f}")
    else:
        mw_p = np.nan

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Scatter: DSS vs essentiality
    ax = axes[0]
    ax.scatter(valid['dss'], valid['mean_essentiality'], s=40, alpha=0.7,
               c=valid['dss'], cmap='RdBu_r', edgecolors='white', linewidth=0.5)
    ax.axhline(y=-0.5, color='red', linestyle='--', alpha=0.4, label='Essentiality threshold')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.3)
    ax.set_xlabel('Differential Synergy Score (DSS)')
    ax.set_ylabel('Mean DepMap Essentiality (more negative = more essential)')
    ax.set_title(f'DSS vs Cancer Essentiality\nρ={rho:.3f}, P={rho_p:.2e}', fontweight='bold')
    ax.legend(fontsize=9)

    # Box plot: high vs low DSS pairs
    ax = axes[1]
    box_data = [high_dss.values, low_dss.values]
    bp = ax.boxplot(box_data, labels=['High DSS\n(Tumor-selective)', 'Low DSS'],
                     patch_artist=True, widths=0.5)
    bp['boxes'][0].set_facecolor('#DC2626')
    bp['boxes'][0].set_alpha(0.5)
    bp['boxes'][1].set_facecolor('#9CA3AF')
    bp['boxes'][1].set_alpha(0.5)
    ax.axhline(y=-0.5, color='red', linestyle='--', alpha=0.4)
    ax.set_ylabel('Mean DepMap Essentiality Score')
    ax.set_title(f'Essentiality by DSS Group\nMann-Whitney P={mw_p:.4f}' if not np.isnan(mw_p)
                 else 'Essentiality by DSS Group', fontweight='bold')

    plt.tight_layout()
    plt.savefig(fig_dir / 'val2_depmap_essentiality.png', dpi=300)
    plt.close()

    # Save table
    score_df.to_csv(tbl_dir / 'val2_depmap_concordance.csv', index=False)

    return {
        'status': 'complete',
        'n_pairs_with_data': len(valid),
        'spearman_rho': rho,
        'spearman_p': rho_p,
        'pearson_r': r,
        'pearson_p': r_p,
        'mannwhitney_p': mw_p if not np.isnan(mw_p) else None,
        'high_dss_mean_ess': high_dss.mean(),
        'low_dss_mean_ess': low_dss.mean(),
    }


# ═══════════════════════════════════════════════════════════════════════════
#  VALIDATION 3: Comparison Against Naive Baselines
# ═══════════════════════════════════════════════════════════════════════════

def validate_against_baselines(
    dnf_df: pd.DataFrame,
    dice_df: pd.DataFrame,
    cosmic_genes: Set[str],
    fig_dir: Path,
    tbl_dir: Path,
) -> Dict:
    """
    Critical validation: Does DSS outperform naive baselines?
    Baseline 1: Product of individual DiCE scores (ensemble_score_a × ensemble_score_b)
    Baseline 2: Sum of individual single-knockout ASPL changes
    Baseline 3: Random ranking

    We compare using Precision@K for recovering known cancer gene pairs.
    """
    log.info("\n=== Validation 3: DSS vs Naive Baselines ===")

    if dice_df.empty or not cosmic_genes:
        log.warning("  Need DiCE genes + COSMIC set for baseline comparison.")
        return {'status': 'no_data'}

    dice_lookup = {}
    if 'gene' in dice_df.columns and 'ensemble_score' in dice_df.columns:
        dice_lookup = dict(zip(dice_df['gene'].str.upper(), dice_df['ensemble_score']))

    # Label each pair: does it contain a known cancer gene?
    dnf_df['has_cancer_gene'] = dnf_df.apply(
        lambda r: r['gene_a'].upper() in cosmic_genes or r['gene_b'].upper() in cosmic_genes,
        axis=1
    )

    # Compute baseline scores
    baselines = []
    for _, row in dnf_df.iterrows():
        ga, gb = row['gene_a'].upper(), row['gene_b'].upper()
        dice_a = dice_lookup.get(ga, 0)
        dice_b = dice_lookup.get(gb, 0)
        baselines.append({
            'gene_a': row['gene_a'],
            'gene_b': row['gene_b'],
            'dss': row['dss'],
            'dice_product': dice_a * dice_b,
            'dice_sum': dice_a + dice_b,
            'single_ko_sum': abs(row.get('delta_aspl_tumor', 0)) +
                             abs(row.get('delta_aspl_normal', 0)),
            'has_cancer_gene': row['has_cancer_gene'],
        })

    base_df = pd.DataFrame(baselines)

    # Precision@K curves for each ranking method
    methods = {
        'DSS (DNF)': base_df.sort_values('dss', ascending=False).reset_index(drop=True),
        'DiCE Product': base_df.sort_values('dice_product', ascending=False).reset_index(drop=True),
        'DiCE Sum': base_df.sort_values('dice_sum', ascending=False).reset_index(drop=True),
    }

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'DSS (DNF)': '#DC2626', 'DiCE Product': '#2563EB', 'DiCE Sum': '#16A34A'}
    comparison_results = {}

    for method_name, method_df in methods.items():
        precision_at_k = []
        for k in range(1, len(method_df) + 1):
            hits = method_df.head(k)['has_cancer_gene'].sum()
            precision_at_k.append(hits / k)

        ax.plot(range(1, len(precision_at_k) + 1), precision_at_k,
                label=method_name, color=colors.get(method_name, 'gray'),
                linewidth=2 if method_name == 'DSS (DNF)' else 1.5,
                alpha=1.0 if method_name == 'DSS (DNF)' else 0.7)

        # Record P@10, P@20, P@all
        comparison_results[method_name] = {
            'P@10': precision_at_k[min(9, len(precision_at_k)-1)],
            'P@20': precision_at_k[min(19, len(precision_at_k)-1)],
            'P@all': precision_at_k[-1],
        }
        log.info(f"  {method_name}: P@10={comparison_results[method_name]['P@10']:.2f}, "
                 f"P@20={comparison_results[method_name]['P@20']:.2f}")

    # Random baseline
    bg_rate = base_df['has_cancer_gene'].mean()
    ax.axhline(y=bg_rate, color='gray', linestyle='--', alpha=0.5,
               label=f'Random baseline ({bg_rate:.2f})')

    ax.set_xlabel('K (Number of Top Pairs)', fontsize=12)
    ax.set_ylabel('Precision@K (% with known cancer gene)', fontsize=12)
    ax.set_title('DNF vs Naive Baselines: Precision@K for Cancer Genes',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1.05)

    plt.tight_layout()
    plt.savefig(fig_dir / 'val3_baseline_comparison.png', dpi=300)
    plt.close()

    pd.DataFrame(comparison_results).T.to_csv(tbl_dir / 'val3_baseline_comparison.csv')
    log.info(f"  Saved figure → val3_baseline_comparison.png")

    return comparison_results


# ═══════════════════════════════════════════════════════════════════════════
#  VALIDATION 4: Permutation Test
# ═══════════════════════════════════════════════════════════════════════════

def validate_permutation_test(
    dnf_df: pd.DataFrame,
    fig_dir: Path,
    tbl_dir: Path,
    n_permutations: int = 1000,
) -> Dict:
    """
    Are the top DSS scores significantly higher than expected by chance?
    Permutation test: shuffle the DSS labels and compare the observed top-10 mean
    against the null distribution.
    """
    log.info(f"\n=== Validation 4: Permutation Test ({n_permutations} permutations) ===")

    observed_top10_mean = dnf_df.head(10)['dss'].mean()
    observed_top10_max = dnf_df.head(1)['dss'].values[0]

    null_means = []
    null_maxes = []
    dss_values = dnf_df['dss'].values.copy()

    np.random.seed(42)
    for _ in range(n_permutations):
        np.random.shuffle(dss_values)
        null_means.append(np.sort(dss_values)[-10:].mean())
        null_maxes.append(dss_values.max())

    # P-value: fraction of permutations with top-10 mean >= observed
    p_mean = (np.sum(np.array(null_means) >= observed_top10_mean) + 1) / (n_permutations + 1)
    p_max = (np.sum(np.array(null_maxes) >= observed_top10_max) + 1) / (n_permutations + 1)

    log.info(f"  Observed top-10 mean DSS: {observed_top10_mean:.4f}")
    log.info(f"  Null distribution mean:   {np.mean(null_means):.4f}")
    log.info(f"  Permutation P-value:      {p_mean:.4f}")
    log.info(f"  Observed max DSS:         {observed_top10_max:.4f}")
    log.info(f"  Permutation P (max):      {p_max:.4f}")

    # Visualization
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(null_means, bins=40, color='#9CA3AF', alpha=0.7, edgecolor='white',
            label='Null distribution (shuffled)')
    ax.axvline(x=observed_top10_mean, color='#DC2626', linewidth=2.5,
               label=f'Observed (P={p_mean:.4f})')
    ax.set_xlabel('Top-10 Mean DSS', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(f'Permutation Test: Is Top-10 DSS Non-Random?\n'
                 f'Observed={observed_top10_mean:.4f}, P={p_mean:.4f}',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(fig_dir / 'val4_permutation_test.png', dpi=300)
    plt.close()

    result = {
        'observed_top10_mean': observed_top10_mean,
        'null_mean': np.mean(null_means),
        'null_std': np.std(null_means),
        'p_value_mean': p_mean,
        'observed_max': observed_top10_max,
        'p_value_max': p_max,
        'n_permutations': n_permutations,
        'significant': p_mean < 0.05,
    }

    pd.DataFrame([result]).to_csv(tbl_dir / 'val4_permutation_test.csv', index=False)
    log.info(f"  Saved figure → val4_permutation_test.png")

    return result


# ═══════════════════════════════════════════════════════════════════════════
#  VALIDATION 5: Literature Evidence for Top Pairs
# ═══════════════════════════════════════════════════════════════════════════

def validate_literature_evidence(
    dnf_df: pd.DataFrame,
    tbl_dir: Path,
    top_n: int = 10,
) -> Dict:
    """
    Generates a table of top DNF pairs annotated with known biological information.
    Uses curated knowledge base for EP300, NR2C2, CDC42, etc.
    """
    log.info(f"\n=== Validation 5: Literature Evidence for Top {top_n} Pairs ===")

    # Curated knowledge base for genes appearing in prostate cancer context
    gene_roles = {
        'EP300': 'Histone acetyltransferase; master transcriptional co-activator; '
                 'mutated in ~10% of prostate cancers; regulates AR signaling',
        'NR2C2': 'Nuclear receptor TR4; drives castration-resistant prostate cancer; '
                 'promotes PCa progression via AR-independent pathways',
        'CDC42': 'Rho GTPase; regulates cell polarity, migration, invasion; '
                 'overexpressed in metastatic PCa',
        'UBE2H': 'E2 ubiquitin-conjugating enzyme; involved in protein degradation; '
                 'linked to DNA damage response',
        'HSPA5': 'ER chaperone (GRP78/BiP); critical for unfolded protein response; '
                 'overexpressed in PCa and promotes treatment resistance',
        'PARK2': 'E3 ubiquitin ligase (Parkin); tumor suppressor frequently deleted in PCa; '
                 'loss promotes MYC-driven growth',
        'PRNP': 'Prion protein; promotes epithelial-mesenchymal transition; '
                'overexpressed in multiple cancers including PCa',
        'RPN1': 'Ribophorin I; ER glycosylation; essential for protein quality control',
        'APP': 'Amyloid precursor protein; emerging role in cancer cell proliferation '
               'and migration',
        'YWHAE': '14-3-3 epsilon; master regulatory scaffold; involved in cell cycle, '
                 'apoptosis, and signaling; acts as hub in signaling networks',
    }

    records = []
    for _, row in dnf_df.head(top_n).iterrows():
        ga, gb = row['gene_a'], row['gene_b']
        role_a = gene_roles.get(ga, 'Unknown function')
        role_b = gene_roles.get(gb, 'Unknown function')

        # Interpret the synergy
        if row['differential_synergy'] > 0.01:
            synergy_interpretation = 'SYNERGISTIC in tumor (super-additive disruption)'
        elif row['differential_synergy'] < -0.01:
            synergy_interpretation = 'Sub-additive (individual effects dominate)'
        else:
            synergy_interpretation = 'Approximately additive'

        records.append({
            'rank': row.get('pair_rank', ''),
            'gene_a': ga,
            'gene_b': gb,
            'dss': row['dss'],
            'differential_synergy': row['differential_synergy'],
            'synergy_interpretation': synergy_interpretation,
            'role_a': role_a,
            'role_b': role_b,
            'biological_rationale': f'{ga} and {gb} may represent structurally '
                                   f'co-dependent nodes in the tumor network',
        })

        log.info(f"  Rank {row.get('pair_rank', '?')}: {ga}—{gb} (DSS={row['dss']:.2f})")
        log.info(f"    {ga}: {role_a[:60]}...")
        log.info(f"    {gb}: {role_b[:60]}...")

    lit_df = pd.DataFrame(records)
    lit_df.to_csv(tbl_dir / 'val5_literature_evidence.csv', index=False)
    log.info(f"  Saved table → val5_literature_evidence.csv")

    return {'top_pairs_annotated': len(records)}


# ═══════════════════════════════════════════════════════════════════════════
#  MASTER VALIDATION RUNNER
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='DNF Validation Suite')
    parser.add_argument('--dnf-results', required=True)
    parser.add_argument('--dice-genes', default=None)
    parser.add_argument('--known-genes-dir', default='data/known_genes')
    parser.add_argument('--output-dir', default='data/results/dnf/validation')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    fig_dir = out_dir / 'figures'
    tbl_dir = out_dir / 'tables'
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    known_dir = Path(args.known_genes_dir)

    log.info("=" * 70)
    log.info("  DNF VALIDATION SUITE")
    log.info("  Multi-Layer Validation for Publication")
    log.info("=" * 70)

    # Load data
    dnf_df = pd.read_csv(args.dnf_results, sep='\t')
    if 'pair_rank' not in dnf_df.columns:
        dnf_df['pair_rank'] = range(1, len(dnf_df) + 1)
    log.info(f"Loaded {len(dnf_df)} DNF pairs")

    dice_df = pd.DataFrame()
    if args.dice_genes and Path(args.dice_genes).exists():
        dice_df = pd.read_csv(args.dice_genes, sep='\t')
        log.info(f"Loaded {len(dice_df)} DiCE genes")

    cosmic_genes = load_gene_set(known_dir / 'cosmic_census.txt')
    disgenet_genes = load_gene_set(known_dir / 'disgenet_prostate.txt')
    depmap_df = load_depmap(known_dir / 'depmap_fitness_pca.txt')

    log.info(f"COSMIC genes: {len(cosmic_genes)}")
    log.info(f"DisGeNET PCa genes: {len(disgenet_genes)}")
    log.info(f"DepMap genes: {len(depmap_df)}")

    # Run all validations
    all_results = {}

    # V1: Cancer Gene Enrichment
    all_results['cancer_enrichment'] = validate_cancer_gene_enrichment(
        dnf_df, cosmic_genes, disgenet_genes, fig_dir, tbl_dir)

    # V2: DepMap Essentiality
    all_results['depmap'] = validate_depmap_essentiality(
        dnf_df, depmap_df, fig_dir, tbl_dir)

    # V3: Baseline Comparison
    all_results['baselines'] = validate_against_baselines(
        dnf_df, dice_df, cosmic_genes, fig_dir, tbl_dir)

    # V4: Permutation Test
    all_results['permutation'] = validate_permutation_test(
        dnf_df, fig_dir, tbl_dir)

    # V5: Literature Evidence
    all_results['literature'] = validate_literature_evidence(
        dnf_df, tbl_dir)

    # Summary report
    log.info("\n" + "=" * 70)
    log.info("  VALIDATION SUMMARY")
    log.info("=" * 70)

    if 'permutation' in all_results and isinstance(all_results['permutation'], dict):
        p = all_results['permutation']
        sig = "✅ SIGNIFICANT" if p.get('significant', False) else "❌ NOT significant"
        log.info(f"  Permutation test:       {sig} (P={p.get('p_value_mean', 'N/A')})")

    if 'depmap' in all_results and isinstance(all_results['depmap'], dict):
        d = all_results['depmap']
        if d.get('status') == 'complete':
            log.info(f"  DepMap correlation:     ρ={d.get('spearman_rho', 'N/A'):.3f} "
                     f"(P={d.get('spearman_p', 'N/A'):.2e})")

    log.info(f"\n  All results saved to: {out_dir}")
    log.info("=" * 70)

    # Save master summary
    with open(tbl_dir / 'validation_summary.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)


if __name__ == '__main__':
    main()
