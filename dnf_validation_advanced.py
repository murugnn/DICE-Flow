#!/usr/bin/env python3
"""
DNF Advanced Validation Suite — Pair-Level Analysis
=====================================================
Addresses the three critical gaps in DNF validation:

  V6: Pair-level SL validation (SynLethDB API + BioGRID physical interactions)
  V7: Pathway / functional coherence (KEGG co-pathway analysis)
  V8: Tumor specificity statistical proof (paired test on ΔASPL)

Usage:
    python dnf_validation_advanced.py \
        --dnf-results data/results/dnf_pairwise_results.txt \
        --biogrid data/processed/networks/biogrid_interactions.txt \
        --kegg-dir data/raw/kegg \
        --output-dir data/results/dnf/validation
"""

import argparse
import logging
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
log = logging.getLogger(__name__)

plt.rcParams.update({
    'figure.dpi': 300, 'savefig.dpi': 300, 'font.size': 11,
    'axes.titlesize': 13, 'axes.labelsize': 12,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.15,
})


# ═══════════════════════════════════════════════════════════════════════════
#  V6: PAIR-LEVEL SYNTHETIC LETHALITY VALIDATION
# ═══════════════════════════════════════════════════════════════════════════

def build_biogrid_interaction_set(biogrid_path: Path) -> Set[frozenset]:
    """Build a set of known physical interactions from BioGRID."""
    interactions = set()
    if not biogrid_path.exists():
        log.warning(f"BioGRID file not found: {biogrid_path}")
        return interactions

    df = pd.read_csv(biogrid_path, sep='\t', nrows=None)
    for _, row in df.iterrows():
        src = str(row.get('source', '')).upper().strip()
        tgt = str(row.get('target', '')).upper().strip()
        if src and tgt:
            interactions.add(frozenset([src, tgt]))

    log.info(f"  BioGRID: Loaded {len(interactions)} unique physical interactions")
    return interactions


def build_synlethdb_set(synlethdb_path: Path) -> Tuple[Set[frozenset], pd.DataFrame]:
    """Load SynLethDB gene_sl_gene.tsv and build a set of known SL pairs."""
    sl_pairs = set()
    if not synlethdb_path.exists():
        log.warning(f"SynLethDB file not found: {synlethdb_path}")
        return sl_pairs, pd.DataFrame()

    df = pd.read_csv(synlethdb_path, sep='\t', low_memory=False)
    # Columns: x_name, y_name are the gene symbols
    for _, row in df.iterrows():
        ga = str(row.get('x_name', '')).upper().strip()
        gb = str(row.get('y_name', '')).upper().strip()
        if ga and gb and ga != 'NAN' and gb != 'NAN':
            sl_pairs.add(frozenset([ga, gb]))

    log.info(f"  SynLethDB: Loaded {len(sl_pairs)} unique SL pairs from {len(df)} records")
    return sl_pairs, df


def validate_synlethdb(
    dnf_df: pd.DataFrame,
    sl_pairs: Set[frozenset],
    sl_raw_df: pd.DataFrame,
    fig_dir: Path,
    tbl_dir: Path,
    top_k: int = 50,
) -> Dict:
    """
    V6b: Cross-reference DNF predictions against SynLethDB known SL pairs.

    Key questions:
    1. Do top DNF pairs appear in SynLethDB? (direct overlap)
    2. Is the overlap statistically enriched vs background? (hypergeometric test)
    3. What are the novel predictions? (pairs NOT in SynLethDB)
    """
    log.info("\n=== V6b: SynLethDB Synthetic Lethality Validation ===")

    if not sl_pairs:
        log.warning("  No SynLethDB data available.")
        return {'status': 'no_data'}

    # Collect all unique genes in our DNF analysis
    dnf_genes = set(dnf_df['gene_a'].str.upper().tolist() +
                    dnf_df['gene_b'].str.upper().tolist())

    # Count how many SL pairs exist in SynLethDB for our gene universe
    # (important for background rate calculation)
    sl_in_universe = set()
    for pair in sl_pairs:
        genes_in_pair = pair & dnf_genes
        if len(genes_in_pair) >= 1:  # At least one gene in our set
            sl_in_universe.add(pair)

    sl_both_in_universe = set()
    for pair in sl_pairs:
        if pair <= dnf_genes:  # Both genes in our set
            sl_both_in_universe.add(pair)

    log.info(f"  Our gene universe: {len(dnf_genes)} genes")
    log.info(f"  SL pairs with ≥1 gene in universe: {len(sl_in_universe)}")
    log.info(f"  SL pairs with BOTH genes in universe: {len(sl_both_in_universe)}")

    # Label all DNF pairs
    dnf_df['synlethdb_hit'] = dnf_df.apply(
        lambda r: frozenset([r['gene_a'].upper(), r['gene_b'].upper()]) in sl_pairs,
        axis=1
    )

    bg_rate = dnf_df['synlethdb_hit'].mean()
    total_hits = dnf_df['synlethdb_hit'].sum()
    log.info(f"  Total DNF pairs matching SynLethDB: {total_hits}/{len(dnf_df)} ({bg_rate*100:.1f}%)")

    # Test at multiple cutoffs
    cutoffs = [10, 20, 50, 100, len(dnf_df)]
    cutoff_results = []

    for k in cutoffs:
        if k > len(dnf_df):
            continue
        top = dnf_df.head(k)
        hits = top['synlethdb_hit'].sum()
        pct = hits / k * 100

        from scipy.stats import binomtest
        try:
            p_val = binomtest(int(hits), k, bg_rate, alternative='greater').pvalue
        except Exception:
            p_val = 1.0

        # Hypergeometric test (more appropriate for overlap)
        # M = total possible pairs, n = known SL pairs, N = our top-K, k = hits
        from scipy.stats import hypergeom
        M = len(dnf_df)  # Total pairs tested
        n_sl = total_hits  # Total SL pairs in our data
        hyper_p = hypergeom.sf(int(hits) - 1, M, n_sl, k) if hits > 0 else 1.0

        cutoff_results.append({
            'top_k': k,
            'sl_hits': int(hits),
            'total': k,
            'hit_rate_pct': round(pct, 1),
            'background_pct': round(bg_rate * 100, 1),
            'fold_enrichment': round(pct / (bg_rate * 100), 2) if bg_rate > 0 else 0,
            'binomial_p': p_val,
            'hypergeom_p': hyper_p,
        })

        marker = "✅" if hyper_p < 0.05 else "  "
        log.info(f"  {marker} Top-{k}: {hits}/{k} known SL pairs "
                 f"({pct:.1f}% vs {bg_rate*100:.1f}% bg, "
                 f"Binom P={p_val:.2e}, Hyper P={hyper_p:.2e})")

    pd.DataFrame(cutoff_results).to_csv(tbl_dir / 'val6b_synlethdb_enrichment.csv', index=False)

    # List the actual SL hits
    sl_hit_pairs = dnf_df[dnf_df['synlethdb_hit']].copy()
    if len(sl_hit_pairs) > 0:
        log.info(f"\n  Known SL pairs found in DNF results:")
        for _, row in sl_hit_pairs.head(20).iterrows():
            log.info(f"    Rank {row.get('pair_rank', '?')}: "
                     f"{row['gene_a']}—{row['gene_b']} (DSS={row['dss']:+.4f})")

        sl_hit_pairs[['gene_a', 'gene_b', 'dss', 'synergy_tumor',
                      'differential_synergy', 'pair_rank']].to_csv(
            tbl_dir / 'val6b_synlethdb_hits.csv', index=False)

    # Novel predictions (not in SynLethDB)
    top_novel = dnf_df.head(top_k)
    novel = top_novel[~top_novel['synlethdb_hit']]
    log.info(f"\n  Top-{top_k}: {len(novel)} novel predicted SL pairs (not in SynLethDB)")
    novel[['gene_a', 'gene_b', 'dss', 'synergy_tumor',
           'differential_synergy']].to_csv(
        tbl_dir / 'val6b_novel_sl_predictions.csv', index=False)

    # Compare DSS: SL-verified vs non-SL
    sl_dss = dnf_df[dnf_df['synlethdb_hit']]['dss']
    non_sl_dss = dnf_df[~dnf_df['synlethdb_hit']]['dss']
    mw_p = np.nan
    if len(sl_dss) >= 3 and len(non_sl_dss) >= 3:
        _, mw_p = stats.mannwhitneyu(sl_dss, non_sl_dss, alternative='two-sided')
        log.info(f"\n  DSS comparison:")
        log.info(f"    Known SL mean DSS:     {sl_dss.mean():.4f} (n={len(sl_dss)})")
        log.info(f"    Non-SL mean DSS:       {non_sl_dss.mean():.4f} (n={len(non_sl_dss)})")
        log.info(f"    Mann-Whitney P:        {mw_p:.4f}")

    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # (A) Enrichment at cutoffs
    ax = axes[0]
    df_plot = pd.DataFrame(cutoff_results)
    x = np.arange(len(df_plot))
    w = 0.35
    ax.bar(x - w/2, df_plot['hit_rate_pct'], w, label='DNF Top-K',
           color='#DC2626', alpha=0.8, edgecolor='white')
    ax.bar(x + w/2, df_plot['background_pct'], w, label='Background',
           color='#9CA3AF', alpha=0.6, edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels([f"Top-{r}" for r in df_plot['top_k']], fontsize=9)
    ax.set_ylabel('% Pairs in SynLethDB')
    ax.set_title('SynLethDB Synthetic Lethality Enrichment', fontweight='bold')
    ax.legend(fontsize=9)

    # (B) DSS distribution: SL vs non-SL
    ax = axes[1]
    if len(sl_dss) > 0:
        ax.hist(non_sl_dss, bins=40, alpha=0.5, color='#9CA3AF',
                edgecolor='white', label=f'Non-SL (n={len(non_sl_dss)})', density=True)
        ax.hist(sl_dss, bins=min(20, max(5, len(sl_dss)//3)), alpha=0.7, color='#DC2626',
                edgecolor='white', label=f'Known SL (n={len(sl_dss)})', density=True)
        ax.legend(fontsize=9)
    ax.set_xlabel('DSS')
    ax.set_ylabel('Density')
    p_label = f"MW P={mw_p:.4f}" if not np.isnan(mw_p) else ""
    ax.set_title(f'DSS: Known SL Pairs vs Others\n{p_label}', fontweight='bold')

    # (C) Cumulative SL hits across ranking
    ax = axes[2]
    cumulative_hits = dnf_df['synlethdb_hit'].cumsum().values
    x_ranks = np.arange(1, len(dnf_df) + 1)
    # Expected by chance
    expected = bg_rate * x_ranks
    ax.plot(x_ranks, cumulative_hits, color='#DC2626', linewidth=2, label='Observed SL hits')
    ax.plot(x_ranks, expected, color='gray', linewidth=1.5, linestyle='--', label='Expected (random)')
    ax.fill_between(x_ranks, expected, cumulative_hits,
                     where=cumulative_hits > expected, alpha=0.15, color='red')
    ax.set_xlabel('Pair Rank (by DSS, descending)')
    ax.set_ylabel('Cumulative SL Hits')
    ax.set_title('Cumulative SL Recovery Curve', fontweight='bold')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(fig_dir / 'val6b_synlethdb.png', dpi=300)
    plt.close()
    log.info(f"  Saved figure → val6b_synlethdb.png")

    return {
        'total_sl_in_db': len(sl_pairs),
        'total_hits_in_dnf': int(total_hits),
        'background_rate': bg_rate,
        'cutoff_results': cutoff_results,
        'n_novel_top50': int(len(novel)),
        'mannwhitney_p': float(mw_p) if not np.isnan(mw_p) else None,
        'sl_mean_dss': float(sl_dss.mean()) if len(sl_dss) > 0 else None,
        'non_sl_mean_dss': float(non_sl_dss.mean()) if len(non_sl_dss) > 0 else None,
    }


def validate_pair_level(
    dnf_df: pd.DataFrame,
    biogrid_interactions: Set[frozenset],
    fig_dir: Path,
    tbl_dir: Path,
    top_k: int = 50,
) -> Dict:
    """
    V6: Do top DNF pairs represent known physical interactions or predicted novel pairs?

    Three checks:
    1. BioGRID physical interaction overlap
    2. Gene pair co-occurrence in the same pathway (via KEGG, handled separately)
    3. Framing of non-overlapping pairs as novel predictions
    """
    log.info("\n=== V6: Pair-Level Interaction Validation ===")

    # Check BioGRID overlap at multiple cutoffs
    cutoffs = [10, 20, 50, 100, len(dnf_df)]
    cutoff_results = []

    # Label all pairs
    dnf_df['biogrid_interaction'] = dnf_df.apply(
        lambda r: frozenset([r['gene_a'].upper(), r['gene_b'].upper()]) in biogrid_interactions,
        axis=1
    )

    bg_rate = dnf_df['biogrid_interaction'].mean()
    log.info(f"  BioGRID background interaction rate: {bg_rate*100:.1f}%")

    for k in cutoffs:
        if k > len(dnf_df):
            continue
        top = dnf_df.head(k)
        hits = top['biogrid_interaction'].sum()
        pct = hits / k * 100

        from scipy.stats import binomtest
        try:
            p_val = binomtest(hits, k, bg_rate, alternative='greater').pvalue
        except Exception:
            p_val = 1.0

        cutoff_results.append({
            'top_k': k,
            'interacting_pairs': hits,
            'total': k,
            'hit_rate_pct': round(pct, 1),
            'background_pct': round(bg_rate * 100, 1),
            'fold_enrichment': round(pct / (bg_rate * 100), 2) if bg_rate > 0 else 0,
            'pvalue': p_val,
        })

        marker = "✅" if p_val < 0.05 else "  "
        log.info(f"  {marker} Top-{k}: {hits}/{k} known interactions "
                 f"({pct:.1f}% vs {bg_rate*100:.1f}% bg, P={p_val:.2e})")

    pd.DataFrame(cutoff_results).to_csv(tbl_dir / 'val6_biogrid_pair_validation.csv', index=False)

    # Identify novel predicted pairs (not in BioGRID)
    top_novel = dnf_df.head(top_k)
    novel_pairs = top_novel[~top_novel['biogrid_interaction']].copy()
    known_pairs = top_novel[top_novel['biogrid_interaction']].copy()

    log.info(f"\n  Top-{top_k} pairs breakdown:")
    log.info(f"    Known interactions (BioGRID): {len(known_pairs)}")
    log.info(f"    Novel predicted pairs:        {len(novel_pairs)}")

    # Save pair tables
    novel_pairs[['gene_a', 'gene_b', 'dss', 'synergy_tumor', 'differential_synergy']].to_csv(
        tbl_dir / 'val6_novel_predicted_pairs.csv', index=False
    )
    known_pairs[['gene_a', 'gene_b', 'dss', 'synergy_tumor', 'differential_synergy']].to_csv(
        tbl_dir / 'val6_known_interaction_pairs.csv', index=False
    )

    # Compare DSS of interacting vs non-interacting pairs
    interacting = dnf_df[dnf_df['biogrid_interaction']]['dss']
    non_interacting = dnf_df[~dnf_df['biogrid_interaction']]['dss']

    if len(interacting) > 5 and len(non_interacting) > 5:
        u_stat, mw_p = stats.mannwhitneyu(interacting, non_interacting, alternative='two-sided')
        log.info(f"\n  DSS of interacting vs non-interacting pairs:")
        log.info(f"    Interacting mean DSS:     {interacting.mean():.4f} (n={len(interacting)})")
        log.info(f"    Non-interacting mean DSS: {non_interacting.mean():.4f} (n={len(non_interacting)})")
        log.info(f"    Mann-Whitney P:           {mw_p:.4f}")
    else:
        mw_p = np.nan

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Bar chart: enrichment at cutoffs
    ax = axes[0]
    df_plot = pd.DataFrame(cutoff_results)
    x = np.arange(len(df_plot))
    w = 0.35
    ax.bar(x - w/2, df_plot['hit_rate_pct'], w, label='DNF Top-K',
           color='#7C3AED', alpha=0.8, edgecolor='white')
    ax.bar(x + w/2, df_plot['background_pct'], w, label='Background',
           color='#9CA3AF', alpha=0.6, edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels([f"Top-{r}" for r in df_plot['top_k']])
    ax.set_ylabel('% Pairs with Known Interaction')
    ax.set_title('BioGRID Physical Interaction Enrichment', fontweight='bold')
    ax.legend()

    # Box plot: DSS of interacting vs non-interacting
    ax = axes[1]
    data = [interacting.values, non_interacting.values]
    bp = ax.boxplot(data, tick_labels=['Known\nInteraction', 'No Known\nInteraction'],
                    patch_artist=True, widths=0.5)
    bp['boxes'][0].set_facecolor('#7C3AED')
    bp['boxes'][0].set_alpha(0.5)
    bp['boxes'][1].set_facecolor('#9CA3AF')
    bp['boxes'][1].set_alpha(0.5)
    ax.set_ylabel('Differential Synergy Score (DSS)')
    p_label = f"P={mw_p:.4f}" if not np.isnan(mw_p) else "N/A"
    ax.set_title(f'DSS by Interaction Status ({p_label})', fontweight='bold')

    plt.tight_layout()
    plt.savefig(fig_dir / 'val6_pair_interaction.png', dpi=300)
    plt.close()
    log.info(f"  Saved figure → val6_pair_interaction.png")

    return {
        'cutoff_results': cutoff_results,
        'n_novel': len(novel_pairs),
        'n_known': len(known_pairs),
        'mannwhitney_p': float(mw_p) if not np.isnan(mw_p) else None,
    }


# ═══════════════════════════════════════════════════════════════════════════
#  V7: PATHWAY / FUNCTIONAL COHERENCE
# ═══════════════════════════════════════════════════════════════════════════

def parse_kegg_pathways(kegg_dir: Path) -> Dict[str, List[str]]:
    """
    Parse KEGG XML files to build a gene → pathway mapping.
    Returns dict: gene_symbol → list of pathway names
    """
    gene_to_pathways = defaultdict(list)
    pathway_to_genes = defaultdict(set)

    if not kegg_dir.exists():
        log.warning(f"KEGG directory not found: {kegg_dir}")
        return gene_to_pathways

    xml_files = list(kegg_dir.glob("*.xml"))
    log.info(f"  Parsing {len(xml_files)} KEGG pathway files...")

    for xml_file in xml_files:
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()

            # Get pathway name
            pathway_name = root.get('title', xml_file.stem)

            # Find gene entries
            for entry in root.findall('.//entry'):
                if entry.get('type') == 'gene':
                    # The 'name' attr contains KEGG IDs like "hsa:7157 hsa:7158"
                    # The 'graphics' child has the actual gene symbols
                    graphics = entry.find('graphics')
                    if graphics is not None:
                        gene_label = graphics.get('name', '')
                        # Symbols can be comma-separated: "TP53, TRP53"
                        for gene in gene_label.replace(',', ' ').replace('...', '').split():
                            gene = gene.strip().upper()
                            if gene and len(gene) > 1 and not gene.startswith('HSA:'):
                                gene_to_pathways[gene].append(pathway_name)
                                pathway_to_genes[pathway_name].add(gene)
        except Exception:
            continue

    log.info(f"  KEGG: {len(gene_to_pathways)} genes across {len(pathway_to_genes)} pathways")
    return gene_to_pathways, pathway_to_genes


def validate_pathway_coherence(
    dnf_df: pd.DataFrame,
    gene_to_pathways: Dict[str, List[str]],
    pathway_to_genes: Dict[str, Set[str]],
    fig_dir: Path,
    tbl_dir: Path,
    top_k: int = 50,
) -> Dict:
    """
    V7: Do top DNF gene pairs share functional pathways?

    For each pair (A, B): Find shared KEGG pathways.
    Compare top-K vs bottom-K pairs for pathway sharing rate.
    """
    log.info("\n=== V7: Pathway / Functional Coherence ===")

    if not gene_to_pathways:
        log.warning("  No KEGG pathway data available.")
        return {'status': 'no_data'}

    # For each pair, compute shared pathways
    pair_pathway_data = []
    for _, row in dnf_df.iterrows():
        ga = row['gene_a'].upper()
        gb = row['gene_b'].upper()

        pw_a = set(gene_to_pathways.get(ga, []))
        pw_b = set(gene_to_pathways.get(gb, []))
        shared = pw_a & pw_b

        pair_pathway_data.append({
            'gene_a': row['gene_a'],
            'gene_b': row['gene_b'],
            'dss': row['dss'],
            'pathways_a': len(pw_a),
            'pathways_b': len(pw_b),
            'shared_pathways': len(shared),
            'has_shared': len(shared) > 0,
            'shared_names': '; '.join(sorted(shared)[:5]),  # Top 5 shared
            'jaccard': len(shared) / len(pw_a | pw_b) if (pw_a | pw_b) else 0,
        })

    pw_df = pd.DataFrame(pair_pathway_data)

    # Compare top vs bottom
    top_pairs = pw_df.head(top_k)
    bottom_pairs = pw_df.tail(top_k)

    top_sharing = top_pairs['has_shared'].mean() * 100
    bottom_sharing = bottom_pairs['has_shared'].mean() * 100
    overall_sharing = pw_df['has_shared'].mean() * 100

    log.info(f"  Top-{top_k} pairs: {top_sharing:.1f}% share a pathway")
    log.info(f"  Bottom-{top_k} pairs: {bottom_sharing:.1f}% share a pathway")
    log.info(f"  Overall: {overall_sharing:.1f}% share a pathway")

    # Statistical test: Is sharing rate higher in top-K?
    from scipy.stats import fisher_exact
    top_shared = top_pairs['has_shared'].sum()
    top_not = len(top_pairs) - top_shared
    bot_shared = bottom_pairs['has_shared'].sum()
    bot_not = len(bottom_pairs) - bot_shared

    table = [[top_shared, top_not], [bot_shared, bot_not]]
    odds_ratio, fisher_p = fisher_exact(table, alternative='greater')
    log.info(f"  Fisher's exact (top > bottom): OR={odds_ratio:.2f}, P={fisher_p:.4f}")

    # Correlation: DSS vs number of shared pathways
    valid_pw = pw_df.dropna(subset=['dss', 'shared_pathways'])
    if len(valid_pw) > 5:
        rho, rho_p = stats.spearmanr(valid_pw['dss'], valid_pw['shared_pathways'])
        log.info(f"  Spearman ρ(DSS, shared_pathways) = {rho:.3f}, P = {rho_p:.2e}")
    else:
        rho, rho_p = 0, 1

    # Show top pairs with their shared pathways
    log.info(f"\n  Top-10 pairs — shared pathways:")
    for _, row in pw_df.head(10).iterrows():
        shared_str = row['shared_names'] if row['shared_names'] else "(none)"
        log.info(f"    {row['gene_a']:>10}—{row['gene_b']:<10} "
                 f"DSS={row['dss']:+.2f}  shared={row['shared_pathways']}  [{shared_str[:50]}]")

    # Save tables
    pw_df.to_csv(tbl_dir / 'val7_pathway_coherence.csv', index=False)
    pw_df.head(50).to_csv(tbl_dir / 'val7_top50_pathway_detail.csv', index=False)

    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1. Sharing rate comparison
    ax = axes[0]
    bars = ax.bar(['Top-50\n(High DSS)', 'Bottom-50\n(Low DSS)', 'Overall'],
                  [top_sharing, bottom_sharing, overall_sharing],
                  color=['#DC2626', '#2563EB', '#9CA3AF'], alpha=0.8,
                  edgecolor='white', width=0.5)
    ax.set_ylabel('% Pairs Sharing ≥1 KEGG Pathway')
    ax.set_title(f'Pathway Sharing: Top vs Bottom Pairs\nFisher P={fisher_p:.4f}',
                 fontweight='bold')
    # Add values on bars
    for bar, val in zip(bars, [top_sharing, bottom_sharing, overall_sharing]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{val:.0f}%', ha='center', fontweight='bold', fontsize=11)

    # 2. DSS vs shared pathways scatter
    ax = axes[1]
    jitter = np.random.normal(0, 0.15, len(pw_df))
    ax.scatter(pw_df['dss'], pw_df['shared_pathways'] + jitter,
               s=20, alpha=0.5, c=pw_df['dss'], cmap='RdBu_r', edgecolors='none')
    ax.set_xlabel('DSS')
    ax.set_ylabel('Number of Shared KEGG Pathways')
    ax.set_title(f'DSS vs Functional Similarity\nρ={rho:.3f}, P={rho_p:.2e}',
                 fontweight='bold')

    # 3. Jaccard similarity distribution
    ax = axes[2]
    top_jaccard = pw_df.head(top_k)['jaccard']
    bot_jaccard = pw_df.tail(top_k)['jaccard']
    ax.hist(top_jaccard, bins=20, alpha=0.7, color='#DC2626', label='Top-50 (High DSS)',
            edgecolor='white', density=True)
    ax.hist(bot_jaccard, bins=20, alpha=0.5, color='#2563EB', label='Bottom-50 (Low DSS)',
            edgecolor='white', density=True)
    ax.set_xlabel('Jaccard Pathway Similarity')
    ax.set_ylabel('Density')
    ax.set_title('Pathway Similarity Distribution', fontweight='bold')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(fig_dir / 'val7_pathway_coherence.png', dpi=300)
    plt.close()
    log.info(f"  Saved figure → val7_pathway_coherence.png")

    return {
        'top_k': top_k,
        'top_sharing_pct': top_sharing,
        'bottom_sharing_pct': bottom_sharing,
        'fisher_p': fisher_p,
        'odds_ratio': odds_ratio,
        'spearman_rho': rho,
        'spearman_p': rho_p,
    }


# ═══════════════════════════════════════════════════════════════════════════
#  V8: TUMOR SPECIFICITY STATISTICAL PROOF
# ═══════════════════════════════════════════════════════════════════════════

def validate_tumor_specificity(
    dnf_df: pd.DataFrame,
    fig_dir: Path,
    tbl_dir: Path,
) -> Dict:
    """
    V8: Prove that tumor network disruption is statistically greater than
    normal network disruption for high-DSS pairs.

    Tests:
    1. Paired Wilcoxon signed-rank test on ΔASPL_tumor vs ΔASPL_normal (all pairs)
    2. Same test on top-K pairs only
    3. Effect size (Cohen's d)
    4. Visualization: paired difference distribution
    """
    log.info("\n=== V8: Tumor Specificity Statistical Proof ===")

    d_tumor = dnf_df['delta_aspl_tumor'].values
    d_normal = dnf_df['delta_aspl_normal'].values
    diff = d_tumor - d_normal  # This IS the DSS

    # 1. All pairs: Paired tests
    wilcox_stat, wilcox_p = stats.wilcoxon(d_tumor, d_normal, alternative='greater')
    ttest_stat, ttest_p = stats.ttest_rel(d_tumor, d_normal)

    # Effect size: Cohen's d
    d_diff = d_tumor - d_normal
    cohens_d = d_diff.mean() / d_diff.std() if d_diff.std() > 0 else 0

    log.info(f"  --- ALL {len(dnf_df)} PAIRS ---")
    log.info(f"  Mean ΔASPL_tumor:  {d_tumor.mean():.4f}")
    log.info(f"  Mean ΔASPL_normal: {d_normal.mean():.4f}")
    log.info(f"  Mean DSS (diff):   {d_diff.mean():.4f}")
    log.info(f"  Wilcoxon signed-rank (tumor > normal): W={wilcox_stat:.0f}, P={wilcox_p:.2e}")
    log.info(f"  Paired t-test:     t={ttest_stat:.3f}, P={ttest_p:.2e}")
    log.info(f"  Cohen's d:         {cohens_d:.3f}")

    results = {
        'all_pairs': {
            'n': len(dnf_df),
            'mean_delta_tumor': d_tumor.mean(),
            'mean_delta_normal': d_normal.mean(),
            'mean_dss': d_diff.mean(),
            'wilcoxon_W': float(wilcox_stat),
            'wilcoxon_p': float(wilcox_p),
            'ttest_t': float(ttest_stat),
            'ttest_p': float(ttest_p),
            'cohens_d': float(cohens_d),
        }
    }

    # 2. Top-K analysis (the high-DSS pairs SHOULD be tumor-specific)
    for k in [10, 20, 50]:
        top = dnf_df.head(k)
        dt = top['delta_aspl_tumor'].values
        dn = top['delta_aspl_normal'].values

        w_s, w_p = stats.wilcoxon(dt, dn, alternative='greater')
        dd = dt - dn
        cd = dd.mean() / dd.std() if dd.std() > 0 else 0

        results[f'top_{k}'] = {
            'n': k,
            'mean_delta_tumor': dt.mean(),
            'mean_delta_normal': dn.mean(),
            'mean_dss': dd.mean(),
            'wilcoxon_p': float(w_p),
            'cohens_d': float(cd),
        }

        sig = "✅" if w_p < 0.05 else "❌"
        log.info(f"\n  --- TOP-{k} PAIRS ---")
        log.info(f"  {sig} Wilcoxon P={w_p:.2e}, Cohen's d={cd:.2f}")
        log.info(f"  Mean ΔASPL_tumor={dt.mean():.4f} vs ΔASPL_normal={dn.mean():.4f}")

    # 3. Fraction of positive DSS pairs (simple but important)
    n_positive = (d_diff > 0).sum()
    n_total = len(d_diff)
    binom_p = stats.binomtest(n_positive, n_total, 0.5, alternative='greater').pvalue
    log.info(f"\n  Pairs with ΔASPL_tumor > ΔASPL_normal: {n_positive}/{n_total} "
             f"({n_positive/n_total*100:.1f}%)")
    log.info(f"  Binomial test (>50%): P={binom_p:.2e}")

    results['fraction_tumor_selective'] = {
        'n_positive': int(n_positive),
        'n_total': n_total,
        'fraction': n_positive / n_total,
        'binomial_p': float(binom_p),
    }

    # 4. Visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # (A) Paired comparison: tumor vs normal disruption for top-50
    ax = axes[0, 0]
    top50 = dnf_df.head(50)
    x = np.arange(len(top50))
    ax.bar(x - 0.2, top50['delta_aspl_tumor'], 0.4, label='Tumor',
           color='#DC2626', alpha=0.7, edgecolor='white')
    ax.bar(x + 0.2, top50['delta_aspl_normal'], 0.4, label='Normal',
           color='#2563EB', alpha=0.7, edgecolor='white')
    ax.set_xlabel('Pair Rank (by DSS)')
    ax.set_ylabel('ΔASPL (Network Disruption)')
    ax.set_title('Top-50 Pairs: Tumor vs Normal Disruption', fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_xlim(-1, 51)

    # (B) DSS distribution with tumor-selective zone
    ax = axes[0, 1]
    ax.hist(d_diff, bins=50, color='#7C3AED', alpha=0.7, edgecolor='white')
    ax.axvline(x=0, color='black', linewidth=2, linestyle='--')
    ax.fill_betweenx([0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 50],
                      0, max(d_diff), alpha=0.1, color='red')
    ax.annotate(f'Tumor-selective\n{n_positive/n_total*100:.0f}% of pairs',
                xy=(max(d_diff)*0.5, ax.get_ylim()[1]*0.8 if ax.get_ylim()[1] > 0 else 30),
                fontsize=10, color='darkred', fontweight='bold', ha='center')
    ax.set_xlabel('DSS (ΔASPL_tumor − ΔASPL_normal)')
    ax.set_ylabel('Number of Pairs')
    ax.set_title(f'DSS Distribution\nBinomial P={binom_p:.2e}', fontweight='bold')

    # (C) Scatter: ΔASPL_tumor vs ΔASPL_normal
    ax = axes[1, 0]
    sc = ax.scatter(d_normal, d_tumor, c=d_diff, cmap='RdBu_r', s=15,
                    alpha=0.7, edgecolors='none')
    plt.colorbar(sc, ax=ax, label='DSS', shrink=0.8)
    # Diagonal
    lims = [min(d_normal.min(), d_tumor.min()),
            max(d_normal.max(), d_tumor.max())]
    ax.plot(lims, lims, 'k--', alpha=0.3, linewidth=1)
    ax.fill_between(lims, lims, [lims[1], lims[1]], alpha=0.05, color='red')
    ax.set_xlabel('ΔASPL Normal')
    ax.set_ylabel('ΔASPL Tumor')
    ax.set_title(f'Paired Disruption Comparison\nWilcoxon P={wilcox_p:.2e}', fontweight='bold')

    # Label top 3
    for _, row in dnf_df.head(3).iterrows():
        ax.annotate(f"{row['gene_a']}—{row['gene_b']}",
                    (row['delta_aspl_normal'], row['delta_aspl_tumor']),
                    fontsize=7, fontweight='bold', color='darkred',
                    textcoords="offset points", xytext=(5, 5))

    # (D) Effect size across ranks
    ax = axes[1, 1]
    window = 50
    rolling_dss = []
    ranks = []
    for start in range(0, len(dnf_df) - window + 1, 10):
        chunk = dnf_df.iloc[start:start + window]
        mean_dss = chunk['dss'].mean()
        rolling_dss.append(mean_dss)
        ranks.append(start + window // 2)

    ax.plot(ranks, rolling_dss, color='#DC2626', linewidth=2)
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    ax.fill_between(ranks, 0, rolling_dss, where=[d > 0 for d in rolling_dss],
                     alpha=0.15, color='red', label='Tumor-selective')
    ax.fill_between(ranks, 0, rolling_dss, where=[d <= 0 for d in rolling_dss],
                     alpha=0.15, color='blue', label='Normal-resilient')
    ax.set_xlabel('Pair Rank (rolling window=50)')
    ax.set_ylabel('Mean DSS')
    ax.set_title('DSS Across Ranking Spectrum', fontweight='bold')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(fig_dir / 'val8_tumor_specificity.png', dpi=300)
    plt.close()
    log.info(f"  Saved figure → val8_tumor_specificity.png")

    # Save numerical results
    pd.DataFrame([results['all_pairs']]).to_csv(tbl_dir / 'val8_tumor_specificity_all.csv', index=False)

    return results


# ═══════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='DNF Advanced Validation (V6-V8)')
    parser.add_argument('--dnf-results', required=True)
    parser.add_argument('--biogrid', default='data/processed/networks/biogrid_interactions.txt')
    parser.add_argument('--synlethdb', default=None, help='SynLethDB gene_sl_gene.tsv file')
    parser.add_argument('--kegg-dir', default='data/raw/kegg')
    parser.add_argument('--output-dir', default='data/results/dnf/validation')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    fig_dir = out_dir / 'figures'
    tbl_dir = out_dir / 'tables'
    fig_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 70)
    log.info("  DNF ADVANCED VALIDATION (V6-V8)")
    log.info("  Pair-Level · Pathway Coherence · Tumor Specificity")
    log.info("=" * 70)

    # Load DNF results
    dnf_df = pd.read_csv(args.dnf_results, sep='\t')
    if 'pair_rank' not in dnf_df.columns:
        dnf_df['pair_rank'] = range(1, len(dnf_df) + 1)
    log.info(f"Loaded {len(dnf_df)} DNF pairs")

    all_results = {}

    # V6: Pair-level validation (BioGRID)
    biogrid_interactions = build_biogrid_interaction_set(Path(args.biogrid))
    all_results['v6_pair_level'] = validate_pair_level(
        dnf_df, biogrid_interactions, fig_dir, tbl_dir)

    # V6b: SynLethDB validation
    if args.synlethdb:
        sl_pairs, sl_raw_df = build_synlethdb_set(Path(args.synlethdb))
        all_results['v6b_synlethdb'] = validate_synlethdb(
            dnf_df, sl_pairs, sl_raw_df, fig_dir, tbl_dir)
    else:
        log.info("\n  SynLethDB file not provided. Skipping V6b.")
        log.info("  Use --synlethdb gene_sl_gene.tsv to enable.")

    # V7: Pathway coherence
    kegg_dir = Path(args.kegg_dir)
    if kegg_dir.exists():
        gene_to_pathways, pathway_to_genes = parse_kegg_pathways(kegg_dir)
        all_results['v7_pathway'] = validate_pathway_coherence(
            dnf_df, gene_to_pathways, pathway_to_genes, fig_dir, tbl_dir)
    else:
        log.warning("KEGG directory not found. Skipping V7.")

    # V8: Tumor specificity proof
    all_results['v8_tumor_specificity'] = validate_tumor_specificity(
        dnf_df, fig_dir, tbl_dir)

    # Summary
    log.info("\n" + "=" * 70)
    log.info("  ADVANCED VALIDATION SUMMARY")
    log.info("=" * 70)

    v6b = all_results.get('v6b_synlethdb', {})
    if v6b and v6b.get('total_hits_in_dnf', 0) > 0:
        log.info(f"  SynLethDB hits:               {v6b['total_hits_in_dnf']}/{len(dnf_df)} pairs")
        if v6b.get('mannwhitney_p') is not None:
            log.info(f"  SL vs non-SL DSS (MW):        P={v6b['mannwhitney_p']:.4f}")

    v8 = all_results.get('v8_tumor_specificity', {}).get('all_pairs', {})
    if v8:
        log.info(f"  Tumor specificity (Wilcoxon): P={v8.get('wilcoxon_p', 'N/A')}")
        log.info(f"  Cohen's d:                    {v8.get('cohens_d', 'N/A')}")

    v7 = all_results.get('v7_pathway', {})
    if v7 and v7 != {'status': 'no_data'}:
        log.info(f"  Pathway sharing (Fisher):     P={v7.get('fisher_p', 'N/A')}")

    log.info(f"\n  All results saved to: {out_dir}")
    log.info("=" * 70)

    # Save master results
    with open(tbl_dir / 'advanced_validation_summary.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)


if __name__ == '__main__':
    main()
