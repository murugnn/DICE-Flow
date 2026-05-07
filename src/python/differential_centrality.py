#!/usr/bin/env python3
"""
Differential Centrality Analysis
Implements DiCE Phase V: Ensemble Ranking based on Absolute Differences
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
import argparse
from scipy.stats import norm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DifferentialCentralityAnalyzer:
    def __init__(self, normal_centrality: str, tumor_centrality: str):
        self.normal_file = Path(normal_centrality)
        self.tumor_file = Path(tumor_centrality)
    
    def analyze(self) -> pd.DataFrame:
        logger.info("Loading centrality data...")
        try:
            normal_df = pd.read_csv(self.normal_file, sep='\t')
            tumor_df = pd.read_csv(self.tumor_file, sep='\t')
        except Exception as e:
            logger.error(f"Error reading files: {e}")
            return pd.DataFrame()

        merged = pd.merge(normal_df, tumor_df, on='gene', suffixes=('_norm', '_tumor'))
        merged['diff_betweenness'] = (merged['betweenness_tumor'] - merged['betweenness_norm']).abs()
        merged['diff_eigenvector'] = (merged['eigenvector_tumor'] - merged['eigenvector_norm']).abs()
        merged['rank_bet'] = merged['diff_betweenness'].rank(pct=True)
        merged['rank_eig'] = merged['diff_eigenvector'].rank(pct=True)

        # Stouffer combined score: rank-percentiles -> z-scores -> sum/sqrt(2).
        # Decorrelated alternative to rank-product (which double-counts when
        # betweenness and eigenvector are correlated, as they typically are).
        # Clip away from {0, 1} since norm.ppf maps them to +/-infinity.
        eps = 1e-6
        z_bet = norm.ppf(merged['rank_bet'].clip(eps, 1 - eps))
        z_eig = norm.ppf(merged['rank_eig'].clip(eps, 1 - eps))
        merged['ensemble_score_stouffer'] = (z_bet + z_eig) / (2 ** 0.5)

        # Always compute the raw ensemble (no hub filtering)
        merged['ensemble_score_raw'] = merged['rank_bet'] * merged['rank_eig']

        # Hub-corrected: same score but applied AFTER dropping top-12% eigenvector
        # AND keeping genes above-average betweenness in at least one condition.
        eig_n = merged['eigenvector_norm'].quantile(0.88)
        eig_t = merged['eigenvector_tumor'].quantile(0.88)
        avg_bet_n = merged['betweenness_norm'].mean()
        avg_bet_t = merged['betweenness_tumor'].mean()
        mask_keep = ((merged['eigenvector_norm'] <= eig_n) & (merged['eigenvector_tumor'] <= eig_t)) & \
                    ((merged['betweenness_norm'] > avg_bet_n) | (merged['betweenness_tumor'] > avg_bet_t))

        hub_df = merged[mask_keep].copy()
        hub_df['ensemble_score_hub_corrected'] = hub_df['rank_bet'] * hub_df['rank_eig']
        merged = merged.merge(
            hub_df[['gene', 'ensemble_score_hub_corrected']],
            on='gene', how='left'
        )

        if getattr(self, 'dice_faithful', False):
            out = hub_df.copy()
            out['ensemble_score'] = out['ensemble_score_hub_corrected']
            out = out.drop(columns=['ensemble_score_hub_corrected'])
            out = out.sort_values('ensemble_score', ascending=False).copy()
            logger.warning(
                "--dice-faithful: dropped %d top-12%% eigenvector hubs (incl. potentially TP53/EP300)",
                len(merged) - len(out),
            )
        else:
            out = merged.sort_values('ensemble_score_raw', ascending=False).copy()
            logger.info(
                "Default mode: kept all %d genes; reporting raw and hub-corrected scores",
                len(out),
            )

        out['final_rank'] = range(1, len(out) + 1)
        out['delta_betweenness'] = out['betweenness_tumor'] - out['betweenness_norm']
        out['delta_eigenvector'] = out['eigenvector_tumor'] - out['eigenvector_norm']

        score_col = 'ensemble_score' if getattr(self, 'dice_faithful', False) else 'ensemble_score_raw'
        logger.info("Ranking complete. Top gene: %s (Score: %.4f)",
                    out.iloc[0]['gene'], out.iloc[0][score_col])

        return out

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--normal', required=True)
    parser.add_argument('--tumor', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--top-n', type=int, default=500)
    parser.add_argument('--dice-faithful', action='store_true',
                        help='Reproduce Pashaei et al. behavior bit-exactly: enable top-12%% '
                             'eigenvector hub removal, drop Stouffer score column.')
    parser.add_argument('--disable-noise', action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.disable_noise:
        args.dice_faithful = False  # legacy: --disable-noise was the inverse

    analyzer = DifferentialCentralityAnalyzer(args.normal, args.tumor)
    analyzer.dice_faithful = bool(args.dice_faithful)
        
    dice_df = analyzer.analyze()
    
    if dice_df.empty:
        logger.error("Analysis produced no results.")
        return

    # Save full
    dice_df.to_csv(args.output, index=False, sep='\t')
    
    # Save top N
    top_output = args.output.replace('.txt', f'_top{args.top_n}.txt')
    dice_df.head(args.top_n).to_csv(top_output, index=False, sep='\t')
    logger.info(f"Saved Top {args.top_n} to {top_output}")

if __name__ == '__main__':
    main()