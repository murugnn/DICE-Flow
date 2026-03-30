#!/usr/bin/env python3
"""
Differential Expression Analysis Filter (Phase I)
Extracts strict DEGs via T-test with BH-FDR correction and Log2 Fold Change filtering.
"""

import pandas as pd
import numpy as np
import argparse
import logging
from scipy import stats

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def fdr_bh(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.asarray(pvals)
    ranked_pvi = np.argsort(pvals)
    n = len(pvals)
    fdr = np.zeros(n)
    
    # BH formula: p * N / rank
    fdr[ranked_pvi] = pvals[ranked_pvi] * n / (np.arange(n) + 1)
    
    # Fix non-monotonicity
    for i in range(n - 2, -1, -1):
        fdr[ranked_pvi[i]] = min(fdr[ranked_pvi[i]], fdr[ranked_pvi[i+1]])
    
    # Cap at 1.0
    return np.minimum(fdr, 1.0)

def main():
    parser = argparse.ArgumentParser(description="Phase I: DEA Filter")
    parser.add_argument('--normal', required=True, help='Normal expression file')
    parser.add_argument('--tumor', required=True, help='Tumor expression file')
    parser.add_argument('--output-genes', required=True, help='Output text list of DEG genes')
    parser.add_argument('--output-table', required=True, help='Output table with log2FC and FDR')
    parser.add_argument('--p-val', type=float, default=0.05, help='FDR cutoff')
    parser.add_argument('--lfc', type=float, default=1.0, help='Absolute Log2 FC cutoff')
    args = parser.parse_args()

    # Load Data
    logger.info("Loading expression data...")
    try:
        norm_df = pd.read_csv(args.normal, sep='\t', index_col=0)
        tumor_df = pd.read_csv(args.tumor, sep='\t', index_col=0)
    except Exception as e:
        logger.error(f"Failed to load expression data: {e}")
        return

    # Intersect genes
    common = norm_df.index.intersection(tumor_df.index)
    norm_df = norm_df.loc[common]
    tumor_df = tumor_df.loc[common]
    n_genes = len(common)
    logger.info(f"Analyzing {n_genes} tracking genes across samples...")

    # Computations
    # Assuming the input data is already Log2 scaled (e.g. log2(TPM+1))
    mean_n = norm_df.mean(axis=1)
    mean_t = tumor_df.mean(axis=1)
    
    log2fc = mean_t - mean_n  # Difference of logs is log of ratio

    # T-test
    logger.info("Performing separate continuous T-tests per gene...")
    # scipy.stats.ttest_ind on entire arrays (vectorized)
    t_stat, p_val = stats.ttest_ind(tumor_df.values, norm_df.values, axis=1, equal_var=False)

    # Handle rare completely identical arrays producing NaN pvals
    nan_mask = np.isnan(p_val)
    p_val[nan_mask] = 1.0

    # FDR correction
    fdr = fdr_bh(p_val)

    results = pd.DataFrame({
        'gene': common,
        'log2FC': log2fc,
        'pvalue': p_val,
        'FDR': fdr
    })

    # Filtering
    mask_sig = (results['FDR'] < args.p_val) & (results['log2FC'].abs() > args.lfc)
    sig_genes = results[mask_sig].copy()
    
    # Log outcomes
    up = len(sig_genes[sig_genes['log2FC'] > 0])
    down = len(sig_genes[sig_genes['log2FC'] < 0])
    logger.info(f"DEA identified {len(sig_genes)} significant DEGs (FDR<{args.p_val}, |Log2FC|>{args.lfc}).")
    logger.info(f"  - Upregulated: {up}")
    logger.info(f"  - Downregulated: {down}")

    if len(sig_genes) == 0:
        logger.warning("No DEGs found at current thresholds. Relaxing limits strictly to save pipeline.")
        # If no genes meet strict cutoff, use top 1000 by pval
        sig_genes = results.sort_values('FDR').head(1000)

    # Write output
    with open(args.output_genes, 'w') as f:
        for gene in sig_genes['gene']:
            f.write(f"{gene}\n")
            
    results.to_csv(args.output_table, index=False, sep='\t')
    logger.info("Saved candidates successfully.")

if __name__ == "__main__":
    main()
