import pandas as pd
import numpy as np
import argparse
import logging
from sklearn.feature_selection import mutual_info_classif

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="DiCE Phase II: Information Gain Filter")
    parser.add_argument('--normal', required=True, help='Normal expression file')
    parser.add_argument('--tumor', required=True, help='Tumor expression file')
    parser.add_argument('--output', required=True, help='Output file for kept gene list')
    args = parser.parse_args()

    # 1. Load Data
    logger.info("Loading expression data...")
    try:
        norm_df = pd.read_csv(args.normal, sep='\t', index_col=0)
        tumor_df = pd.read_csv(args.tumor, sep='\t', index_col=0)
    except Exception as e:
        logger.error(f"Failed to load expression data: {e}")
        return

    # Ensure genes match
    common_genes = norm_df.index.intersection(tumor_df.index)
    norm_df = norm_df.loc[common_genes]
    tumor_df = tumor_df.loc[common_genes]
    
    logger.info(f"Analyzing {len(common_genes)} common genes...")

    # 2. Prepare for Scikit-Learn (Features=Genes, Samples=Rows)
    # Transpose: rows are samples, cols are genes
    X_norm = norm_df.T
    X_tumor = tumor_df.T
    
    X = pd.concat([X_norm, X_tumor])
    
    # Create Labels: 0 for Normal, 1 for Tumor
    y = np.concatenate([np.zeros(len(X_norm)), np.ones(len(X_tumor))])

    # 3. Calculate Information Gain (Mutual Information)
    logger.info("Calculating Information Gain (this may take a moment)...")
    ig_scores = mutual_info_classif(X, y, discrete_features=False, random_state=42)

    # 4. Filter Logic (Keep genes > Mean IG)
    # DiCE Paper Phase II: "Genes are selected if IG > Average IG"
    threshold = np.mean(ig_scores)
    logger.info(f"IG Threshold (Mean): {threshold:.6f}")

    results = pd.DataFrame({'gene': common_genes, 'ig': ig_scores})
    kept_genes = results[results['ig'] > threshold]['gene'].tolist()

    logger.info(f"Filtering complete: Kept {len(kept_genes)} / {len(common_genes)} genes.")

    # 5. Save
    with open(args.output, 'w') as f:
        for gene in kept_genes:
            f.write(f"{gene}\n")
    logger.info(f"Saved candidate gene list to {args.output}")

if __name__ == "__main__":
    main()