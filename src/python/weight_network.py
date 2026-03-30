#!/usr/bin/env python3
"""
Network Weighting using Gene Expression
Implements DiCE weighting: Weight = 1 - |Pearson correlation|
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import logging
from pathlib import Path
from typing import Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class NetworkWeighter:
    """Weight network edges using gene expression correlations"""
    
    def __init__(self, network_file: str, expression_file: str, 
                 scale_factor: int = 1000):
        """
        Args:
            network_file: Merged network file (source, target)
            expression_file: Gene expression matrix (genes x samples)
            scale_factor: Multiply weights by this for integer conversion
        """
        self.network_file = Path(network_file)
        self.expression_file = Path(expression_file)
        self.scale_factor = scale_factor
        self.epsilon = 1e-6  # Minimum weight to avoid zero
    
    def weight_network(self) -> pd.DataFrame:
        """Generate weighted network"""
        
        # Load network
        logger.info(f"Loading network from {self.network_file}")
        network_df = pd.read_csv(self.network_file, sep='\t')
        logger.info(f"Loaded {len(network_df)} edges")
        
        # Load expression data
        logger.info(f"Loading expression data from {self.expression_file}")
        expr_df = self._load_expression()
        logger.info(f"Loaded expression for {len(expr_df)} genes across {expr_df.shape[1]} samples")
        
        # Get available genes
        available_genes = set(expr_df.index)
        
        # Filter network to only genes with expression data
        network_df['has_source'] = network_df['source'].isin(available_genes)
        network_df['has_target'] = network_df['target'].isin(available_genes)
        
        initial_count = len(network_df)
        network_df = network_df[network_df['has_source'] & network_df['has_target']]
        logger.info(f"Filtered to {len(network_df)} edges with expression data "
                   f"({initial_count - len(network_df)} removed)")
        
        # Calculate weights
        logger.info("Calculating edge weights...")
        weights = []
        
        for idx, row in network_df.iterrows():
            if idx % 10000 == 0:
                logger.info(f"Processed {idx}/{len(network_df)} edges")
            
            source = row['source']
            target = row['target']
            
            weight = self._calculate_weight(expr_df, source, target)
            weights.append(weight)
        
        network_df['weight'] = weights
        
        # Convert to integers for Duan algorithm
        network_df['weight_int'] = (network_df['weight'] * self.scale_factor).astype(int)
        
        # Ensure no zero weights
        network_df.loc[network_df['weight_int'] == 0, 'weight_int'] = 1
        
        logger.info(f"Weight statistics:")
        logger.info(f"  Mean: {network_df['weight'].mean():.4f}")
        logger.info(f"  Std: {network_df['weight'].std():.4f}")
        logger.info(f"  Min: {network_df['weight'].min():.4f}")
        logger.info(f"  Max: {network_df['weight'].max():.4f}")
        
        return network_df[['source', 'target', 'weight', 'weight_int']]
    
    def _load_expression(self) -> pd.DataFrame:
        """Load expression matrix"""
        
        # Try to detect file format
        file_ext = self.expression_file.suffix.lower()
        
        if file_ext == '.csv':
            df = pd.read_csv(self.expression_file, index_col=0)
        elif file_ext in ['.txt', '.tsv', '.tab']:
            df = pd.read_csv(self.expression_file, sep='\t', index_col=0)
        elif file_ext in ['.gct']:
            # GCT format (skip first two lines)
            df = pd.read_csv(self.expression_file, sep='\t', skiprows=2, index_col=1)
            df = df.iloc[:, 1:]  # Remove Description column
        else:
            # Try tab-separated as default
            df = pd.read_csv(self.expression_file, sep='\t', index_col=0)
        
        # Clean gene names (remove version numbers like .1, .2)
        df.index = df.index.str.split('.').str[0]
        
        return df
    
    def _calculate_weight(self, expr_df: pd.DataFrame, 
                         gene1: str, gene2: str) -> float:
        """
        Calculate DiCE weight: 1 - |Pearson correlation|
        
        High correlation -> Low weight (short distance)
        Low correlation -> High weight (long distance)
        """
        
        try:
            expr1 = expr_df.loc[gene1].values
            expr2 = expr_df.loc[gene2].values
            
            # Calculate Pearson correlation
            corr, _ = pearsonr(expr1, expr2)
            
            # Handle NaN correlations
            if np.isnan(corr):
                corr = 0.0
            
            # DiCE formula: Weight = 1 - |r|
            weight = 1.0 - abs(corr)
            
            # Add epsilon to avoid zero weights
            weight = max(weight, self.epsilon)
            
            return weight
            
        except Exception as e:
            logger.warning(f"Error calculating weight for {gene1}-{gene2}: {e}")
            return 1.0  # Default weight


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Weight network using gene expression correlations'
    )
    parser.add_argument('--network', required=True, 
                       help='Merged network file')
    parser.add_argument('--expression', required=True, 
                       help='Gene expression matrix file')
    parser.add_argument('--output', required=True, 
                       help='Output weighted network file')
    parser.add_argument('--scale', type=int, default=1000,
                       help='Scale factor for integer weights (default: 1000)')
    
    args = parser.parse_args()
    
    # Weight network
    weighter = NetworkWeighter(args.network, args.expression, args.scale)
    weighted_df = weighter.weight_network()
    
    # Save
    weighted_df.to_csv(args.output, index=False, sep='\t')
    logger.info(f"Saved weighted network to {args.output}")
    
    # Also save integer-only version for C++ engine
    int_output = args.output.replace('.txt', '_int.txt')
    weighted_df[['source', 'target', 'weight_int']].to_csv(
        int_output, index=False, sep='\t', header=False
    )
    logger.info(f"Saved integer-weighted network to {int_output}")


if __name__ == '__main__':
    main()
