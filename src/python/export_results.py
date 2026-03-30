"""
Visualization Module for Cytoscape Export
Fixed: Correct Labels and Margins for Gene Names
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

class NetworkVisualizer:
    """Create visualizations and Cytoscape files"""
    
    def __init__(self, data_dir: str = "data"):
        self.data_dir = Path(data_dir)
        self.processed_dir = self.data_dir / "processed"
        self.results_dir = self.data_dir / "results"
        
        # Output directories
        self.cytoscape_dir = self.results_dir / "cytoscape"
        self.figures_dir = self.results_dir / "figures"
        
        self.cytoscape_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
    
    def prepare_cytoscape_network(self, network_file: str, output_name: str = "network"):
        logger.info(f"Preparing {output_name} for Cytoscape...")
        network_path = self.processed_dir / "weighted" / network_file
        try:
            df = pd.read_csv(network_path, sep='\t', header=None, names=['Source', 'Target', 'Weight'])
            output_path = self.cytoscape_dir / f"{output_name}_edges.csv"
            df.to_csv(output_path, index=False)
            logger.info(f"Saved Cytoscape network to {output_path}")
            return output_path
        except Exception as e:
            logger.error(f"Failed to load network {network_path}: {e}")
            return None
    
    def prepare_node_attributes(self, dice_df: pd.DataFrame, output_name: str = "nodes"):
        logger.info("Preparing node attributes for Cytoscape...")
        node_attrs = dice_df.copy()
            
        rename_map = {
            'gene': 'GeneName',
            'ensemble_score': 'EnsembleScore',
            'final_rank': 'DiCERank',
            'delta_betweenness': 'DeltaBetweenness',
            'delta_eigenvector': 'DeltaEigenvector',
            'betweenness_norm': 'Betweenness_Normal',
            'betweenness_tumor': 'Betweenness_Tumor',
        }
        
        # Only rename columns that actually exist
        available_cols = {k: v for k, v in rename_map.items() if k in node_attrs.columns}
        node_attrs = node_attrs.rename(columns=available_cols)
        
        # Add Regulation Direction
        if 'Betweenness_Tumor' in node_attrs.columns and 'Betweenness_Normal' in node_attrs.columns:
            node_attrs['Regulation'] = np.where(
                node_attrs['Betweenness_Tumor'] > node_attrs['Betweenness_Normal'], 
                'Increased_Centrality', 
                'Decreased_Centrality'
            )
        
        output_path = self.cytoscape_dir / f"{output_name}_attributes.csv"
        node_attrs.to_csv(output_path, index=False)
        logger.info(f"Saved node attributes to {output_path}")
        return output_path
    
    def plot_dice_rankings(self, dice_df: pd.DataFrame, top_n: int = 20):
        """Plot top DiCE genes with correct margins"""
        # Sort descending so top rank (1) is at the top of the bar chart
        top_genes = dice_df.head(top_n).iloc[::-1] 
        
        plt.figure(figsize=(12, 8)) # Wider figure
        
        # Plot bars
        plt.barh(top_genes['gene'], top_genes['ensemble_score'], color='teal')
        
        # Labeling
        plt.xlabel('DiCE Importance Score (Higher is Better)', fontsize=12)
        plt.title(f'Top {top_n} Rewired Cancer Genes', fontsize=14)
        plt.grid(axis='x', alpha=0.3)
        
        # Increase font size for gene names
        plt.yticks(fontsize=10)
        
        # CRITICAL FIX: Adjust margins so names are not cut off
        plt.tight_layout()
        
        output_path = self.figures_dir / f"top_{top_n}_dice_genes.png"
        plt.savefig(output_path, dpi=300)
        logger.info(f"Saved ranking plot to {output_path}")
        plt.close()

    def plot_centrality_scatter(self, dice_df: pd.DataFrame):
        """Scatter plot of Normal vs Tumor Betweenness"""
        if 'betweenness_norm' not in dice_df.columns:
            logger.warning("Missing raw centrality columns, skipping scatter plot.")
            return

        plt.figure(figsize=(8, 8))
        
        eps = 1e-9
        x = np.log10(dice_df['betweenness_norm'] + eps)
        y = np.log10(dice_df['betweenness_tumor'] + eps)
        
        sns.scatterplot(x=x, y=y, alpha=0.3, edgecolor=None, color='grey')
        
        # Highlight top 10
        top = dice_df.head(10)
        top_x = np.log10(top['betweenness_norm'] + eps)
        top_y = np.log10(top['betweenness_tumor'] + eps)
        
        sns.scatterplot(x=top_x, y=top_y, color='red', s=60, label='Top Candidates')
        
        # Add labels with better offset
        for i, row in top.iterrows():
            plt.text(np.log10(row['betweenness_norm']+eps)+0.1, 
                     np.log10(row['betweenness_tumor']+eps), 
                     row['gene'], fontsize=9, color='darkred')

        plt.plot([x.min(), x.max()], [x.min(), x.max()], 'k--', alpha=0.5)
        plt.title("Differential Centrality (Rewiring Analysis)")
        plt.xlabel("Log10 Betweenness (Normal)")
        plt.ylabel("Log10 Betweenness (Tumor)")
        plt.legend()
        plt.tight_layout()
        
        output_path = self.figures_dir / "centrality_scatter.png"
        plt.savefig(output_path, dpi=300)
        logger.info(f"Saved scatter plot to {output_path}")
        plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', default='data')
    args = parser.parse_args()
    
    viz = NetworkVisualizer(args.data_dir)
    
    # Load DiCE Results
    dice_path = viz.results_dir / "dice_genes.txt"
    if not dice_path.exists():
        logger.error(f"Results not found at {dice_path}")
        return
        
    dice_df = pd.read_csv(dice_path, sep='\t')
    
    # Export for Cytoscape
    viz.prepare_cytoscape_network("network_tumor.txt", "network_tumor")
    viz.prepare_node_attributes(dice_df, "node_attributes")
    
    # Generate Figures
    viz.plot_dice_rankings(dice_df)
    viz.plot_centrality_scatter(dice_df)
    
    logger.info("Visualization complete.")

if __name__ == "__main__":
    main()