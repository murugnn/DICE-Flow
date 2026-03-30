#!/usr/bin/env python3
"""
Post-Analysis Utilities
Enrichment analysis and visualization for DiCE results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DiCEVisualizer:
    """Visualization tools for DiCE results"""
    
    def __init__(self, dice_file: str, output_dir: str = 'figures'):
        self.dice_df = pd.read_csv(dice_file, sep='\t')
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        sns.set_style("whitegrid")
        sns.set_palette("husl")
    
    def plot_top_genes(self, n: int = 20):
        """Bar plot of top N DiCE genes"""
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        top_n = self.dice_df.head(n)
        
        ax.barh(range(len(top_n)), top_n['ensemble_score'])
        ax.set_yticks(range(len(top_n)))
        ax.set_yticklabels(top_n['gene'])
        ax.set_xlabel('Ensemble Score', fontsize=12)
        ax.set_title(f'Top {n} DiCE Genes', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        
        plt.tight_layout()
        output = self.output_dir / f'top_{n}_genes.png'
        plt.savefig(output, dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot: {output}")
        plt.close()
    
    def plot_centrality_changes(self):
        """Scatter plot of betweenness vs eigenvector changes"""
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Color by ensemble score
        scatter = ax.scatter(
            self.dice_df['delta_betweenness'],
            self.dice_df['delta_eigenvector'],
            c=self.dice_df['ensemble_score'],
            cmap='viridis',
            alpha=0.6,
            s=50
        )
        
        # Annotate top 10 genes
        top_10 = self.dice_df.head(10)
        for _, row in top_10.iterrows():
            ax.annotate(
                row['gene'],
                (row['delta_betweenness'], row['delta_eigenvector']),
                fontsize=8,
                xytext=(5, 5),
                textcoords='offset points'
            )
        
        ax.set_xlabel('Δ Betweenness Centrality', fontsize=12)
        ax.set_ylabel('Δ Eigenvector Centrality', fontsize=12)
        ax.set_title('Centrality Changes: Tumor vs Normal', 
                    fontsize=14, fontweight='bold')
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Ensemble Score', fontsize=10)
        
        plt.tight_layout()
        output = self.output_dir / 'centrality_scatter.png'
        plt.savefig(output, dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot: {output}")
        plt.close()
    
    def plot_rank_distribution(self):
        """Distribution of ensemble scores"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Histogram
        ax1.hist(self.dice_df['ensemble_score'], bins=50, 
                edgecolor='black', alpha=0.7)
        ax1.set_xlabel('Ensemble Score', fontsize=12)
        ax1.set_ylabel('Frequency', fontsize=12)
        ax1.set_title('Distribution of Ensemble Scores', 
                     fontsize=14, fontweight='bold')
        ax1.axvline(self.dice_df['ensemble_score'].mean(), 
                   color='red', linestyle='--', label='Mean')
        ax1.legend()
        
        # Rank plot
        ax2.plot(range(1, len(self.dice_df) + 1), 
                self.dice_df['ensemble_score'].values)
        ax2.set_xlabel('Rank', fontsize=12)
        ax2.set_ylabel('Ensemble Score', fontsize=12)
        ax2.set_title('Ensemble Score by Rank', 
                     fontsize=14, fontweight='bold')
        ax2.set_yscale('log')
        
        plt.tight_layout()
        output = self.output_dir / 'score_distribution.png'
        plt.savefig(output, dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot: {output}")
        plt.close()
    
    def plot_heatmap(self, n: int = 50):
        """Heatmap of centrality measures for top genes"""
        
        top_n = self.dice_df.head(n)
        
        # Prepare data for heatmap
        data = top_n[['betweenness_norm', 'betweenness_tumor',
                      'eigenvector_norm', 'eigenvector_tumor']].T
        data.columns = top_n['gene']
        
        fig, ax = plt.subplots(figsize=(16, 4))
        
        sns.heatmap(data, cmap='RdBu_r', center=0, 
                   cbar_kws={'label': 'Centrality Value'},
                   linewidths=0.5, ax=ax)
        
        ax.set_ylabel('Centrality Measure', fontsize=12)
        ax.set_xlabel('Gene', fontsize=12)
        ax.set_title(f'Centrality Heatmap: Top {n} DiCE Genes',
                    fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        output = self.output_dir / f'centrality_heatmap_top{n}.png'
        plt.savefig(output, dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot: {output}")
        plt.close()
    
    def generate_all_plots(self):
        """Generate all standard plots"""
        logger.info("Generating visualizations...")
        
        self.plot_top_genes(20)
        self.plot_centrality_changes()
        self.plot_rank_distribution()
        self.plot_heatmap(50)
        
        logger.info(f"All plots saved to {self.output_dir}/")


class EnrichmentAnalyzer:
    """Gene set enrichment analysis placeholder"""
    
    def __init__(self, dice_file: str):
        self.dice_df = pd.read_csv(dice_file, sep='\t')
    
    def export_for_enrichr(self, output_file: str, top_n: int = 500):
        """Export gene list for Enrichr analysis"""
        
        top_genes = self.dice_df.head(top_n)['gene'].tolist()
        
        with open(output_file, 'w') as f:
            for gene in top_genes:
                f.write(f"{gene}\n")
        
        logger.info(f"Exported {len(top_genes)} genes for Enrichr: {output_file}")
        logger.info("Upload to: https://maayanlab.cloud/Enrichr/")
    
    def export_for_david(self, output_file: str, top_n: int = 500):
        """Export gene list for DAVID analysis"""
        
        self.export_for_enrichr(output_file, top_n)
        logger.info("Use this list at: https://david.ncifcrf.gov/")
    
    def export_for_gprofiler(self, output_file: str, top_n: int = 500):
        """Export gene list for g:Profiler"""
        
        self.export_for_enrichr(output_file, top_n)
        logger.info("Use this list at: https://biit.cs.ut.ee/gprofiler/")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='DiCE post-analysis utilities')
    parser.add_argument('--dice-file', required=True,
                       help='DiCE genes file (dice_genes_top500.txt)')
    parser.add_argument('--output-dir', default='figures',
                       help='Output directory for figures')
    parser.add_argument('--visualize', action='store_true',
                       help='Generate all visualizations')
    parser.add_argument('--export-enrichr', type=str,
                       help='Export gene list for Enrichr')
    parser.add_argument('--top-n', type=int, default=500,
                       help='Number of top genes to export')
    
    args = parser.parse_args()
    
    # Visualization
    if args.visualize:
        viz = DiCEVisualizer(args.dice_file, args.output_dir)
        viz.generate_all_plots()
    
    # Export for enrichment
    if args.export_enrichr:
        enrichment = EnrichmentAnalyzer(args.dice_file)
        enrichment.export_for_enrichr(args.export_enrichr, args.top_n)
        enrichment.export_for_david(
            args.export_enrichr.replace('.txt', '_david.txt'), 
            args.top_n
        )


if __name__ == '__main__':
    main()
