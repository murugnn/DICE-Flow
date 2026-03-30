#!/usr/bin/env python3
"""
Network Merger - Updated with DiCE Candidate Filtering
"""

import pandas as pd
import logging
import argparse
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class NetworkMerger:
    def __init__(self, kegg_file: str, biogrid_file: str, candidate_file: str = None):
        self.kegg_file = Path(kegg_file)
        self.biogrid_file = Path(biogrid_file)
        self.candidate_set = None
        
        if candidate_file:
            logger.info(f"Loading candidate genes from {candidate_file}")
            with open(candidate_file, 'r') as f:
                self.candidate_set = set(line.strip() for line in f if line.strip())
            logger.info(f"Loaded {len(self.candidate_set)} candidate genes for filtering.")

    def merge(self) -> pd.DataFrame:
        logger.info("Loading KEGG interactions...")
        kegg_df = pd.read_csv(self.kegg_file, sep='\t')
        
        logger.info("Loading BioGRID interactions...")
        biogrid_df = pd.read_csv(self.biogrid_file, sep='\t')
        
        # --- FILTERING LOGIC ---
        if self.candidate_set:
            k_len_before = len(kegg_df)
            kegg_df = kegg_df[kegg_df['source'].isin(self.candidate_set) & 
                              kegg_df['target'].isin(self.candidate_set)]
            logger.info(f"Filtered KEGG: {k_len_before} -> {len(kegg_df)} edges")

            b_len_before = len(biogrid_df)
            biogrid_df = biogrid_df[biogrid_df['source'].isin(self.candidate_set) & 
                                    biogrid_df['target'].isin(self.candidate_set)]
            logger.info(f"Filtered BioGRID: {b_len_before} -> {len(biogrid_df)} edges")
        # -----------------------

        kegg_edges = set()
        for _, row in kegg_df.iterrows():
            kegg_edges.add((row['source'], row['target']))
        
        merged_edges = []
        
        # Add KEGG (Directed)
        for edge in kegg_edges:
            merged_edges.append({'source': edge[0], 'target': edge[1], 'source_db': 'KEGG', 'directed': True})
        
        # Add BioGRID (Bidirectional if not in KEGG)
        biogrid_added = 0
        for _, row in biogrid_df.iterrows():
            gene_a, gene_b = row['source'], row['target']
            
            if (gene_a, gene_b) not in kegg_edges and (gene_b, gene_a) not in kegg_edges:
                merged_edges.append({'source': gene_a, 'target': gene_b, 'source_db': 'BioGRID', 'directed': False})
                merged_edges.append({'source': gene_b, 'target': gene_a, 'source_db': 'BioGRID', 'directed': False})
                biogrid_added += 2
        
        merged_df = pd.DataFrame(merged_edges)
        if not merged_df.empty:
            merged_df = merged_df[merged_df['source'] != merged_df['target']]
            merged_df = merged_df.drop_duplicates(subset=['source', 'target'])
        
        return merged_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kegg', required=True)
    parser.add_argument('--biogrid', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--candidates', help="Optional path to filtered gene list")
    
    args = parser.parse_args()
    
    merger = NetworkMerger(args.kegg, args.biogrid, args.candidates)
    merged_df = merger.merge()
    
    merged_df.to_csv(args.output, index=False, sep='\t')
    logger.info(f"Saved merged network with {len(merged_df)} edges to {args.output}")

if __name__ == '__main__':
    main()