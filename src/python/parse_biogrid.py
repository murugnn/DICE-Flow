import pandas as pd
import argparse
import logging
from pathlib import Path
import sys

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

class BioGRIDParser:
    def __init__(self, input_path, output_path):
        self.input_path = Path(input_path)
        self.output_path = Path(output_path)
        
    def parse(self, organism_id=9606):
        """
        Parse BioGRID TAB2 file.
        Default organism_id 9606 is Homo sapiens.
        """
        if not self.input_path.exists():
            logger.error(f"Input file not found: {self.input_path}")
            sys.exit(1)

        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Parsing BioGRID from {self.input_path}")
        
        # BioGRID Tab2 Columns (0-based index):
        # 0: BioGRID Interaction ID
        # 1: Entrez Gene Interactor A
        # 2: Entrez Gene Interactor B
        # 7: Official Symbol Interactor A
        # 8: Official Symbol Interactor B
        # 11: Experimental System
        # 12: Experimental System Type
        # 15: Organism ID Interactor A
        # 16: Organism ID Interactor B
        
        # We will load by index to avoid header name issues
        use_cols = [0, 1, 2, 7, 8, 11, 12, 15, 16]
        col_names = [
            'interaction_id', 'entrez_a', 'entrez_b', 
            'symbol_a', 'symbol_b', 
            'system', 'system_type', 
            'org_a', 'org_b'
        ]

        interactions = []
        chunk_size = 100000
        
        try:
            # Skip the header row (row 0) and use our own names
            # header=0 tells pandas row 0 is header, but we replace it with `names`
            # However, BioGRID usually has a header line. We'll rely on skipping it if it exists.
            
            # First, detect if header exists and starts with #
            with open(self.input_path, 'r') as f:
                first_line = f.readline()
                header_row = 0 if first_line.startswith('#') else None

            processed_count = 0
            
            # Read in chunks
            with pd.read_csv(self.input_path, sep='\t', chunksize=chunk_size, 
                           usecols=use_cols, names=col_names, header=header_row,
                           low_memory=False, on_bad_lines='skip') as reader:
                
                for chunk in reader:
                    # Normalize organism IDs to integers (handle NaNs)
                    chunk['org_a'] = pd.to_numeric(chunk['org_a'], errors='coerce').fillna(0).astype(int)
                    chunk['org_b'] = pd.to_numeric(chunk['org_b'], errors='coerce').fillna(0).astype(int)
                    
                    # Filter for specific organism (both interactors must belong to it)
                    mask = (chunk['org_a'] == organism_id) & (chunk['org_b'] == organism_id)
                    
                    # Filter for physical interactions only (optional, but recommended for signaling)
                    # Common physical types: Affinity Capture-Western, Two-hybrid, Co-crystal structure, etc.
                    # We usually keep everything or filter 'genetic'. Let's keep all for now 
                    # or strictly physical if desired.
                    # mask &= (chunk['system_type'] == 'physical') 
                    
                    filtered = chunk[mask].copy()
                    
                    if not filtered.empty:
                        # Select relevant columns for output
                        # We want: Source, Target, Type, Weight (default 1.0)
                        out_df = pd.DataFrame({
                            'source': filtered['symbol_a'],
                            'target': filtered['symbol_b'],
                            'interaction': filtered['system'],
                            'type': filtered['system_type']
                        })
                        interactions.append(out_df)
                    
                    processed_count += len(chunk)
                    print(f"Processed {processed_count} rows...", end='\r')

            if not interactions:
                logger.warning("No interactions found for the specified organism.")
                return pd.DataFrame()

            final_df = pd.concat(interactions, ignore_index=True)
            
            # Remove self-loops and duplicates
            final_df = final_df[final_df['source'] != final_df['target']]
            final_df.drop_duplicates(subset=['source', 'target'], inplace=True)
            
            # Save
            final_df.to_csv(self.output_path, sep='\t', index=False)
            logger.info(f"\nSaved {len(final_df)} unique interactions to {self.output_path}")
            
            return final_df

        except Exception as e:
            logger.error(f"Failed to parse BioGRID: {e}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Parse BioGRID TAB2 file')
    parser.add_argument('--input', required=True, help='Path to raw BioGRID tab2 file')
    parser.add_argument('--output', required=True, help='Path to output processed file')
    parser.add_argument('--organism', type=int, default=9606, help='NCBI Taxonomy ID (default: 9606 for human)')
    
    args = parser.parse_args()
    
    parser = BioGRIDParser(args.input, args.output)
    parser.parse(organism_id=args.organism)

if __name__ == '__main__':
    main()