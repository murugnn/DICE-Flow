#!/usr/bin/env python3
"""
KEGG XML Parser
Extracts directed interactions from KEGG pathway XML files
"""

import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
from typing import Set, Tuple
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class KEGGParser:
    """Parse KEGG KGML (XML) files to extract directed protein interactions"""
    
    def __init__(self, kegg_dir: str):
        self.kegg_dir = Path(kegg_dir)
        self.interactions = set()
        
    def parse_all_files(self) -> pd.DataFrame:
        """Parse all XML files in KEGG directory"""
        logger.info(f"Parsing KEGG files from {self.kegg_dir}")
        
        xml_files = list(self.kegg_dir.glob("*.xml"))
        logger.info(f"Found {len(xml_files)} XML files")
        
        for xml_file in xml_files:
            try:
                self._parse_file(xml_file)
            except Exception as e:
                logger.warning(f"Failed to parse {xml_file}: {e}")
        
        logger.info(f"Extracted {len(self.interactions)} unique directed interactions")
        
        # Convert to DataFrame
        df = pd.DataFrame(list(self.interactions), columns=['source', 'target'])
        df['direction'] = 'KEGG'
        
        return df
    
    def _parse_file(self, xml_file: Path):
        """Parse a single KGML file to extract Gene Symbols"""
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            # Build entry ID to Gene SYMBOL mapping
            entry_map = {}
            for entry in root.findall('entry'):
                entry_id = entry.get('id')
                entry_type = entry.get('type')
                
                if entry_type in ['gene', 'protein']:
                    # LOOK HERE: We grab the <graphics> child tag
                    graphics = entry.find('graphics')
                    if graphics is not None:
                        # Name is often "TP53, P53, LFS1..." -> take the first one
                        name_str = graphics.get('name', '')
                        if name_str:
                            # Split by comma, take first item, strip '...', uppercase
                            symbol = name_str.split(',')[0].strip('.').upper()
                            entry_map[entry_id] = [symbol]
            
            # Extract relations (directed interactions)
            for relation in root.findall('relation'):
                entry1 = relation.get('entry1')
                entry2 = relation.get('entry2')
                rel_type = relation.get('type')
                
                # Only keep activation, inhibition, expression, etc.
                if rel_type not in ['PPrel', 'GErel', 'PCrel']:
                    continue
                
                if entry1 in entry_map and entry2 in entry_map:
                    for gene1 in entry_map[entry1]:
                        for gene2 in entry_map[entry2]:
                            # Filter out self-loops immediately
                            if gene1 == gene2: 
                                continue
                                
                            subtypes = relation.findall('subtype')
                            
                            # Check subtypes for directionality logic
                            # (We store all as directed, but logic determines implication)
                            self.interactions.add((gene1, gene2))
                            
        except Exception as e:
            logger.warning(f"Error parsing {xml_file.name}: {e}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Parse KEGG XML files')
    parser.add_argument('--kegg-dir', required=True, help='Directory containing KEGG XML files')
    parser.add_argument('--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    # Parse KEGG
    kegg_parser = KEGGParser(args.kegg_dir)
    df = kegg_parser.parse_all_files()
    
    # Save
    df.to_csv(args.output, index=False, sep='\t')
    logger.info(f"Saved {len(df)} interactions to {args.output}")


if __name__ == '__main__':
    main()
