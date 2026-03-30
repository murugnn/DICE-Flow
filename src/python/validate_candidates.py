import pandas as pd
import argparse
import logging

logging.basicConfig(level=logging.INFO)

def load_network_degree(path):
    """Calculate degree (number of connections) for every gene in a network file"""
    try:
        # Read Source, Target columns
        df = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1], names=['Source', 'Target'])
        
        # Stack source and target to get a list of all nodes involved in edges
        all_nodes = pd.concat([df['Source'], df['Target']])
        
        # Count frequency of each node = Degree
        degree_counts = all_nodes.value_counts()
        return degree_counts
    except Exception as e:
        logging.error(f"Error loading {path}: {e}")
        return pd.Series()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--top-genes', required=True, help='Path to dice_genes_top500.txt')
    parser.add_argument('--normal-net', required=True, help='Path to processed network_normal.txt')
    parser.add_argument('--tumor-net', required=True, help='Path to processed network_tumor.txt')
    args = parser.parse_args()

    # 1. Load Top Candidates
    candidates = pd.read_csv(args.top_genes, sep='\t')
    top_20 = candidates.head(20)['gene'].tolist()
    
    print(f"Validating Top 20 Candidates: {top_20}")

    # 2. Calculate Raw Degrees directly from network files
    print("Calculating network degrees (this might take a moment)...")
    deg_normal = load_network_degree(args.normal_net)
    deg_tumor = load_network_degree(args.tumor_net)

    # 3. Build Validation Table
    print("\n" + "="*65)
    print(f"{'GENE':<10} | {'Normal Deg':<12} | {'Tumor Deg':<12} | {'Change':<10} | {'Status'}")
    print("="*65)

    for gene in top_20:
        d_norm = deg_normal.get(gene, 0)
        d_tumor = deg_tumor.get(gene, 0)
        
        diff = d_tumor - d_norm
        
        # Status Check
        if d_tumor == 0 and d_norm == 0:
            status = "ORPHAN (Error)" # Should not happen if centrality > 0
        elif d_tumor < 5 and d_norm < 5:
            status = "Low Confidence" # Too few edges to be a reliable hub
        elif diff > 0:
            status = "Gained Hubness"
        elif diff < 0:
            status = "Lost Hubness"
        else:
            status = "Stable"

        print(f"{gene:<10} | {d_norm:<12} | {d_tumor:<12} | {diff:<+10} | {status}")
    print("="*65)
    print("NOTE: If 'Tumor Deg' is 0 but it's a top candidate, you have a math artifact (division by zero).")

if __name__ == "__main__":
    main()