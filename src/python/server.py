import pandas as pd
import networkx as nx
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from typing import List, Dict, Any, Optional
import logging

# --- CONFIGURATION ---
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DiCE-API")

app = FastAPI(title="DiCE-Duan API", version="2.1")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], 
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- GLOBAL STATE ---
DATA_STORE = {
    "G_normal": None,
    "G_tumor": None,
    "knockout_df": None,
    "top_genes": [],
    "log2fc_map": {},
    "essential_genes": set()
}

BASE_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = BASE_DIR / "data/processed/weighted"
RESULTS_DIR = BASE_DIR / "data/results"
TABLES_DIR = RESULTS_DIR / "validation/tables"

def load_graph(file_path: Path) -> nx.Graph:
    if not file_path.exists():
        logger.error(f"Graph file not found: {file_path}")
        return nx.Graph()
    df = pd.read_csv(file_path, sep='\t')
    # Load with weight (1 - correlation)
    G = nx.from_pandas_edgelist(df, 'source', 'target', edge_attr='weight')
    return G

@app.on_event("startup")
async def startup_event():
    logger.info("Loading DiCE-Duan Data...")
    
    # 1. Load Graphs
    DATA_STORE["G_normal"] = load_graph(PROCESSED_DIR / "network_normal.txt")
    DATA_STORE["G_tumor"] = load_graph(PROCESSED_DIR / "network_tumor.txt")
    
    # 2. Load Knockout Results
    ko_path = RESULTS_DIR / "knockout_analysis.txt"
    if ko_path.exists():
        DATA_STORE["knockout_df"] = pd.read_csv(ko_path, sep='\t')
        
    # 3. Load Top Genes 
    dice_path = RESULTS_DIR / "dice_genes_top500.txt"
    if dice_path.exists():
        df_top = pd.read_csv(dice_path, sep='\t')
        DATA_STORE["top_genes"] = df_top['gene'].tolist()

    # 4. Load Biological Annotations (DEA)
    dea_path = RESULTS_DIR / "dea_results_table.txt"
    if dea_path.exists():
        dea_df = pd.read_csv(dea_path, sep='\t')
        DATA_STORE["log2fc_map"] = dict(zip(dea_df['gene'], dea_df['log2FC']))
        
    # 5. Load Essentiality Annotations (DepMap)
    essential_path = TABLES_DIR / "09_fitness_overlap_genes.csv"
    if essential_path.exists():
        ess_df = pd.read_csv(essential_path)
        DATA_STORE["essential_genes"] = set(ess_df['gene'].tolist())

    logger.info("Server Ready.")

# --- HELPER FUNCTIONS ---

def build_subgraph_json(G: nx.Graph, nodes: list):
    """Creates a D3/Recharts compatible JSON for a subset of nodes"""
    subgraph = G.subgraph(nodes)
    
    log2fc_map = DATA_STORE["log2fc_map"]
    essential_genes = DATA_STORE["essential_genes"]
    
    response_nodes = []
    for n in subgraph.nodes():
        # Add a "group" or "rank" to help coloring
        # We assume the input list 'nodes' is already sorted by rank
        rank = nodes.index(n) + 1 if n in nodes else 999
        log2fc = float(log2fc_map.get(n, 0.0))
        is_essential = n in essential_genes
        
        response_nodes.append({
            "id": n, 
            "label": n, 
            "rank": rank,
            "log2fc": log2fc,
            "is_essential": is_essential
        })
        
    response_links = []
    for u, v, d in subgraph.edges(data=True):
        response_links.append({
            "source": u,
            "target": v,
            "weight": round(d.get('weight', 1.0), 4)
        })
        
    return {"nodes": response_nodes, "links": response_links}

def get_path_details(G, path):
    log2fc_map = DATA_STORE["log2fc_map"]
    essential_genes = DATA_STORE["essential_genes"]
    
    nodes = [{
        "id": n, 
        "label": n,
        "log2fc": float(log2fc_map.get(n, 0.0)),
        "is_essential": n in essential_genes
    } for n in path]
    links = []
    total_cost = 0.0
    
    for i in range(len(path) - 1):
        u, v = path[i], path[i+1]
        if G.has_edge(u, v):
            weight = G[u][v].get('weight', 1.0)
            total_cost += weight
            links.append({
                "source": u,
                "target": v,
                "weight": round(weight, 4)
            })
        
    return {
        "found": True,
        "total_cost": round(total_cost, 4),
        "nodes": nodes,
        "links": links
    }

# --- ENDPOINTS ---

@app.get("/api/network/initial")
async def get_initial_network(limit: int = 150):
    """
    Returns the 'Skeleton Network': The Top N DiCE genes (default 150).
    Input: ?limit=150
    """
    top_genes = DATA_STORE["top_genes"]
    if not top_genes:
        return {"nodes": [], "links": []}
    
    # 1. Slice the list to the requested limit (e.g., top 150)
    subset_genes = top_genes[:limit]
    
    # 2. Use Tumor graph 
    G = DATA_STORE["G_tumor"]
    
    # 3. Build Subgraph
    payload = build_subgraph_json(G, subset_genes)
    
    return payload

@app.get("/api/knockout")
async def get_knockout():
    """Returns top virtual knockout candidates"""
    df = DATA_STORE["knockout_df"]
    if df is None or df.empty:
        return []
    
    df['vitality_score'] = pd.to_numeric(df['vitality_score'], errors='coerce')
    sorted_df = df.sort_values(by="vitality_score", ascending=False).head(50)
    
    records = sorted_df.to_dict(orient="records")
    essential_genes = DATA_STORE["essential_genes"]
    
    for row in records:
        row["is_essential"] = row.get("gene") in essential_genes
        
    return records

@app.get("/api/path")
async def get_path(source: str, target: str):
    """
    Returns full graph data (nodes/edges) for paths in BOTH networks.
    """
    src = source.upper().strip()
    tgt = target.upper().strip()
    
    G_norm = DATA_STORE["G_normal"]
    G_tumor = DATA_STORE["G_tumor"]
    
    res = {
        "normal": {"found": False, "nodes": [], "links": [], "total_cost": 0},
        "tumor": {"found": False, "nodes": [], "links": [], "total_cost": 0},
        "stats": {"cost_diff": 0, "status": "Unknown"}
    }
    
    # Calc Normal Path
    try:
        if src in G_norm and tgt in G_norm:
            p = nx.shortest_path(G_norm, src, tgt, weight="weight")
            res["normal"] = get_path_details(G_norm, p)
    except nx.NetworkXNoPath:
        pass

    # Calc Tumor Path
    try:
        if src in G_tumor and tgt in G_tumor:
            p = nx.shortest_path(G_tumor, src, tgt, weight="weight")
            res["tumor"] = get_path_details(G_tumor, p)
    except nx.NetworkXNoPath:
        pass

    # Calculate Stats
    if res["normal"]["found"] and res["tumor"]["found"]:
        diff = res["normal"]["total_cost"] - res["tumor"]["total_cost"]
        res["stats"]["cost_diff"] = round(diff, 4)
        if diff > 0.1:
            res["stats"]["status"] = "Shortcut Created (Gain of Function)"
        elif diff < -0.1:
            res["stats"]["status"] = "Path Lengthened (Loss of Efficiency)"
        else:
            res["stats"]["status"] = "Stable"
            
    return res

@app.get("/api/rewiring")
async def get_rewiring(gene: str):
    """
    Returns the 'Ego Network' (neighbors) for a gene.
    Classifies edges as 'maintained', 'lost', or 'gained'.
    """
    gene = gene.upper().strip()
    G_n = DATA_STORE["G_normal"]
    G_t = DATA_STORE["G_tumor"]
    
    if gene not in G_t and gene not in G_n:
        return {"nodes": [], "links": []}
        
    neighbors_n = set(G_n.neighbors(gene)) if gene in G_n else set()
    neighbors_t = set(G_t.neighbors(gene)) if gene in G_t else set()
    
    all_neighbors = neighbors_n.union(neighbors_t)
    
    log2fc_map = DATA_STORE["log2fc_map"]
    essential_genes = DATA_STORE["essential_genes"]
    
    nodes = [{
        "id": gene, 
        "type": "hub", 
        "log2fc": float(log2fc_map.get(gene, 0.0)),
        "is_essential": gene in essential_genes
    }]
    links = []
    
    # Limit to top 30 neighbors
    sorted_neighbors = list(all_neighbors)
    sorted_neighbors = sorted_neighbors[:30] 
    
    for n in sorted_neighbors:
        status = "maintained"
        weight = 0
        
        if n in neighbors_n and n not in neighbors_t:
            status = "lost"
            weight = G_n[gene][n]['weight']
        elif n in neighbors_t and n not in neighbors_n:
            status = "gained"
            weight = G_t[gene][n]['weight']
        else:
            status = "maintained"
            weight = (G_n[gene][n]['weight'] + G_t[gene][n]['weight']) / 2
            
        nodes.append({
            "id": n, 
            "type": "neighbor", 
            "status": status,
            "log2fc": float(log2fc_map.get(n, 0.0)),
            "is_essential": n in essential_genes
        })
        links.append({
            "source": gene,
            "target": n,
            "status": status,
            "weight": round(weight, 4)
        })
        
    return {"nodes": nodes, "links": links}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)