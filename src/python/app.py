import streamlit as st
import pandas as pd
import networkx as nx
from pyvis.network import Network
import tempfile
import numpy as np

# --- CONFIGURATION ---
st.set_page_config(layout="wide", page_title="DiCE-Duan Explorer")

# --- DATA LOADING --- 
@st.cache_data
def load_graph(file_path):
    """
    Loads a weighted network safely.
    Fixes applied:
    1. skips header row
    2. only reads first 3 cols (ignores weight_int)
    3. forces string types for genes
    """
    try:
        # 1. Read CSV: Explicitly use columns 0, 1, 2 (Source, Target, Weight)
        df = pd.read_csv(file_path, sep='\t', 
                         header=0,   # Skip the first row (header labels)
                         usecols=[0, 1, 2], # Only read Source, Target, Weight
                         names=['Source', 'Target', 'Weight'], # Standardize names
                         dtype={'Source': str, 'Target': str, 'Weight': float},
                         na_values=['nan', 'NaN'])
        
        # 2. Drop incomplete rows
        df = df.dropna(subset=['Source', 'Target'])
        
        # 3. Filter out "Numeric" gene names (Garbage Data Fix)
        # Sometimes huge files have glitches where a number appears as a gene name
        def is_valid_gene(name):
            try:
                # If it can be cast to a float (e.g., "0.412"), it's likely noise
                float(name) 
                return False
            except ValueError:
                # If it fails to cast (e.g., "PARK2"), it's a real gene name
                return True

        df = df[df['Source'].apply(is_valid_gene) & df['Target'].apply(is_valid_gene)]
        
        # 4. Create Graph
        G = nx.from_pandas_edgelist(df, 'Source', 'Target', edge_attr='Weight', create_using=nx.Graph())
        return G
    except Exception as e:
        st.error(f"Error loading {file_path}: {e}")
        return None

# --- VISUALIZATION FUNCTION ---
def visualize_subgraph(G, nodes_to_highlight, title, highlight_color="#FF4B4B"):
    """Creates an interactive PyVis graph"""
    nt = Network(height="450px", width="100%", bgcolor="#222222", font_color="white")
    
    subgraph = G.subgraph(nodes_to_highlight)
    
    for node in subgraph.nodes():
        node_id = str(node) # Force string to prevent Pyvis crash
        color = highlight_color if node in nodes_to_highlight else "#97c2fc"
        size = 25 if node in nodes_to_highlight else 15
        nt.add_node(node_id, label=node_id, color=color, size=size)
    
    for u, v, data in subgraph.edges(data=True):
        weight = data.get('Weight', 0)
        # Visual trick: Thicker line = Stronger Traffic
        # Assuming Weight is correlation (0-1) or Cost (0-1). 
        # Adjust logic if needed. Here we assume generic "Importance".
        try:
            thickness = 1 + (float(weight) * 3)
        except:
            thickness = 1
        
        nt.add_edge(str(u), str(v), value=thickness, color="#555555", title=f"Weight: {weight:.4f}")

    nt.from_nx(subgraph)
    nt.toggle_physics(True)
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp:
        nt.save_graph(tmp.name)
        return tmp.name

# --- MAIN UI ---
st.title("🧬 DiCE-Duan: Cancer Network Rewiring Explorer")

# 1. Load Data
st.sidebar.header("Data Sources")
# Ensure these paths match your server location
normal_path = "data/processed/weighted/network_normal.txt"
tumor_path = "data/processed/weighted/network_tumor.txt"

status_text = st.sidebar.empty()
status_text.text("Loading networks (this may take 30s)...")

G_norm = load_graph(normal_path)
G_tumor = load_graph(tumor_path)

if G_norm is None or G_tumor is None:
    st.stop()

status_text.success(f"Online: Normal ({len(G_norm)} nodes) | Tumor ({len(G_tumor)} nodes)")

# --- TABS ---
tab1, tab2 = st.tabs(["🚀 Feature 1: The Path Switch", "🕸️ Feature 2: Rewiring Proof"])

# === FEATURE 1: SHORTEST PATH SWITCH ===
with tab1:
    st.markdown("### Trace Signaling Routes")
    col1, col2 = st.columns(2)
    with col1: source_gene = st.text_input("Source Gene", value="EGFR").upper().strip()
    with col2: target_gene = st.text_input("Target Gene", value="MYC").upper().strip()

    if st.button("Trace Signal"):
        # Check existence
        if source_gene not in G_norm:
            st.error(f"Source '{source_gene}' not found in Normal network.")
        elif target_gene not in G_norm:
            st.error(f"Target '{target_gene}' not found in Normal network.")
        else:
            try:
                # Calculate Paths
                path_norm = nx.shortest_path(G_norm, source=source_gene, target=target_gene, weight='Weight')
                cost_norm = nx.shortest_path_length(G_norm, source=source_gene, target=target_gene, weight='Weight')
                
                # Check tumor existence
                if source_gene in G_tumor and target_gene in G_tumor:
                    path_tumor = nx.shortest_path(G_tumor, source=source_gene, target=target_gene, weight='Weight')
                    cost_tumor = nx.shortest_path_length(G_tumor, source=source_gene, target=target_gene, weight='Weight')
                    delta_val = cost_norm - cost_tumor
                else:
                    path_tumor = []
                    cost_tumor = 0
                    delta_val = 0
                    st.warning("Path broken in Tumor network (Gene missing).")

                # METRICS
                m1, m2 = st.columns(2)
                m1.metric("Normal Path Cost", f"{cost_norm:.4f}")
                if path_tumor:
                    m2.metric("Tumor Path Cost", f"{cost_tumor:.4f}", delta=f"{delta_val:.4f}")

                # VISUALS
                v1, v2 = st.columns(2)
                with v1:
                    st.info(f"Normal: {' → '.join(path_norm)}")
                    h_n = visualize_subgraph(G_norm, path_norm, "Normal", "#00CCFF")
                    st.components.v1.html(open(h_n).read(), height=450)
                
                with v2:
                    if path_tumor:
                        st.warning(f"Tumor: {' → '.join(path_tumor)}")
                        h_t = visualize_subgraph(G_tumor, path_tumor, "Tumor", "#FF4B4B")
                        st.components.v1.html(open(h_t).read(), height=450)
            except nx.NetworkXNoPath:
                st.error("No path exists between these genes.")

# === FEATURE 2: NEIGHBORHOOD REWIRING ===
with tab2:
    st.markdown("### Prove Functional Rewiring (Traffic Change)")
    st.caption("Validates if a gene has gained stronger/faster connections in cancer.")
    
    hub_gene = st.text_input("Inspect Candidate Gene", value="PARK2").upper().strip()
    
    if st.button("Analyze Hub"):
        if hub_gene not in G_norm:
            st.error(f"Gene '{hub_gene}' not found in network.")
        else:
            # Get neighbors from Normal graph
            nbrs_raw = list(G_norm.neighbors(hub_gene))
            
            # --- SAFE LOOP (Comparing Intersection) ---
            weights_norm = []
            weights_tumor = []
            valid_neighbors = []
            
            for n in nbrs_raw:
                # Only include neighbors that exist in BOTH graphs to calculate stats
                if n in G_tumor and G_tumor.has_edge(hub_gene, n):
                    w_n = G_norm[hub_gene][n]['Weight']
                    w_t = G_tumor[hub_gene][n]['Weight']
                    
                    weights_norm.append(w_n)
                    weights_tumor.append(w_t)
                    valid_neighbors.append(n)
            
            if not weights_norm:
                st.warning("No valid overlapping neighbors found (Check data consistency).")
            else:
                avg_n = np.mean(weights_norm)
                avg_t = np.mean(weights_tumor)
                
                # Display Stats
                st.write(f"**{hub_gene}** Connectivity Profile:")
                c1, c2, c3 = st.columns(3)
                c1.metric("Valid Neighbors", len(valid_neighbors), help="Neighbors present in both networks")
                
                # Logic: If weight is Correlation, Higher is Better. If Cost, Lower is Better.
                # Just showing raw values here.
                c2.metric("Avg Weight (Normal)", f"{avg_n:.4f}")
                c3.metric("Avg Weight (Tumor)", f"{avg_t:.4f}", delta=f"{avg_t - avg_n:.4f}")

                st.write("---")
                
                # Visuals
                v1, v2 = st.columns(2)
                # Show max 20 neighbors to keep UI fast
                nodes_show = [hub_gene] + valid_neighbors[:20] 
                
                with v1:
                    st.caption("Normal Traffic (First 20)")
                    h_n = visualize_subgraph(G_norm, nodes_show, "Norm", "#00CCFF")
                    st.components.v1.html(open(h_n).read(), height=450)
                with v2:
                    st.caption("Tumor Traffic (First 20)")
                    h_t = visualize_subgraph(G_tumor, nodes_show, "Tumor", "#FF4B4B")
                    st.components.v1.html(open(h_t).read(), height=450)