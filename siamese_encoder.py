import os
import copy
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv
from torch_geometric.utils import dropout_edge
import pandas as pd
import numpy as np
import networkx as nx
from scipy.linalg import orthogonal_procrustes

# =============================================================================
# 1. GRAPH ENCODER ARCHITECTURE
# =============================================================================
class GNNEncoder(nn.Module):
    """
    Base Graph Neural Network Encoder (fθ).
    Maps structural neighborhood features into a dense latent space.
    """
    def __init__(self, in_channels, hidden_channels, out_channels):
        super(GNNEncoder, self).__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)
        self.activation = nn.PReLU()

    def forward(self, x, edge_index, edge_weight=None):
        x = self.conv1(x, edge_index, edge_weight)
        x = self.activation(x)
        x = self.conv2(x, edge_index, edge_weight)
        return x

# =============================================================================
# 2. SELF-SUPERVISED SIAMESE MODEL (BGRL)
# =============================================================================
class BGRLSiamese(nn.Module):
    """
    Bootstrapped Graph Latent Representation (BGRL) Model.
    Uses an 'Online' network that learns to predict the output of a 'Target' network.
    The Target network is an Exponential Moving Average (EMA) of the Online network.
    """
    def __init__(self, encoder, hidden_dim, pred_dim):
        super(BGRLSiamese, self).__init__()
        # Online network (Encoder + Predictor)
        self.online_encoder = encoder
        self.predictor = nn.Sequential(
            nn.Linear(hidden_dim, pred_dim),
            nn.PReLU(),
            nn.Linear(pred_dim, hidden_dim)
        )
        
        # Target network (EMA of Encoder)
        self.target_encoder = copy.deepcopy(self.online_encoder)
        
        # Freeze target network (no gradient tracking)
        for param in self.target_encoder.parameters():
            param.requires_grad = False

    def update_target_network(self, tau=0.99):
        """EMA update for the target network."""
        for online_param, target_param in zip(self.online_encoder.parameters(), self.target_encoder.parameters()):
            target_param.data = tau * target_param.data + (1.0 - tau) * online_param.data

    def forward_online(self, x, edge_index, edge_weight=None):
        return self.predictor(self.online_encoder(x, edge_index, edge_weight))

    def forward_target(self, x, edge_index, edge_weight=None):
        with torch.no_grad():
            return self.target_encoder(x, edge_index, edge_weight)
            
    def loss(self, p, z):
        """Negative Cosine Similarity Loss"""
        p = F.normalize(p, p=2, dim=1)
        z = F.normalize(z, p=2, dim=1)
        return 2 - 2 * (p * z).sum(dim=-1).mean()

# =============================================================================
# 3. ALIGNMENT & DISPLACEMENT UTILITIES (AIM 2)
# =============================================================================
def align_and_score(emb_normal, emb_tumor, gene_names, anchor_genes=None):
    """
    Solves the cross-condition latent alignment problem using Orthogonal Procrustes.
    Returns the structural displacement (delta) for each gene.
    """
    # If no anchors provided, use all genes (less robust but works for testing)
    if anchor_genes is None:
        anchor_indices = list(range(len(gene_names)))
    else:
        anchor_indices = [i for i, g in enumerate(gene_names) if g in anchor_genes]
    
    A_anchor = emb_normal[anchor_indices]
    B_anchor = emb_tumor[anchor_indices]
    
    # Fit Orthogonal Map (R) that aligns Tumor to Normal based on stable anchors
    R, scale = orthogonal_procrustes(B_anchor, A_anchor)
    
    # Align the entire Tumor latent space
    emb_tumor_aligned = np.dot(emb_tumor, R)
    
    # Calculate Displacement (Euclidean distance in aligned latent space)
    displacement = np.linalg.norm(emb_normal - emb_tumor_aligned, axis=1)
    
    results = pd.DataFrame({
        'gene': gene_names,
        'displacement_score': displacement
    }).sort_values('displacement_score', ascending=False)
    
    return results, emb_normal, emb_tumor_aligned

# =============================================================================
# 4. DATA LOADER UTILITY (Bridges DiCE Pipeline to PyTorch)
# =============================================================================
def load_dice_network_to_pyg(file_path):
    """Converts a DiCE weighted network file into a PyTorch Geometric Data object."""
    df = pd.read_csv(file_path, sep='\t')
    
    # Create unified node vocabulary
    nodes = sorted(list(set(df['source']).union(set(df['target']))))
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    
    # Map edges to indices
    src_idx = df['source'].map(node_to_idx).values
    tgt_idx = df['target'].map(node_to_idx).values
    edge_index = torch.tensor(np.vstack((src_idx, tgt_idx)), dtype=torch.long)
    
    # Edge weights (invert DiCE distance back to similarity weight if necessary)
    if 'weight' in df.columns:
        edge_weight = torch.tensor(df['weight'].values, dtype=torch.float)
    else:
        edge_weight = torch.ones(edge_index.shape[1], dtype=torch.float)
        
    # Dummy structural features (e.g., node degree) since Aim 2 uses topology
    # In a full model, this could be the HKS (Heat Kernel Signature) from Aim 1
    degrees = torch.bincount(edge_index[0], minlength=len(nodes)).float().view(-1, 1)
    # Normalize features
    x = (degrees - degrees.mean()) / (degrees.std() + 1e-8)
    
    data = Data(x=x, edge_index=edge_index, edge_attr=edge_weight)
    return data, nodes

# =============================================================================
# 5. TRAINING LOOP
# =============================================================================
def train_bgrl(model, data, optimizer, epochs=100, p_drop=0.2):
    """Trains the Siamese encoder on a single graph using self-supervision."""
    model.train()
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        
        # Data Augmentation: View 1 (Drop edges)
        edge_index_1, mask_1 = dropout_edge(data.edge_index, p=p_drop)
        edge_weight_1 = data.edge_attr[mask_1] if data.edge_attr is not None else None
        
        # Data Augmentation: View 2 (Drop edges differently)
        edge_index_2, mask_2 = dropout_edge(data.edge_index, p=p_drop)
        edge_weight_2 = data.edge_attr[mask_2] if data.edge_attr is not None else None

        # Online predictions
        p1 = model.forward_online(data.x, edge_index_1, edge_weight_1)
        p2 = model.forward_online(data.x, edge_index_2, edge_weight_2)
        
        # Target representations
        z1 = model.forward_target(data.x, edge_index_1, edge_weight_1)
        z2 = model.forward_target(data.x, edge_index_2, edge_weight_2)
        
        # Symmetric Loss
        loss = model.loss(p1, z2) + model.loss(p2, z1)
        loss.backward()
        optimizer.step()
        
        # Update Target Network (EMA)
        model.update_target_network()
        
        if (epoch + 1) % 20 == 0:
            print(f"Epoch {epoch+1:03d} | BGRL Loss: {loss.item():.4f}")

# =============================================================================
# 6. EXECUTION SCRIPT
# =============================================================================
if __name__ == "__main__":
    print("=== DiCE-FLOW: Siamese GNN Encoder Initialization ===")
    
    # 1. File Paths (Adjust to point to your actual DiCE pipeline output)
    normal_net_path = "data/processed/weighted/network_normal.txt"
    tumor_net_path = "data/processed/weighted/network_tumor.txt"
    
    # For testing, we create dummy files if real ones don't exist
    if not os.path.exists(normal_net_path):
        print("Mocking data for test...")
        os.makedirs("data/processed/weighted", exist_ok=True)
        pd.DataFrame({'source': ['A','B','C','D'], 'target': ['B','C','D','A'], 'weight': [0.1, 0.2, 0.9, 0.5]}).to_csv(normal_net_path, sep='\t', index=False)
        pd.DataFrame({'source': ['A','B','C','E'], 'target': ['B','E','D','A'], 'weight': [0.9, 0.1, 0.2, 0.8]}).to_csv(tumor_net_path, sep='\t', index=False)

    # 2. Load Data
    data_normal, genes_n = load_dice_network_to_pyg(normal_net_path)
    data_tumor, genes_t = load_dice_network_to_pyg(tumor_net_path)
    
    # Ensure graph node alignment (simplification: assume they share the same gene universe)
    all_genes = sorted(list(set(genes_n).union(set(genes_t))))
    
    print(f"Loaded Normal Graph: {data_normal.num_nodes} nodes, {data_normal.num_edges} edges")
    print(f"Loaded Tumor Graph: {data_tumor.num_nodes} nodes, {data_tumor.num_edges} edges")

    # 3. Initialize Model
    # Feature dim = 1 (just degrees for now), Latent dim = 64
    encoder = GNNEncoder(in_channels=1, hidden_channels=128, out_channels=64)
    siamese_model = BGRLSiamese(encoder, hidden_dim=64, pred_dim=32)
    optimizer = torch.optim.Adam(siamese_model.parameters(), lr=0.005, weight_decay=1e-5)

    # 4. Train the Siamese Model (Shared weights over both conditions)
    print("\n--- Training on Normal Graph ---")
    train_bgrl(siamese_model, data_normal, optimizer, epochs=40)
    
    print("\n--- Training on Tumor Graph ---")
    train_bgrl(siamese_model, data_tumor, optimizer, epochs=40)

    # 5. Extract Embeddings (Using the EMA Target Encoder for stability)
    siamese_model.eval()
    with torch.no_grad():
        emb_norm = siamese_model.forward_target(data_normal.x, data_normal.edge_index).numpy()
        emb_tum = siamese_model.forward_target(data_tumor.x, data_tumor.edge_index).numpy()

    # 6. Perform Orthogonal Procrustes Alignment (Aim 2)
    print("\n--- Performing Cross-Condition Alignment ---")
    # In a real scenario, use genes known to NOT be related to cancer as anchors
    # Here, we pass None to use all genes, or define a list: anchors = ["GAPDH", "ACTB"]
    displacement_df, aligned_norm, aligned_tum = align_and_score(emb_norm, emb_tum, genes_n, anchor_genes=None)
    
    print("\nTop 5 Genes by Structural Latent Displacement (Rewiring):")
    print(displacement_df.head(5))
    
    # Save the output
    os.makedirs("data/results/dice_flow", exist_ok=True)
    displacement_df.to_csv("data/results/dice_flow/aim2_latent_displacement.csv", index=False)
    print("\nResults saved to data/results/dice_flow/aim2_latent_displacement.csv")