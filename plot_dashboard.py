#!/usr/bin/env python3

"""
Visualization Dashboard for nf-core/manifold pipeline.
Generates a summary plot containing:
1. Geometric Manifold Projection (PHATE)
2. Causal Gene Regulatory Network (Top Drivers)
"""

import scanpy as sc
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import os
import sys

def main():
    # Define input path (output from the pipeline)
    input_file = "results/causality/test_dataset_causal.h5ad"
    output_img = "dashboard_results.png"

    print(f"Loading data from: {input_file}...")

    # Validate file existence
    if not os.path.exists(input_file):
        print(f"Error: File not found at {input_file}")
        print("Please ensure the pipeline completed successfully and 'results' directory exists.")
        sys.exit(1)

    try:
        adata = sc.read_h5ad(input_file)
    except Exception as e:
        print(f"Error reading .h5ad file: {e}")
        sys.exit(1)

    # Initialize Dashboard (2 Subplots)
    fig, ax = plt.subplots(1, 2, figsize=(18, 8))
    plt.subplots_adjust(wspace=0.25)

    # --- PANEL 1: GEOMETRIC PROJECTION (PHATE) ---
    print("Plotting Geometric Manifold (PHATE)...")
    
    if 'X_phate' in adata.obsm:
        phate_coords = adata.obsm['X_phate']
        
        # Color by condition if available, else use single color
        if 'condition' in adata.obs:
            groups = adata.obs['condition'].unique()
            for group in groups:
                mask = adata.obs['condition'] == group
                ax[0].scatter(
                    phate_coords[mask, 0], 
                    phate_coords[mask, 1], 
                    label=group, 
                    s=60, 
                    alpha=0.7, 
                    edgecolors='w', 
                    linewidth=0.5
                )
            ax[0].legend(title="Condition", loc='best')
        else:
            ax[0].scatter(
                phate_coords[:, 0], 
                phate_coords[:, 1], 
                c='teal', 
                s=60, 
                alpha=0.7, 
                edgecolors='w', 
                linewidth=0.5
            )
            
        ax[0].set_title("Manifold Geometry (PHATE)", fontsize=16, fontweight='bold')
        ax[0].set_xlabel("PHATE 1", fontsize=12)
        ax[0].set_ylabel("PHATE 2", fontsize=12)
        ax[0].grid(True, linestyle='--', alpha=0.3)
        
    else:
        ax[0].text(0.5, 0.5, "PHATE coordinates not found", ha='center', fontsize=14)
        ax[0].axis('off')

    # --- PANEL 2: CAUSAL REGULATORY NETWORK ---
    print("Plotting Causal Regulatory Network...")

    # Check if centrality scores exist (computed in step 3 of pipeline)
    if 'causal_centrality' in adata.var:
        # Select top 15 master regulators based on centrality
        top_genes = adata.var.sort_values('causal_centrality', ascending=False).head(15).index.tolist()
        
        # Subset data for visualization
        adata_top = adata[:, top_genes]
        
        # Compute correlation matrix for edge visualization
        # (Note: In a full production script, we would visualize the causal graph edges directly)
        if hasattr(adata_top.X, "toarray"):
            matrix = adata_top.X.toarray()
        else:
            matrix = adata_top.X
            
        corr_matrix = np.corrcoef(matrix.T)
        
        # Initialize NetworkX graph
        G = nx.Graph()
        
        # Add nodes with size proportional to centrality
        for i, gene in enumerate(top_genes):
            score = adata.var.loc[gene, 'causal_centrality']
            # Scale node size for visibility
            G.add_node(gene, size=score * 3000)

        # Add edges based on correlation threshold
        threshold = 0.3
        rows, cols = np.where(np.abs(corr_matrix) > threshold)
        
        for r, c in zip(rows, cols):
            if r < c: # Avoid duplicates and self-loops
                weight = corr_matrix[r, c]
                G.add_edge(top_genes[r], top_genes[c], weight=weight)

        # Compute layout
        pos = nx.spring_layout(G, k=0.5, seed=42)
        
        # Draw Nodes
        node_sizes = [nx.get_node_attributes(G, 'size')[n] for n in G.nodes()]
        nx.draw_networkx_nodes(
            G, pos, 
            ax=ax[1], 
            node_size=node_sizes, 
            node_color='tomato', 
            alpha=0.9, 
            edgecolors='black'
        )
        
        # Draw Labels
        nx.draw_networkx_labels(
            G, pos, 
            ax=ax[1], 
            font_size=11, 
            font_weight='bold', 
            font_color='black'
        )
        
        # Draw Edges
        nx.draw_networkx_edges(
            G, pos, 
            ax=ax[1], 
            width=1.5, 
            alpha=0.4, 
            edge_color='gray'
        )
        
        ax[1].set_title("Inferred Regulatory Network (Top Drivers)", fontsize=16, fontweight='bold')
        ax[1].axis('off')

    else:
        ax[1].text(0.5, 0.5, "Causality metadata not found", ha='center', fontsize=14)
        ax[1].axis('off')

    # --- SAVE OUTPUT ---
    plt.tight_layout()
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Success! Dashboard saved to: {os.path.abspath(output_img)}")

if __name__ == "__main__":
    main()
