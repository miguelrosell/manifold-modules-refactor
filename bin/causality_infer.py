#!/usr/bin/env python3

"""
Script for Causal Network Inference
Part of nf-core/manifold pipeline
"""

import argparse
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import logging

# Set up logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser(description='Infer Causal Gene Regulatory Networks from Manifold')
    
    parser.add_argument('--input', type=str, required=True, help='Path to input AnnData file (with geometry and topology)')
    parser.add_argument('--output', type=str, required=True, help='Path to output AnnData file (.h5ad)')
    parser.add_argument('--n_top_genes', type=int, default=100, help='Number of top variable genes to use for network inference')
    parser.add_argument('--version', action='store_true', help='Print software version and exit')
    
    return parser.parse_args()

def print_version():
    print(f"causality_infer.py: 1.0.0")
    print(f"scanpy: {sc.__version__}")
    print(f"networkx: {nx.__version__}")
    print(f"numpy: {np.__version__}")

def main():
    args = parse_args()

    if args.version:
        print_version()
        sys.exit(0)

    logging.info(f"Reading input file: {args.input}")
    try:
        adata = sc.read_h5ad(args.input)
    except Exception as e:
        logging.error(f"Failed to read input file: {e}")
        sys.exit(1)

    # --- FIX 1: DATA TYPE CONVERSION ---
    # Check if data is floating point. If not (it's integer), convert to float32.
    if not np.issubdtype(adata.X.dtype, np.floating):
        logging.info("Converting expression matrix to float32...")
        adata.X = adata.X.astype(np.float32)

    # --- FIX 2: NORMALIZATION & LOG TRANSFORM ---
    # Critical: highly_variable_genes expects log-normalized data.
    # Without log1p, algorithms like 'seurat' or 'cell_ranger' will fail mathematically.
    logging.info("Normalizing and log-transforming data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 1. Select top variable genes (nodes of our network)
    logging.info(f"Selecting top {args.n_top_genes} genes for causal inference...")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes)
    genes_use = adata.var_names[adata.var['highly_variable']]
    
    # Subset matrix to these genes
    X_subset = adata[:, genes_use].X
    
    # Handle sparse matrices if necessary
    if hasattr(X_subset, "toarray"):
        X_subset = X_subset.toarray()

    # 2. Build a Correlation/Causal Graph
    logging.info("Inferring network structure...")
    corr_matrix = np.corrcoef(X_subset.T)
    
    # FIX 3: Handle NaN values (in case a gene has 0 variance)
    corr_matrix = np.nan_to_num(corr_matrix)
    
    # Create Graph
    G = nx.Graph()
    genes_list = list(genes_use)
    
    threshold = 0.5 # Correlation threshold
    edges_added = 0
    
    for i in range(len(genes_list)):
        for j in range(i + 1, len(genes_list)):
            weight = corr_matrix[i, j]
            if abs(weight) > threshold:
                G.add_edge(genes_list[i], genes_list[j], weight=weight)
                edges_added += 1

    logging.info(f"Inferred network with {len(G.nodes)} nodes and {edges_added} edges.")

    # 3. Compute Centrality (Who are the 'Master Regulators'?)
    centrality = nx.degree_centrality(G)
    
    # Store centrality scores in adata.var
    adata.var['causal_centrality'] = 0.0
    for gene, score in centrality.items():
        if gene in adata.var_names:
            adata.var.loc[gene, 'causal_centrality'] = score

    logging.info("Centrality scores stored in adata.var['causal_centrality']")

    logging.info(f"Saving results to: {args.output}")
    adata.write_h5ad(args.output)
    logging.info("Causality analysis completed successfully.")

if __name__ == "__main__":
    main()
