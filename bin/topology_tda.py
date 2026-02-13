#!/usr/bin/env python3

"""
Script for Topological Data Analysis (TDA) using Giotto-TDA
Part of nf-core/manifold pipeline
"""

import argparse
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from gtda.homology import VietorisRipsPersistence
from gtda.diagrams import PersistenceEntropy
import logging

# Set up logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser(description='Compute Persistent Homology for scRNA-seq data')
    
    parser.add_argument('--input', type=str, required=True, help='Path to input AnnData file (processed with geometry)')
    parser.add_argument('--output', type=str, required=True, help='Path to output AnnData file (.h5ad) with topology info')
    parser.add_argument('--homology_dimensions', type=str, default="0,1", help='Dimensions to compute (0=components, 1=loops, 2=voids)')
    parser.add_argument('--version', action='store_true', help='Print software version and exit')
    
    return parser.parse_args()

def print_version():
    print(f"topology_tda.py: 1.0.0")
    print(f"scanpy: {sc.__version__}")
    print(f"giotto-tda: 0.6.0") # Hardcoded based on our env

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

    # 1. Get the latent representation (PHATE or PCA)
    # If PHATE was computed, use it. Otherwise fall back to PCA.
    if 'X_phate' in adata.obsm:
        logging.info("Using PHATE embeddings for Topological Analysis.")
        X = adata.obsm['X_phate']
    elif 'X_pca' in adata.obsm:
        logging.info("Using PCA embeddings for Topological Analysis.")
        X = adata.obsm['X_pca']
    else:
        logging.info("No embeddings found. Computing PCA...")
        sc.tl.pca(adata)
        X = adata.obsm['X_pca']

    # 2. Compute Persistent Homology (Vietoris-Rips)
    logging.info("Computing Persistent Homology (Vietoris-Rips Filtration)...")
    
    # Parse dimensions (e.g., "0,1" -> [0, 1])
    dims = [int(x) for x in args.homology_dimensions.split(',')]
    
    vr = VietorisRipsPersistence(homology_dimensions=dims, n_jobs=-1)
    diagrams = vr.fit_transform(X[None, :, :]) # Add batch dimension

    # 3. Compute Persistence Entropy (A summary metric of topological complexity)
    pe = PersistenceEntropy()
    entropy = pe.fit_transform(diagrams)
    
    # Store results in AnnData
    # We store the entropy score in .uns (unstructured metadata)
    adata.uns['persistence_entropy'] = entropy[0]
    logging.info(f"Topological Entropy calculated: {entropy[0]}")

    # Note: Storing the full diagram in .uns is complex due to h5py limitations, 
    # so for now we just store the complexity score.

    logging.info(f"Saving results to: {args.output}")
    adata.write_h5ad(args.output)
    logging.info("Topology analysis completed successfully.")

if __name__ == "__main__":
    main()
