#!/usr/bin/env python3

"""
Script for Geometric Manifold Learning (PHATE / Diffusion Maps)
Part of nf-core/manifold pipeline
"""

import argparse
import sys
import scanpy as sc
import phate
import pandas as pd
import logging

# Set up logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser(description='Compute Manifold Embeddings for scRNA-seq data')
    
    parser.add_argument('--input', type=str, required=True, help='Path to input AnnData file (.h5ad)')
    parser.add_argument('--output', type=str, required=True, help='Path to output AnnData file (.h5ad)')
    parser.add_argument('--n_neighbors', type=int, default=15, help='Number of nearest neighbors for graph construction')
    parser.add_argument('--method', type=str, default='phate', choices=['phate', 'diffmap'], help='Manifold learning method to use')
    parser.add_argument('--version', action='store_true', help='Print software version and exit')
    
    return parser.parse_args()

def print_version():
    """
    Prints versions of key libraries for Nextflow provenance tracking.
    """
    print(f"geometry_manifold.py: 1.0.0")
    print(f"scanpy: {sc.__version__}")
    print(f"phate: {phate.__version__}")
    print(f"pandas: {pd.__version__}")

def main():
    args = parse_args()

    # Dynamic versioning check (CRITICAL for nf-core compliance)
    if args.version:
        print_version()
        sys.exit(0)

    logging.info(f"Reading input file: {args.input}")
    try:
        adata = sc.read_h5ad(args.input)
    except Exception as e:
        logging.error(f"Failed to read input file: {e}")
        sys.exit(1)

    logging.info(f"Preprocessing: Computing neighbors (k={args.n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors)

    if args.method == 'phate':
        logging.info("Computing PHATE embedding (Potential of Heat-diffusion for Affinity-based Transition Embedding)...")
        
        # Initialize PHATE operator
        phate_op = phate.PHATE(n_jobs=-1, knn=args.n_neighbors)
        # Fit and transform
        data_phate = phate_op.fit_transform(adata.to_df())
        
        # Store in AnnData object (following standard convention)
        adata.obsm['X_phate'] = data_phate
        logging.info("PHATE coordinates stored in adata.obsm['X_phate']")

    elif args.method == 'diffmap':
        logging.info("Computing Diffusion Maps...")
        sc.tl.diffmap(adata)
        logging.info("Diffusion Map coordinates stored in adata.obsm['X_diffmap']")

    logging.info(f"Saving results to: {args.output}")
    adata.write_h5ad(args.output)
    logging.info("Process completed successfully.")

if __name__ == "__main__":
    main()
