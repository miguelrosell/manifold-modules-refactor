import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import os

def main():
    print("Generating dummy dataset with pre-computed PCA...")

    # 1. Create synthetic count data
    # We generate random integers to simulate raw gene expression counts
    # Shape: 100 cells x 50 genes
    n_cells = 100
    n_genes = 50
    X = np.random.randint(0, 100, size=(n_cells, n_genes))

    # 2. Create metadata
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    # 3. Initialize AnnData object
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # 4. Preprocessing (Simulating upstream pipeline steps)
    # The 'scdownstream' pipeline typically performs normalization and log transformation
    # before running manifold learning.
    print("Normalizing and log-transforming...")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # 5. Compute PCA (CRITICAL STEP)
    # Your PHATE and Diffmap modules expect 'X_pca' to be present in adata.obsm.
    # We compute 20 components to simulate a real analysis.
    print("Computing PCA (n_comps=20)...")
    sc.tl.pca(adata, n_comps=20)

    # 6. Save the dataset
    output_file = "test_dataset.h5ad"
    print(f"Saving to {output_file}...")
    try:
        adata.write(output_file)
        print("Done! Dataset is ready for testing.")
    except Exception as e:
        print(f"Error saving file: {e}")

if __name__ == "__main__":
    main()
