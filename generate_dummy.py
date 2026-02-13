import scanpy as sc
import numpy as np
import pandas as pd
import anndata

# 1. Crear matriz de conteo falsa (100 células, 50 genes)
counts = np.random.randint(0, 100, size=(100, 50))

# 2. Crear objeto AnnData
adata = anndata.AnnData(X=counts)

# 3. Añadir metadatos falsos (necesario para que parezca real)
adata.obs_names = [f"cell_{i}" for i in range(100)]
adata.var_names = [f"Gene_{i}" for i in range(50)]
adata.obs['batch'] = 'sample_1'
adata.obs['condition'] = np.random.choice(['control', 'treated'], 100)

# 4. Guardar como .h5ad
print("Generando test_dataset.h5ad...")
adata.write("test_dataset.h5ad")
print("¡Hecho! Listo para la ciencia.")
