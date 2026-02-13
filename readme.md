# nf-core/manifold 🧬

**A Nextflow pipeline for Geometric Deep Learning, Topological Data Analysis (TDA), and Causal Inference in scRNA-seq data.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-%E2%89%A520.10.0-0db7ed.svg)](https://www.docker.com/)

## 🚀 Introduction

**nf-core/manifold** is a specialized bioinformatics pipeline designed to go beyond standard clustering. It reconstructs the continuous trajectory of cells (Manifold Learning), analyzes the complexity of that shape (Topology), and infers the gene regulatory networks driving the process (Causality).

The pipeline performs three main steps:
1.  **Geometry:** Dimensionality reduction using PHATE to preserve global structure.
2.  **Topology:** Persistent Homology calculation (Vietoris-Rips filtration) to detect cycles and voids.
3.  **Causality:** Inference of Gene Regulatory Networks (GRN) based on manifold-based correlation and centrality.

`![Pipeline Results](dashboard_results.png)`

## 🛠️ Usage

To run the pipeline with Docker (actually recommended):

```bash
nextflow run main.nf \
    --input your_dataset.h5ad \
    --outdir results \
    -profile docker

```

### 📂 Input Format

* **Format:** `.h5ad` (AnnData object from Scanpy/AnnData).
* **Requirements:** A count matrix (raw or preprocessed).

### 🖥️ Output Structure

The results will be stored in the `results/` directory:

```bash
results/
├── geometry/          # Contains PHATE embeddings
│   └── *_geometric.h5ad
├── topology/          # Contains Topological features (Persistence Entropy)
│   └── *_topology.h5ad
└── causality/         # Contains Inferred Gene Networks & Centrality Scores
    └── *_causal.h5ad

```

## 📊 Visualization 

You can visualize the final results using the provided Python script:

```bash
# Using the pipeline's Docker image
docker run --rm \
    -u $(id -u):$(id -g) \
    -v "$(pwd):/app" \
    -w /app \
    -e MPLCONFIGDIR=/tmp \
    -e NUMBA_CACHE_DIR=/tmp \
    nf-core/manifold:dev \
    python plot_dashboard.py

```

## 📜 Citations

This pipeline uses the following key libraries:

* [Nextflow](https://www.nextflow.io/)
* [Scanpy](https://scanpy.readthedocs.io/)
* [PHATE](https://phate.readthedocs.io/)
* [Giotto-TDA](https://giotto-ai.github.io/gtda-docs/)
* [NetworkX](https://networkx.org/)

---

Created by Miguel Rosell Hidalgo as part of a bioinformatics portfolio.

```

