# Manifold Learning & TDA Modules for nf-core/scdownstream

This repository contains the refactored modular components for Manifold Learning, Topological Data Analysis (TDA), and Causal Inference. 

The codebase has been restructured from a standalone pipeline into **Nextflow DSL2 modules** and a **Subworkflow** to facilitate integration into `nf-core/scdownstream`.

## 📂 Repository Structure

The logic is divided into three main computational blocks:

1.  **Geometry (Manifold Learning):**
    * `PHATE`: Potential of Heat-diffusion for Affinity-based Transition Embedding.
    * `DIFFMAP`: Diffusion Maps (via Scanpy).
2.  **Topology (TDA):**
    * `TOPOLOGY`: Persistent Homology using **Ripser**. Computes H0 and H1 persistence diagrams on the manifold embeddings.
3.  **Causality (Inference):**
    * `CAUSALITY`: Trajectory inference and graph abstraction using **PAGA** and **Leiden** clustering.

### Directory Layout

.
├── bin/                    # Standalone Python CLI scripts
│   ├── run_phate.py
│   ├── run_diffmap.py
│   ├── run_topology.py     # TDA implementation
│   └── run_causality.py    # PAGA/Leiden implementation
├── modules/local/          # Nextflow DSL2 process definitions
├── subworkflows/local/     # Orchestrating subworkflow
└── tests/                  # nf-test suite


## 🛠️ Dependencies

The Docker environment includes:
* **Core:** `scanpy`, `pandas`, `numpy`, `scikit-learn`
* **Geometry:** `phate`
* **Topology:** `ripser`, `persim`
* **Graph/Causality:** `networkx`, `python-igraph`, `leidenalg`

## 🧪 Testing (nf-test)

This repository includes a full unit test suite using **nf-test**.

### Prerequisites
* Nextflow
* Docker
* nf-test (Download: `curl -fsSL https://code.askimed.com/install/nf-test | bash`)

### Running the Tests
To verify all modules (Phate, Diffmap, Topology, Causality), run:

```bash
./nf-test test --profile docker

Current Status: ✅ All 4 modules passing.
Usage (Subworkflow)

The logic is orchestrated in subworkflows/local/manifold_learning.nf.

Input:

    Channel containing [ meta, h5ad ] (Pre-computed PCA required).

    Methods list (e.g., 'phate,diffmap').

Logic:

    Computes requested Geometries (PHATE/Diffmap).

    Computes Topology (TDA) on top of the geometry.

    Infers Causality/Trajectory from the topological structure.

Refactored by Miguel Rosell for scdownstream integration.
