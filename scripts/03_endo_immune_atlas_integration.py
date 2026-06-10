"""
03_endo_immune_atlas_integration.py

Integrates HVGs from all tissue types and batch corrects using Harmony.

Input:  data/processed/GSE179640_qc.h5ad
Output: data/processed/GSE179640_integration.h5ad
"""
# -- Load packages
import os
import scanpy as sc
import numpy as np
import harmonypy
import anndata

# -- Set paths
input_path = snakemake.input[0]
output_path = snakemake.output[0]

# -- Load object
combined = sc.read_h5ad(input_path)

# -- Normalize counts
sc.pp.normalize_total(combined)
sc.pp.log1p(combined)

# -- Select HVGs
sc.pp.highly_variable_genes(combined, n_top_genes=snakemake.params.n_top_genes, batch_key="sample_id")
combined = combined[:, combined.var["highly_variable"]].copy()

# -- Dimensionality Reduction
sc.tl.pca(combined)

# -- Batch Correction - Harmony
ho = harmonypy.run_harmony(
    combined.obsm["X_pca"],
    combined.obs,
    "sample_id"
)

combined.obsm["X_pca_harmony"] = ho.Z_corr

# -- Nearest Neighbors
sc.pp.neighbors(
    combined,
    use_rep="X_pca_harmony",
    n_neighbors=snakemake.params.n_neighbors
)

# -- UMAP
sc.tl.umap(combined)

# -- Save integration object for integration
combined.write_h5ad(output_path)