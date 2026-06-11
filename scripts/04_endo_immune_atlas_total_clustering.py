"""
04_endo_immune_atlas_total_clustering.py

Uses Leiden clustering and manual cell labeling
to annotate cell clusters. Output are labeled 
clusters in an object. 

Input:  data/processed/GSE179640_integration.h5ad
Output: data/processed/GSE179640_total_clustering.h5ad
"""
# -- Load packages
import os
import pandas as pd 
import scanpy as sc
import numpy as np
import anndata as ad
from matplotlib import pyplot as plt

# -- Set paths
input_path = snakemake.input[0]
output_path = snakemake.output[0]

# -- Load data
combined = sc.read_h5ad(input_path)

# -- Cluster

sc.tl.leiden(combined, 
			 resolution = snakemake.params.resolution,
			 flavor="igraph",
                 )

# DEG
sc.tl.rank_genes_groups(
    combined,
    groupby="leiden_res_0.01",
    method="wilcoxon"
)

markers = sc.get.rank_genes_groups_df(
    combined,
    group=None
)

markers = markers[
    (markers["logfoldchanges"] > snakemake.params.logfoldthreshold) &
    (markers["pvals_adj"] < snakemake.params.padj_threshold)]

# -- Marker Gene List
marker_genes = {
    "Myeloid":     ["CD68", "CD14", "CD86", "FCGR1A", "LYZ",
                    "S100A8", "S100A9", "ITGAM", "CSF1R", "MPEG1"],

    "Lymphoid":    ["CD3D", "CD3E", "CD3G", "CD8A", "CD4",
                    "GZMA", "GZMB", "NKG7", "NCAM1", "KLRD1"],

    "Stromal":     ["COL1A1", "COL1A2", "PDGFRA", "DCN", "LUM",
                    "VIM", "ACTA2", "THY1", "POSTN"],

    "Epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19",
                    "WFDC2", "MUC1", "CLDN3", "CLDN4", "FXYD3"],

    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "EMCN",
                    "ENG", "CLEC14A", "MCAM", "FLT1", "KDR"]
}

# -- Label clusters
combined.obs["cluster_label"] = combined.obs["leiden_res_0.01"].map(
    {
        "0": "Myeloid",
        "1": "Lymphoid",
        "2": "Stromal",
        "3": "Epithelial",
        "4": "Endothelial"
    }
)

# -- Save clustering object for subsetting
combined.write_h5ad(output_path)