"""
02_endo_immune_atlas_qc.py

Runs QC on GSE179640 processed object and returns
a cleaned object - removes doublets!

Input:  data/processed/GSE179640_raw.h5ad
Output: data/processed/GSE179640_qc.h5ad
"""

# -- Load libraries
import os
import scanpy as sc
import pandas as pd
import numpy as np

# -- Set paths
input_path = snakemake.input[0]
output_path = snakemake.output[0]

# -- Load data
combined = sc.read_h5ad(input_path)

# --- Quality Control
# hemoglobin genes (hemos)
combined.var["hemos"] = combined.var_names.str.contains("^HB[^(P)]")

# ribosomal genes (ribos)
combined.var["ribos"] = combined.var_names.str.startswith(("RPS", "RPL"))

# mitochondrial genes (mitos)
combined.var["mt"] = combined.var_names.str.contains("MT-")


# -- Calculate QC
sc.pp.calculate_qc_metrics(combined,
                           qc_vars = ["mt"],
                           percent_top = None,
                           inplace = True,
                           log1p= False)

# -- Filter low quality cells
sc.pp.filter_cells(combined, min_genes = snakemake.params.min_genes)
sc.pp.filter_cells(combined, min_UMI = snakemake.params.min_UMI)
sc.pp.filter_cells(combined, max_UMI = snakemake.params.max_UMI)


# -- Remove cells with high mt %
combined.obs["outlier_mt"] = combined.obs["pct_counts_mt"] > snakemake.params.mt_threshold
combined = combined[~combined.obs["outlier_mt"]].copy()

# -- Doublet removal
sc.pp.scrublet(combined, batch_key="sample_id")
combined = combined[~combined.obs["predicted_doublet"]].copy()


# -- Filter genes
sc.pp.filter_genes(combined, min_cells = snakemake.params.min_cells)

# -- Save object
combined.write_h5ad(output_path)