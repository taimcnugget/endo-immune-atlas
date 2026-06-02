"""
01_endo_immune_atlas_data_collection.py

Loads GSE179640 h5 files, labels samples with metadata,
concatenates into a single AnnData object and writes to disk.

Input:  GEO supplementary files (downloaded manually)
Output: data/processed/GSE179640_raw.h5ad
"""


# -- Load packages
import anndata as ad
import scanpy as sc
import os
import numpy as np
import pandas as pd

# -- Config from snakemake
output_path = snakemake.output[0]

# -- Set up metadata
# -- Samples to exclude
EXCLUDE = {"EOR", "EcPA"}

# -- Metadata map keyed by sample_id extracted from filename
metadata = {
    # controls
    "GSM6102532_C01_Ctrl": {"patient_id": "C01", "tissue_type": "healthy",          "condition": "control",         "lesion_site": None},
    "GSM6102533_C02_Ctrl": {"patient_id": "C02", "tissue_type": "healthy",          "condition": "control",         "lesion_site": None},
    "GSM6102534_C03_Ctrl": {"patient_id": "C03", "tissue_type": "healthy",          "condition": "control",         "lesion_site": None},

    # eutopic
    "GSM6102537_E01_EuE":  {"patient_id": "E01", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102540_E02_EuE":  {"patient_id": "E02", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102543_E03_EuE":  {"patient_id": "E03", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102546_E04_EuE":  {"patient_id": "E04", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102549_E05_EuE":  {"patient_id": "E05", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102551_E06_EuE":  {"patient_id": "E06", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102554_E07_EuE":  {"patient_id": "E07", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102555_E08_EuE":  {"patient_id": "E08", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},
    "GSM6102560_E09_EuE":  {"patient_id": "E09", "tissue_type": "eutopic",          "condition": "endometriosis",   "lesion_site": None},

    # ectopic peritoneal
    "GSM6595248_E01_EcP":  {"patient_id": "E01", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},
    "GSM6595250_E02_EcP":  {"patient_id": "E02", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},
    "GSM6595252_E03_EcP":  {"patient_id": "E03", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},
    "GSM6595254_E04_EcP":  {"patient_id": "E04", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},
    "GSM6595256_E05_EcP":  {"patient_id": "E05", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},
    "GSM6102550_E06_EcP":  {"patient_id": "E06", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},
    "GSM6102553_E07_EcP":  {"patient_id": "E07", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},
    "GSM6595258_E09_EcP":  {"patient_id": "E09", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "peritoneal"},

    # ectopic ovarian
    "GSM6102552_E07_EcO":  {"patient_id": "E07", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "ovarian"},
    "GSM6102556_E09_EcO":  {"patient_id": "E09", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "ovarian"},
    "GSM6595261_E10_EcO":  {"patient_id": "E10", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "ovarian"},
    "GSM6102562_E11_EcO":  {"patient_id": "E11", "tissue_type": "ectopic",          "condition": "endometriosis",   "lesion_site": "ovarian"},
}

# -- Collect h5 files
h5_files = []

for root, dirs, files in os.walk(snakemake.input[0]):
    for f in files:
        if f.endswith(".h5"):
            h5_files.append(os.path.join(root, f))
            
# -- Set up object
adatas = []

for file in h5_files:
    sample_id = os.path.basename(file).replace("_filtered_feature_bc_matrix.h5", "")

# -- Exclude organoids and unwanted tissue types (Adjacent Peritoneium)

    if any(excl in sample_id for excl in EXCLUDE):
      print(f"Skipping: {sample_id}")
      continue

# -- Exclude samples without any metadata

    if sample_id not in metadata:
      print(
          f"Skipping {sample_id} - No metadata found."
          )
      continue

# -- Load data
    adata = sc.read_10x_h5(file)
    adata.var_names_make_unique()

# -- Add metadata to the object
    meta = metadata[sample_id]
    adata.obs["sample_id"] = sample_id
    adata.obs["patient_id"] = meta['patient_id']
    adata.obs["tissue_type"] = meta["tissue_type"]
    adata.obs["condition"] = meta["condition"]
    adata.obs["lesion_site"] = meta["lesion_site"]
    adata.obs["dataset"] = "GSE179640"

# -- Append sample name to cell barcodes
    adata.obs_names = [
        f"{sample_id}_{bc}"
        for bc in adata.obs_names
    ]

    print(f"{sample_id} loaded - {adata.n_obs} total cells")
    adatas.append(adata)

# -- Concat
combined = ad.concat(adatas, join = "outer")
combined.var_names_make_unique()

# -- Write output
combined.write_h5ad(output_path)