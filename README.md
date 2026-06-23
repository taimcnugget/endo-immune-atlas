# endo-immune-atlas

A reproducible Snakemake pipeline for spatial immune and immunosenescence profiling in endometriosis using public single-cell and spatial transcriptomics data.

---

## Background

Endometriosis affects approximately 1 in 10 women of reproductive age and is characterized by the growth of endometrial-like tissue outside the uterus. Despite its prevalence, the immune landscape of endometriotic lesions remains poorly characterized at single-cell resolution. This project integrates scRNA-seq and spatial transcriptomics to map immune cell states, with a focused analysis on immunosenescence  the acquisition of a senescent phenotype by immune cells  within ectopic, eutopic, and healthy endometrial tissue.

---

## Biological Questions

1. How do immune cell states differ across ectopic, eutopic, and healthy endometrium?
2. Are immunosenescent populations enriched in endometriotic lesions?
3. How do senescence-high immune cells spatially organize within the lesion microenvironment?
4. What ligand-receptor interactions are enriched between senescent immune cells and other cell types in the lesion niche?

---

## Datasets

| Accession | Type | Description |
|---|---|---|
| GSE179640 | scRNA-seq | Ectopic, eutopic, and healthy endometrium (Tan et al.) |
| GSM6690475 / GSM6690476 | Visium ST | Ectopic endometriosis lesion tissue sections |

---

## Pipeline Overview

This project is built as a modular Snakemake pipeline. Each analysis step is an independent rule with defined inputs, outputs, and parameters controlled via `config/config.yaml`.

```
endo-immune-atlas/
 config/
    config.yaml
 workflow/
    Snakefile
    rules/
        qc.smk
        integration.smk
        total_clustering.smk
        subset_immune_cells.smk
        immunosenescence.smk
        spatial_transcriptomics.smk
        neighborhood_analysis.smk
        niche_discovery.smk
        cell_cell_communication.smk
        network_analysis.smk
 scripts/
    01_endo_immune_atlas_data_collection.py
    02_endo_immune_atlas_qc.py
    03_endo_immune_atlas_integration.py
    04_endo_immune_atlas_total_clustering.py
    05_endo_immune_atlas_subset_clustering.py
    06_endo_immune_atlas_immunosenescence.py
    07_endo_immune_atlas_spatial_transcriptomics.py
    08_endo_immune_atlas_neighborhood_analysis.py
    09_endo_immune_atlas_niche_discovery.py
    10_endo_immune_atlas_cellular_communication.py
    11_endo_immune_atlas_network_analysis.py
 notebooks/          # exploratory analysis only
 data/
    raw/            # GEO downloads  not tracked by git
    processed/      # pipeline intermediates  not tracked by git
 results/
    preprocessing/  # QC plots, UMAPs, cluster markers
    figures/        # final analysis figures
 envs/
    endo_env.yaml
 .gitignore
 README.md
```

### DAG

```mermaid
flowchart TD
    gse[(scRNA-seq<br/>GSE179640)]
    gsm[(Visium ST<br/>GSM6690475/76)]

    subgraph preprocess["Single-Cell Atlas Construction"]
        qc[qc]
        integration[integration]
        cluster[clustering_annotation]
        subset[subset_immune_cells]
    end

    subgraph cellstate["Cell-State Analysis"]
        immuno[immunosenescence]
    end

    subgraph spatial["Spatial Analysis"]
        spatial_map[spatial_transcriptomics]
        neighborhood[neighborhood_analysis]
        niche[niche_discovery]
    end

    subgraph systems["Systems Biology"]
        ccc[cell_cell_communication]
        network[network_analysis]
    end

    subgraph outputs["Outputs"]
        preprocess_figures([Preprocessing Reports])
        final_figures([Final Figures & Reports])
    end

    gse --> qc --> integration --> cluster --> subset --> immuno

    gsm --> spatial_map

    immuno --> spatial_map
    spatial_map --> neighborhood
    neighborhood --> niche
    niche --> ccc
    ccc --> network

    qc --> preprocess_figures
    integration --> preprocess_figures
    cluster --> preprocess_figures
    subset --> preprocess_figures

    immuno --> final_figures
    spatial_map --> final_figures
    neighborhood --> final_figures
    niche --> final_figures
    ccc --> final_figures
    network --> final_figures

    classDef dataset fill:#D3D3D3,stroke:#000,color:#000;
    classDef scrna fill:#BFD7ED,stroke:#000,color:#000;
    classDef state fill:#A5DFCD,stroke:#000,color:#000;
    classDef spatial fill:#F7B267,stroke:#000,color:#000;
    classDef systems fill:#7D82B8,stroke:#000,color:#000;
    classDef output fill:#EF767A,stroke:#000,color:#000,font-weight:bold;

    class gse,gsm dataset;
    class qc,integration,cluster,subset scrna;
    class immuno state;
    class spatial_map,neighborhood,niche spatial;
    class ccc,network systems;
    class preprocess_figures,final_figures output;
```

---

## Notebook Outline

| Notebook | Tools | Description |
|---|---|---|
| 01_data_collection | AnnData | Data loading, AnnData generation, sample labeling |
| 02_qc | Scanpy | Quality control, cell filtering, doublet removal |
| 03_integration | Scanpy, harmonypy | Dataset integration, batch correction, PCA |
| 04_total_clustering | Scanpy | Broad cell type annotation |
| 05_subset_clustering | Scanpy, CellTypist | Immune subsetting, fine-resolution clustering/annotation |
| 06_immunosenescence | Scanpy | Senescence and dysfunction scoring, cell state classification, DEG analysis|
| 07_spatial_transcriptomics | cell2location | Spatial deconvolution, spatial immune mapping within lesions |
| 08_neighborhood analysis | Squidpy | Identification of spatially enriched cellular neighborhoods and cell-state co-localization patterns |
| 09_niche_discovery | CellCharter | Discovery/characterization of multicellular tissue niches |
| 10_cellular_communication | LIANA+ | Ligand-receptor inference between neighboring cell populations and tissue niches |
| 11_network_analysis | NetworkX, igraph, pandas | Construction and analysis of cell-state interaction networks to identify key signaling hubs and microenvironmental programs |

---

### Setup

```bash
# clone the repo
git clone https://github.com/taimcnugget/endo-immune-atlas.git
cd endo-immune-atlas

# create the conda environment
conda env create -f envs/endo_env.yaml
conda activate endo_pipeline

# download data from GEO
# GSE179640: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640
# GSM6690475: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6690475
# GSM6690476: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6690476
# place downloaded files in data/raw/ following the structure in config/config.yaml

# dry run to verify pipeline
snakemake -n

# run the full pipeline
snakemake --cores 4
```

---

## Status

Current phase: Spatial microenvironment analysis

Completed:
- Single-cell immune atlas construction
- Immune cell annotation and lineage characterization
- Immunosenescence and dysfunction scoring
- Cell-state co-occurrence analysis
- Differential expression analysis of senescent immune populations

In progress:
- Spatial deconvolution of endometriosis lesions
- Identification of spatially organized immune microenvironments

Next steps:
- Cellular neighborhood analysis
- Niche discovery
- Cell-cell communication inference
- Network characterization of senescent and dysfunctional tissue ecosystems

---

## Future Directions

*To be completed after analysis. This section will document follow-up experiments and biological questions arising from the results.*

---

## Author

Tailynn Y. McCarty, PhD
Computational Immunology | Women's Health  
[LinkedIn](www.linkedin.com/in/tailynn)  [GitHub](github.com/taimcnugget)
