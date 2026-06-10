import os
BASE = os.environ.get("PIPELINE_BASE", config["paths"]["base"])

rule integration:
    input:
        os.path.join(BASE, config["paths"]["qc_object"])
    output:
        os.path.join(BASE, config["paths"]["integration_object"])
    params:
        n_top_genes  = config["integration"]["n_top_genes"],
        n_pcs        = config["integration"]["n_pcs"],
        n_neighbors  = config["integration"]["n_neighbors"]

    log:
        os.path.join(BASE, config["paths"]["logs"], "integration.log")
    script:
        "../../scripts/03_endo_immune_atlas_integration.py"