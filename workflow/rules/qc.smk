import os
BASE = os.environ.get("PIPELINE_BASE", config["paths"]["base"])

rule qc:
    input:
        os.path.join(BASE, config["paths"]["raw_object"])
    output:
        os.path.join(BASE, config["paths"]["qc_object"])
    params:
        min_genes    = config["qc"]["min_genes"],
        min_counts   = config["qc"]["min_UMI"],
        max_counts   = config["qc"]["max_UMI"],
        mt_threshold = config["qc"]["mt_threshold"],
        min_cells    = config["qc"]["min_cells"]
    log:
        os.path.join(BASE, config["paths"]["logs"], "qc.log")
    script:
        "../../scripts/02_endo_immune_atlas_qc.py"