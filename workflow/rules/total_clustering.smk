import os
BASE = os.environ.get("PIPELINE_BASE", config["paths"]["base"])

rule qc:
    input:
        os.path.join(BASE, config["paths"]["integration_object"])
    output:
        os.path.join(BASE, config["paths"]["total_clustering_object"])
    params:
        resolution   = config["total_clustering"]["resolution"],
        logfoldthreshold   = config["total_clustering"]["logfoldthreshold"],
        padj_threshold   = config["total_clustering"]["padj_threshold"]
    log:
        os.path.join(BASE, config["paths"]["logs"], "total_clustering.log")
    script:
        "../../scripts/04_endo_immune_atlas_total_clustering.py"