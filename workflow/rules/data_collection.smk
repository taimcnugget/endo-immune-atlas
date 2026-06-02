import os
BASE = os.environ.get("PIPELINE_BASE", config["paths"]["base"])
rule data_collection:
    output:
        os.path.join(BASE, config["paths"]["raw_object"])
    log:
        os.path.join(BASE, config["paths"]["logs"], "data_collection.log")
    script:
        "../../scripts/01_endo_immune_atlas_data_collection.py"