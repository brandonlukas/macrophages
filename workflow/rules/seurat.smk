rule findallmarkers:
    input:
        rules.prepare_cells.output,
    output:
        "results/seurat/rna/findallmarkers.csv",
    params:
        assay="RNA",
    script:
        "../scripts/seurat/findallmarkers.R"


rule findallmarkers_z:
    input:
        rules.prepare_cells.output,
    output:
        "results/seurat/z/findallmarkers.csv",
    params:
        assay="Z",
    script:
        "../scripts/seurat/findallmarkers.R"


rule find_de:
    input:
        rules.prepare_cells.output,
    output:
        "results/seurat/rna/find_de.csv",
    params:
        assay="RNA",
    script:
        "../scripts/seurat/find_de.R"


rule find_de_z:
    input:
        rules.prepare_cells.output,
    output:
        "results/seurat/z/find_de.csv",
    params:
        assay="Z",
    script:
        "../scripts/seurat/find_de.R"
