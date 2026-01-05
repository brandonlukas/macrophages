rule infer_tree:
    input:
        config["inputs"]["cells"],
    output:
        "results/lamian/infer_tree.rds",
    params:
        seed_1=42,
        seed_2=12345,
    script:
        "../scripts/lamian/infer_tree.R"


rule evaluate_uncertainty:
    input:
        rules.infer_tree.output,
    output:
        "results/lamian/evaluate_uncertainty.rds",
    script:
        "../scripts/lamian/evaluate_uncertainty.R"


rule prepare_cells:
    input:
        cells=rules.extract_acts.output,
        res=rules.infer_tree.output,
        aggregate_acts=rules.aggregate_acts.output,
    output:
        "results/seurat/cells.rds",
    script:
        "../scripts/seurat/prepare_cells.R"


rule perform_tde:
    input:
        cells=rules.prepare_cells.output,
        res=rules.infer_tree.output,
    output:
        "results/lamian/rna/tde.rds",
    params:
        assay="RNA",
        test_type="time",
        nfeatures=5000,
    threads: 24
    script:
        "../scripts/lamian/perform_test.R"


rule perform_xde:
    input:
        cells=rules.prepare_cells.output,
        res=rules.infer_tree.output,
    output:
        "results/lamian/rna/xde.rds",
    params:
        assay="RNA",
        test_type="variable",
        nfeatures=5000,
    threads: 24
    script:
        "../scripts/lamian/perform_test.R"


rule perform_tde_z:
    input:
        cells=rules.prepare_cells.output,
        res=rules.infer_tree.output,
    output:
        "results/lamian/z/tde.rds",
    params:
        assay="Z_avg",
        test_type="time",
    threads: 8
    script:
        "../scripts/lamian/perform_test.R"


rule perform_xde_z:
    input:
        cells=rules.prepare_cells.output,
        res=rules.infer_tree.output,
    output:
        "results/lamian/z/xde.rds",
    params:
        assay="Z_avg",
        test_type="variable",
    threads: 8
    script:
        "../scripts/lamian/perform_test.R"


rule export_branches:
    input:
        rules.infer_tree.output,
    output:
        "results/lamian/branches.csv",
    script:
        "../scripts/lamian/export_branches.R"
