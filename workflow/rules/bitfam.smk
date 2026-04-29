rule download_targets:
    output:
        "downloads/{factor}.tsv",
    params:
        url="https://chip-atlas.dbcls.jp/data/mm10/target/{factor}.5.tsv",
    shell:
        """
        wget {params.url} -O {output}
        """


rule whitelist_targets:
    input:
        rules.download_targets.output,
    output:
        "results/bitfam/network/whitelist_targets/{factor}.csv",
    params:
        whitelist_regex=config["params"]["bitfam"]["network"]["whitelist_regex"],
    script:
        "../scripts/bitfam/whitelist_targets.py"


rule chip_atlas:
    input:
        expand(rules.whitelist_targets.output, factor=atlas_factors),
    output:
        "results/bitfam/network/chip_atlas.csv",
    run:
        import pandas as pd

        df = pd.concat([pd.read_csv(f) for f in input])
        df.to_csv(output[0], index=False)


rule zbtb46:
    input:
        config["resources"]["zbtb46"],
    output:
        "results/bitfam/network/zbtb46.csv",
    script:
        "../scripts/bitfam/zbtb46.py"


rule full_network:
    input:
        chip_atlas=rules.chip_atlas.output,
        zbtb46=rules.zbtb46.output,
    output:
        "results/bitfam/network/full.csv",
    run:
        import pandas as pd

        df = pd.concat([pd.read_csv(f) for f in input])
        df.to_csv(output[0], index=False)


rule pruned_network:
    input:
        cells=config["inputs"]["cells"],
        network=rules.full_network.output,
    output:
        "results/bitfam/network/pruned.csv",
    script:
        "../scripts/bitfam/pruned_network.R"


rule whitelist_experiments:
    input:
        rules.pruned_network.output,
    output:
        "results/bitfam/network/whitelist_experiments.tsv",
    params:
        cache_dir="results/bitfam/network/whitelist_experiments_cache",
    script:
        "../scripts/bitfam/whitelist_experiments.py"


rule run_bitfam:
    input:
        cells=config["inputs"]["cells"],
        network=rules.pruned_network.output,
    output:
        "results/bitfam/results/res.rds",
    params:
        min_targets=10,
    script:
        "../scripts/bitfam/run_bitfam.R"


rule extract_acts:
    input:
        cells=config["inputs"]["cells"],
        res=rules.run_bitfam.output,
    output:
        "results/bitfam/results/acts.rds",
    script:
        "../scripts/bitfam/extract_acts.R"


# For multiple runs with different random seeds ----
rule run_bitfam_runs:
    input:
        cells=config["inputs"]["cells"],
        network=rules.pruned_network.output,
    output:
        "results/bitfam_runs/run_{run}/res.rds",
    params:
        min_targets=10,
    script:
        "../scripts/bitfam_runs/run_bitfam.R"


rule extract_acts_runs:
    input:
        cells=config["inputs"]["cells"],
        res=rules.run_bitfam_runs.output,
    output:
        "results/bitfam_runs/run_{run}/acts.rds",
    script:
        "../scripts/bitfam/extract_acts.R"


rule aggregate_acts:
    input:
        run0=rules.extract_acts.output,
        run1=rules.extract_acts_runs.output[0].format(run=1),
        run2=rules.extract_acts_runs.output[0].format(run=2),
        run3=rules.extract_acts_runs.output[0].format(run=3),
    output:
        "results/bitfam_runs/acts.csv",
    script:
        "../scripts/bitfam_runs/aggregate_acts.R"
