rule run_pgd:
    input:
        rules.rds_to_hdf5.output,
        rules.infer_tree.output,
    output:
        "results/pgd_celloracle/cells.h5ad",
    script:
        "../scripts/pgd/run_pgd.py"


rule pgd_extract_metadata:
    input:
        rules.run_pgd.output,
    output:
        "results/pgd_celloracle/cells.metadata.parquet",
    script:
        "../scripts/celloracle/extract_metadata.py"


rule pgd_prepare_oracle:
    input:
        rules.run_pgd.output,
        rules.pruned_network.output,
    output:
        "results/pgd_celloracle/cells.celloracle.oracle",
    script:
        "../scripts/celloracle/prepare_oracle.py"


rule pgd_prepare_gradient:
    input:
        rules.pgd_prepare_oracle.output[0],
    output:
        "results/pgd_celloracle/cells.celloracle.gradient",
    params:
        n_grid=config["params"]["celloracle"]["n_grid"],
        min_mass=config["params"]["celloracle"]["min_mass"],
    script:
        "../scripts/celloracle/prepare_gradient.py"


rule pgd_run_simulation:
    input:
        oracle=rules.pgd_prepare_oracle.output[0],
        gradient=rules.pgd_prepare_gradient.output[0],
        branches=rules.export_branches.output[0],
    output:
        perturb_scores="results/pgd_celloracle/simulations/{factor}/perturb_scores.parquet",
        mc_counts="results/pgd_celloracle/simulations/{factor}/mc_counts.parquet",
        mc_transitions="results/pgd_celloracle/simulations/{factor}/mc_transitions.parquet",
    params:
        n_grid=config["params"]["celloracle"]["n_grid"],
        min_mass=config["params"]["celloracle"]["min_mass"],
    threads: 4
    script:
        "../scripts/celloracle/run_simulation.py"


rule pgd_aggregate_perturb_scores:
    input:
        expand(rules.pgd_run_simulation.output.perturb_scores, factor=factor_list),
    output:
        "results/pgd_celloracle/perturb_scores.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_perturb_scores.py"


rule pgd_aggregate_mc_transitions:
    input:
        metadata=rules.pgd_extract_metadata.output[0],
        mc_transitions=expand(
            rules.pgd_run_simulation.output.mc_transitions,
            factor=factor_list,
        ),
    output:
        "results/pgd_celloracle/mc_transitions/{step}.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_mc_transitions.py"


# ------------------------------


rule pgd_run_simulation_overexpress:
    input:
        oracle=rules.pgd_prepare_oracle.output[0],
        gradient=rules.pgd_prepare_gradient.output[0],
        branches=rules.export_branches.output[0],
    output:
        perturb_scores="results/pgd_celloracle_overexpress/simulations/{factor}/perturb_scores.parquet",
        mc_counts="results/pgd_celloracle_overexpress/simulations/{factor}/mc_counts.parquet",
        mc_transitions="results/pgd_celloracle_overexpress/simulations/{factor}/mc_transitions.parquet",
    params:
        n_grid=config["params"]["celloracle"]["n_grid"],
        min_mass=config["params"]["celloracle"]["min_mass"],
    threads: 4
    script:
        "../scripts/celloracle/run_simulation_overexpress.py"


rule pgd_aggregate_perturb_scores_overexpress:
    input:
        expand(
            rules.pgd_run_simulation_overexpress.output.perturb_scores,
            factor=factor_list,
        ),
    output:
        "results/pgd_celloracle_overexpress/perturb_scores.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_perturb_scores.py"


rule pgd_aggregate_mc_transitions_overexpress:
    input:
        metadata=rules.pgd_extract_metadata.output[0],
        mc_transitions=expand(
            rules.pgd_run_simulation_overexpress.output.mc_transitions,
            factor=factor_list,
        ),
    output:
        "results/pgd_celloracle_overexpress/mc_transitions/{step}.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_mc_transitions.py"
