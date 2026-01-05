rule rds_to_hdf5:
    input:
        rules.prepare_cells.output,
        rules.pruned_network.output,
    output:
        "results/celloracle/cells.h5ad",
    params:
        conda_env=config["params"]["celloracle"]["conda_env"],
    script:
        "../scripts/celloracle/rds_to_hdf5.R"


rule extract_metadata:
    input:
        rules.rds_to_hdf5.output,
    output:
        "results/celloracle/cells.metadata.parquet",
    script:
        "../scripts/celloracle/extract_metadata.py"


rule prepare_oracle:
    input:
        rules.rds_to_hdf5.output,
        rules.pruned_network.output,
    output:
        "results/celloracle/cells.celloracle.oracle",
    script:
        "../scripts/celloracle/prepare_oracle.py"


rule prepare_gradient:
    input:
        rules.prepare_oracle.output[0],
    output:
        "results/celloracle/cells.celloracle.gradient",
    params:
        n_grid=config["params"]["celloracle"]["n_grid"],
        min_mass=config["params"]["celloracle"]["min_mass"],
    script:
        "../scripts/celloracle/prepare_gradient.py"


rule run_simulation:
    input:
        oracle=rules.prepare_oracle.output[0],
        gradient=rules.prepare_gradient.output[0],
        branches=rules.export_branches.output[0],
    output:
        perturb_scores="results/celloracle/simulations/{factor}/perturb_scores.parquet",
        mc_counts="results/celloracle/simulations/{factor}/mc_counts.parquet",
        mc_transitions="results/celloracle/simulations/{factor}/mc_transitions.parquet",
    params:
        n_grid=config["params"]["celloracle"]["n_grid"],
        min_mass=config["params"]["celloracle"]["min_mass"],
    threads: 4
    script:
        "../scripts/celloracle/run_simulation.py"


rule aggregate_perturb_scores:
    input:
        expand(rules.run_simulation.output.perturb_scores, factor=factor_list),
    output:
        "results/celloracle/perturb_scores.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_perturb_scores.py"


rule aggregate_mc_transitions:
    input:
        metadata=rules.extract_metadata.output[0],
        mc_transitions=expand(
            rules.run_simulation.output.mc_transitions,
            factor=factor_list,
        ),
    output:
        "results/celloracle/mc_transitions/{step}.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_mc_transitions.py"


# ------------------------------


rule run_simulation_overexpress:
    input:
        oracle=rules.prepare_oracle.output[0],
        gradient=rules.prepare_gradient.output[0],
        branches=rules.export_branches.output[0],
    output:
        perturb_scores="results/celloracle_overexpress/simulations/{factor}/perturb_scores.parquet",
        mc_counts="results/celloracle_overexpress/simulations/{factor}/mc_counts.parquet",
        mc_transitions="results/celloracle_overexpress/simulations/{factor}/mc_transitions.parquet",
    params:
        n_grid=config["params"]["celloracle"]["n_grid"],
        min_mass=config["params"]["celloracle"]["min_mass"],
    threads: 4
    script:
        "../scripts/celloracle/run_simulation_overexpress.py"


rule aggregate_perturb_scores_overexpress:
    input:
        expand(
            rules.run_simulation_overexpress.output.perturb_scores, factor=factor_list
        ),
    output:
        "results/celloracle_overexpress/perturb_scores.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_perturb_scores.py"


rule aggregate_mc_transitions_overexpress:
    input:
        metadata=rules.extract_metadata.output[0],
        mc_transitions=expand(
            rules.run_simulation_overexpress.output.mc_transitions,
            factor=factor_list,
        ),
    output:
        "results/celloracle_overexpress/mc_transitions/{step}.parquet",
    params:
        factor_list=factor_list,
    script:
        "../scripts/celloracle/aggregate_mc_transitions.py"
