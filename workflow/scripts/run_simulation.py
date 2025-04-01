import pandas as pd
import celloracle as co
from celloracle.applications import Oracle_development_module


def run_simulation(inputs, outputs, wildcards, params):
    oracle = co.load_hdf5(inputs["oracle"])
    gradient = co.load_hdf5(inputs["gradient"])
    goi = wildcards["factor"]
    perturb_condition = {goi: 0.0}

    # In silico TF perturbation
    oracle.simulate_shift(
        perturb_condition=perturb_condition,
        n_propagation=3,
    )

    # Get transition probability
    oracle.estimate_transition_prob(n_neighbors=200)

    # Calculate embedding
    oracle.calculate_embedding_shift()

    # Markov simulation
    oracle.run_markov_chain_simulation(calculate_randomized=False)

    # Calculated simulated vectors
    oracle.calculate_p_mass(n_grid=params["n_grid"])
    oracle.calculate_mass_filter(min_mass=params["min_mass"])

    # Compare simulation vectors with development vectors
    dev = Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=gradient)
    dev.load_perturb_simulation_data(oracle_object=oracle)
    dev.calculate_inner_product()
    dev.calculate_digitized_ip()

    df_perturb_scores = pd.concat(
        [
            dev.inner_product_df,
            pd.DataFrame(
                dev.gridpoints_coordinates[~dev.mass_filter_simulation],
                columns=["x", "y"],
            ),
        ],
        axis=1,
    )
    df_mc_counts = oracle.count_cells_in_mc_resutls(oracle.cluster_column_name)
    df_mc_transitions = oracle.markvov_transition_id

    # Save results
    df_perturb_scores.to_parquet(outputs["perturb_scores"])
    df_mc_counts.to_parquet(outputs["mc_counts"])
    df_mc_transitions.columns = df_mc_transitions.columns.astype(str)
    df_mc_transitions.to_parquet(outputs["mc_transitions"])


snakemake = snakemake  # type: ignore
run_simulation(snakemake.input, snakemake.output, snakemake.wildcards, snakemake.params)
