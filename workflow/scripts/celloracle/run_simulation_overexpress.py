import pandas as pd
import numpy as np
import celloracle as co
from celloracle.applications import Oracle_development_module
import scanpy as sc


# new
def parse_branch_clusters(branch_name: str) -> list[int]:
    """
    Extract the ordered cluster IDs from a branch string, e.g.
    'backbone 12,6,9,3,4,1,11' -> [12, 6, 9, 3, 4, 1, 11]
    'branch: 12,6,9,3,10'      -> [12, 6, 9, 3, 10]
    """
    last_token = branch_name.split()[-1]  # e.g. '12,6,9,3,4,1,11'
    return [int(x) for x in last_token.split(",")]


def dev_module(oracle, gradient, cell_idx=None):
    dev = Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=gradient)
    dev.load_perturb_simulation_data(oracle_object=oracle, cell_idx_use=cell_idx)
    dev.calculate_inner_product()
    dev.calculate_digitized_ip()
    df = pd.concat(
        [
            dev.inner_product_df,
            pd.DataFrame(
                dev.gridpoints_coordinates[~dev.mass_filter_simulation],
                columns=["x", "y"],
            ),
            pd.DataFrame(dev.flow[~dev.mass_filter_simulation], columns=["dx", "dy"]),
            pd.DataFrame(
                dev.ref_flow[~dev.mass_filter_simulation], columns=["ref_dx", "ref_dy"]
            ),
        ],
        axis=1,
    )
    return df


def fetch_perturb_condition(oracle, goi, safe_range_fold=2.0):
    # based on _is_perturb_condition_valid in celloracle/trajectory/oracle_utility.py
    actual_values = sc.get.obs_df(
        oracle.adata, keys=[goi], layer="imputed_count"
    ).values
    min_ = actual_values.min()
    max_ = actual_values.max()
    range_ = max_ - min_

    upper_limit = range_ * (safe_range_fold - 1) + max_
    return {goi: upper_limit}


def run_simulation(inputs, outputs, wildcards, params):
    oracle = co.load_hdf5(inputs["oracle"])
    gradient = co.load_hdf5(inputs["gradient"])
    df_branches = pd.read_csv(inputs["branches"])
    df_branches = df_branches.merge(
        oracle.adata.obs[["clusterid"]],
        left_on="cell_id",
        right_index=True,
        sort=False,
    ).sort_index()
    goi = wildcards["factor"]

    # Overexpression setup
    perturb_condition = fetch_perturb_condition(oracle, goi)

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
    df_all = dev_module(oracle, gradient)

    def do_branch(branch_):
        cell_ids = df_branches.query("branch == @branch_")["cell_id"].tolist()
        cell_idx = np.where(oracle.adata.obs_names.isin(cell_ids))[0]
        df_branch = dev_module(oracle, gradient, cell_idx=cell_idx).assign(
            branch=branch_
        )
        return df_branch

    # new
    def do_subbranches(branch_):
        def do_subbranch(i, j):
            i, j = str(i), str(j)
            cell_ids = (
                df_branches.query("branch == @branch_")
                .query("clusterid == @i or clusterid == @j")["cell_id"]
                .tolist()
            )
            cell_idx = np.where(oracle.adata.obs_names.isin(cell_ids))[0]
            df_subbranch = dev_module(oracle, gradient, cell_idx=cell_idx).assign(
                branch=branch_,
                subbranch=f"{i}_{j}",
                c_start=i,
                c_end=j,
            )
            return df_subbranch

        clusterids = parse_branch_clusters(branch_)
        df_subbranches = []
        for i, j in zip(clusterids[:-1], clusterids[1:]):
            df_subbranch = do_subbranch(i, j)
            df_subbranches.append(df_subbranch)
        return pd.concat(df_subbranches)

    branches = df_branches["branch"].unique()
    df_ps_branches = [do_branch(branch) for branch in branches]
    df_ps_subbranches = [do_subbranches(branch) for branch in branches]  # new
    df_perturb_scores = pd.concat([df_all, *df_ps_branches, *df_ps_subbranches])

    df_mc_counts = oracle.count_cells_in_mc_resutls(oracle.cluster_column_name)
    df_mc_transitions = oracle.markvov_transition_id

    # Save results
    df_perturb_scores.to_parquet(outputs["perturb_scores"])
    df_mc_counts.to_parquet(outputs["mc_counts"])
    df_mc_transitions.columns = df_mc_transitions.columns.astype(str)
    df_mc_transitions.to_parquet(outputs["mc_transitions"])


snakemake = snakemake  # type: ignore
run_simulation(snakemake.input, snakemake.output, snakemake.wildcards, snakemake.params)
