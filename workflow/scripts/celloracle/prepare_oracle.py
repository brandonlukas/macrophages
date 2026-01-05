import scanpy as sc
import celloracle as co
import numpy as np
import pandas as pd


def prepare_oracle(inputs, output_file):
    # Load data
    adata = sc.read_h5ad(inputs[0])
    adata.X = adata.layers["counts_RNA"].copy()
    adata.obs["clusterid"] = adata.obs["clusterid"].astype(str)

    # Load the base GRN
    base_GRN = pd.read_csv(inputs[1])
    TF_dict = base_GRN.groupby("target")["source"].agg(list).to_dict()

    # Initialize Oracle object
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name="clusterid",
        embedding_name="X_umap",
    )
    oracle.import_TF_data(TFdict=TF_dict)

    # KNN imputation
    oracle.perform_PCA()
    oracle.knn_imputation()

    # GRN calculation
    links = oracle.get_links(cluster_name_for_GRN_unit="clusterid")

    # Make predictive models for simulation
    links.filter_links()
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(use_cluster_specific_TFdict=True)

    # Save oracle object
    oracle.to_hdf5(output_file)


def _knn_imputation_tutorial_style(oracle: co.Oracle) -> co.Oracle:
    # KNN imputation
    oracle.perform_PCA()
    n_comps = np.where(
        np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002)
    )[0][0]
    n_comps = min(n_comps, 50)
    n_cell = oracle.adata.shape[0]
    k = int(0.025 * n_cell)
    oracle.knn_imputation(
        n_pca_dims=n_comps,
        k=k,
        balanced=True,
        b_sight=k * 8,
        b_maxl=k * 4,
        n_jobs=4,
    )
    return oracle


snakemake = snakemake  # type: ignore
prepare_oracle(snakemake.input, snakemake.output[0])
