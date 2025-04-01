import scanpy as sc
import celloracle as co
import numpy as np
import pandas as pd


def prepare_oracle(input_file, output_file, params):
    # Load data
    adata = sc.read_h5ad(input_file)
    adata.X = adata.layers["counts_RNA"].copy()
    adata.obs["clusterid"] = adata.obs["clusterid"].astype(str)

    # Load the base GRN
    base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()

    # Add Zbtb46 as a TF
    # Meredith et al. (2012); J Exp Med; PMID: 22851594
    zbtb46_targets = list(
        set(
            str(x)
            for x in pd.read_excel(params["zbtb46_targets"]).to_numpy().flatten()
            if isinstance(x, str)
        )
    )
    zbtb46 = pd.DataFrame({"gene_short_name": zbtb46_targets, "Zbtb46": 1.0})
    base_GRN = base_GRN.merge(
        zbtb46,
        how="left",
        on="gene_short_name",
    ).replace(np.nan, 0.0)

    # Initialize Oracle object
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name="clusterid",
        embedding_name="X_umap",
    )
    oracle.import_TF_data(TF_info_matrix=base_GRN)

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

    # GRN calculation
    links = oracle.get_links(cluster_name_for_GRN_unit="clusterid")

    # Make predictive models for simulation
    links.filter_links()
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(use_cluster_specific_TFdict=True)

    # Save oracle object
    oracle.to_hdf5(output_file)


snakemake = snakemake  # type: ignore
prepare_oracle(snakemake.input[0], snakemake.output[0], snakemake.params)
