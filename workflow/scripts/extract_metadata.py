import anndata as ad
import pandas as pd


def extract_metadata(input_file, output_file):
    adata = ad.read_h5ad(input_file, backed="r")
    df = (
        pd.DataFrame(
            adata.obsm["X_umap"],
            index=adata.obs_names,
            columns=["umap_1", "umap_2"],
        )
        .merge(adata.obs, left_index=True, right_index=True)
        .reset_index(drop=True)
    )
    df.to_parquet(output_file)


snakemake = snakemake  # type: ignore
extract_metadata(snakemake.input[0], snakemake.output[0])
