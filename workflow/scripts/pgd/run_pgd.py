import scanpy as sc
import torch
import pgdiffusion as pgd
import rpy2.robjects as ro
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.rlike.container import NamedList

input_file = "/mnt/box/data/timkoh/macrophages/from_mac/pgd_cells.h5ad"
adata2 = sc.read_h5ad(input_file)


def run_pgd(anndata_input, lamian_branches, output_file):
    adata = sc.read_h5ad(anndata_input)
    trajectories = load_orders_from_lamian(lamian_branches)

    # Get embeddings (PCA, UMAP, etc.)
    X = torch.tensor(adata.obsm["X_pca"], dtype=torch.float32)

    # Build pseudotime graph
    edge_index = pgd.build_graph(adata, trajectories, neighbors_per_side=50)

    # Apply diffusion
    X_diffused = pgd.diffuse(
        X,
        torch.tensor(edge_index, dtype=torch.long),
        alpha=0.6,
        n_steps=1,
    )

    # Store results
    adata.obsm["X_pseudotime"] = X_diffused.cpu().numpy()

    # Visualize
    sc.pp.neighbors(adata, use_rep="X_pseudotime")
    # sc.tl.umap(adata)
    adata.obsm["X_umap"] = adata2.obsm["X_umap"]  # use precomputed UMAP
    adata.write_h5ad(output_file)  # save


def load_orders_from_lamian(path: str) -> dict:
    # Load Lamian infer_tree result
    with localconverter(default_converter + numpy2ri.converter + pandas2ri.converter):
        result = ro.r(f"readRDS('{path}')")

    orders = _named_list_to_dict(_named_list_to_dict(result)["order"])
    return orders


def _named_list_to_dict(named_list: NamedList) -> dict:
    return dict(
        zip(
            list(map(str, named_list._NamedList__names)),
            list(named_list._NamedList__list),
        )
    )


snakemake = snakemake  # type: ignore
run_pgd(snakemake.input[0], snakemake.input[1], snakemake.output[0])

# Unfortunately the UMAP doesn't look good if done on PC but looks great on my Mac
# It was generated using same code
