import celloracle as co
from celloracle.applications import Gradient_calculator


def prepare_gradient(input_file, output_file, params):
    oracle = co.load_hdf5(input_file)
    gradient = Gradient_calculator(
        oracle_object=oracle,
        pseudotime_key="pseudotime",
    )
    gradient.calculate_p_mass(n_grid=params["n_grid"])
    gradient.calculate_mass_filter(min_mass=params["min_mass"])
    gradient.transfer_data_into_grid()
    gradient.calculate_gradient()
    gradient.to_hdf5(output_file)


snakemake = snakemake  # type: ignore
prepare_gradient(snakemake.input[0], snakemake.output[0], snakemake.params)

# input_file = "results/celloracle/cells.celloracle.oracle"
# params = {
#     "n_grid": 60,
#     "min_mass": 8,
# }
# oracle = co.load_hdf5(input_file)

# import scanpy as sc
# import matplotlib.pyplot as plt

# goi = "Zbtb46"
# sc.tl.pca(oracle.adata)
# sc.pp.neighbors(oracle.adata)
# sc.tl.draw_graph(oracle.adata)
# sc.pl.draw_graph(
#     oracle.adata,
#     color=[goi, oracle.cluster_column_name],
#     layer="imputed_count",
#     use_raw=False,
#     cmap="viridis",
# )
# plt.savefig("test.png", dpi=300, bbox_inches="tight")

# oracle.calculate_p_mass(n_grid=60)
# oracle.suggest_mass_thresholds(n_suggestion=12)

# gradient = Gradient_calculator(
#     oracle_object=oracle,
#     pseudotime_key="pseudotime",
# )
# gradient.calculate_p_mass(n_grid=params["n_grid"])
# gradient.calculate_mass_filter(min_mass=params["min_mass"])
# gradient.transfer_data_into_grid()
# gradient.calculate_gradient()
