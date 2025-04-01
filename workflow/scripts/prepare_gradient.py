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
    gradient.transfer_data_into_grid(args={"method": "knn", "n_knn": 50})
    gradient.calculate_gradient()
    gradient.to_hdf5(output_file)


snakemake = snakemake  # type: ignore
prepare_gradient(snakemake.input[0], snakemake.output[0], snakemake.params)
