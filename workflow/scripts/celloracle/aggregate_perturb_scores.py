import pandas as pd
from tqdm import tqdm


def aggregate_perturb_scores(inputs, output_file, params):
    df_list = []
    iterable = list(zip(inputs, params["factor_list"]))
    for input_file, factor in tqdm(iterable):
        df = pd.read_parquet(input_file)
        df["factor"] = factor
        df_list.append(df)

    df = pd.concat(df_list)
    df = df.reset_index(drop=True)
    df.to_parquet(output_file)


snakemake = snakemake  # type: ignore
aggregate_perturb_scores(snakemake.input, snakemake.output[0], snakemake.params)
