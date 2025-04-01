import pandas as pd
from tqdm import tqdm


def aggregate_mc_transitions(inputs, output_file, wildcards, params):
    metadata = pd.read_parquet(inputs["metadata"])

    df_list = []
    iterable = list(zip(inputs["mc_transitions"], params["factor_list"]))
    for input_file, factor in tqdm(iterable):
        idx = pd.read_parquet(input_file)[wildcards["step"]].values
        df = metadata.iloc[idx].copy()
        df["factor"] = factor
        df_list.append(df)

    df = pd.concat(df_list)
    df = df.reset_index(drop=True)
    df.to_parquet(output_file)


snakemake = snakemake  # type: ignore
aggregate_mc_transitions(
    snakemake.input, snakemake.output[0], snakemake.wildcards, snakemake.params
)
