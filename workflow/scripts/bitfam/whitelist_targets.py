import pandas as pd


def whitelist_targets(input_file, output_file, params, wildcards):
    df = (
        pd.read_csv(input_file, sep="\t", index_col=0)
        .filter(regex=params.whitelist_regex)
        .assign(source=wildcards.factor, target=lambda x: x.index)
        .melt(id_vars=["source", "target"], var_name="evidence", value_name="weight")
        .query("weight > 0")
        .groupby(by=["source", "target"])
        .agg({"weight": "sum", "evidence": list})
        .reset_index()
    )
    df.to_csv(output_file, index=False)


snakemake = snakemake  # type: ignore
whitelist_targets(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.params,
    snakemake.wildcards,
)
