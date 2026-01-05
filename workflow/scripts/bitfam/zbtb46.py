import pandas as pd


def zbtb46(input_file, output_file):
    targets = list(
        set(
            str(x)
            for x in pd.read_excel(input_file).values.flatten()
            if isinstance(x, str)
        )
    )
    df = pd.DataFrame(
        {
            "source": "Zbtb46",
            "target": targets,
            "weight": 1,
            "evidence": "PMID: 22851594",
        }
    )
    df.to_csv(output_file, index=False)


snakemake = snakemake  # type: ignore
zbtb46(snakemake.input[0], snakemake.output[0])
