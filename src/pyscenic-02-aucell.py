import os
import glob
import pickle
from pathlib import Path
import pandas as pd
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.rlike.container import NamedList
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell


def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]


def namedlist_to_dict(named_list: NamedList) -> dict:
    """
    Convert an rpy2 NamedList into a standard Python dictionary.
    """
    # NamedList supports tuple unpacking: (names, values)
    names = named_list._NamedList__names
    values = [value for value in named_list]
    return {str(name): value for name, value in zip(names, values)}


def load_expression_data(filename: str | Path) -> pd.DataFrame:
    """
    Load expression data from an RDS file using rpy2 and convert it to a pandas DataFrame.

    Parameters:
    filename (str or Path): Path to the RDS file containing the Seurat object.

    Returns:
    pd.DataFrame: DataFrame with genes as rows and cells as columns.
    """
    script = """
local({{
  library(Seurat)
  cells <- readRDS("{filename}")
  list(
      counts = cells@assays[["RNA"]]@counts,
      cells = colnames(cells),
      genes = rownames(cells)
  )
}})
"""
    with localconverter(
        default_converter
        + numpy2ri.converter
        + pandas2ri.converter
        + anndata2ri.converter
    ):
        results = r(script.format(filename=filename))

    results = namedlist_to_dict(results)
    df = pd.DataFrame(
        results["counts"].toarray(), index=results["genes"], columns=results["cells"]
    )
    return df


if __name__ == "__main__":
    ex_matrix = load_expression_data("resources/all_cells.rds").T
    tf_names = load_tf_names("resources/scenic_resources/allTFs_mm.txt")
    db_fnames = glob.glob(
        "resources/scenic_resources/mm10*.genes_vs_motifs.rankings.feather"
    )

    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    adjacencies = pd.read_parquet("results/pyscenic/adjacencies.parquet")
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(
            dbs,
            modules,
            "resources/scenic_resources/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
        )

    # This is necessary for some reason ----
    output_path = Path("results/pyscenic/aucell")
    output_path.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path / "motifs.csv")

    # ----
    df = load_motifs(output_path / "motifs.csv")

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)
    with open(output_path / "regulons.p", "wb") as f:
        pickle.dump(regulons, f)

    # Compute the AUCell matrix (cells x regulons).
    auc_mtx = aucell(ex_matrix, regulons, seed=42)
    auc_mtx.to_parquet(output_path / "auc_mtx.parquet")
