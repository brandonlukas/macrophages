from pathlib import Path
import pandas as pd
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.rlike.container import NamedList
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2


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
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

    output_path = Path("results/pyscenic/adjacencies.parquet")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adjacencies.to_parquet(output_path)
