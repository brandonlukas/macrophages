library(tidyverse)
library(Seurat)
library(scCustomize)
library(reticulate)

use_condaenv("celloracle_env", required = TRUE)

cells <- readRDS("resources/all_cells.rds")
res <- readRDS("resources/res.rds")

# Reformat pseudotime to remove duplicates
res$pseudotime <- res[["pseudotime"]] %>%
  enframe() %>%
  distinct() %>%
  deframe()

# Add pseudotime and clusterid to cells
cells$pseudotime <- res$pseudotime[Cells(cells)]
cells$clusterid <- res$clusterid[Cells(cells)]

genes <- VariableFeatures(cells)
genes_of_interest <- c("Zbtb46")
genes <- c(genes, genes_of_interest)

cells_diet <- DietSeurat(
  cells,
  layers = c("counts", "data"),
  features = genes,
  dimreducs = "umap"
)
cells_diet <- UpdateSeuratObject(cells_diet)

as.anndata(
  x = cells_diet,
  file_path = "./resources",
  file_name = "all_cells.h5ad"
)
