library(tidyverse)
library(Seurat)

cells <-
  "resources/all_cells.rds" %>%
  readRDS()

DimPlot(cells, group.by = "cell_type")
cells <- CreateSeuratObject(
  counts = cells[["RNA"]]@counts,
  meta.data = cells@meta.data
)

cells <- cells %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

DimPlot(cells, group.by = "cell_type")

dir.create("inputs", showWarnings = FALSE)
saveRDS(cells, "inputs/cells.rds")
