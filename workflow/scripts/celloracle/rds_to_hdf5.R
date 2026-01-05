library(tidyverse)
library(Seurat)
library(scCustomize)
library(reticulate)

input_file <- snakemake@input[[1]]
network_file <- snakemake@input[[2]]
output_file <- snakemake@output[[1]]
conda_env <- snakemake@params[["conda_env"]]

use_condaenv(conda_env, required = TRUE)

cells <- readRDS(input_file)

grn <- read_csv(network_file)
tf_list <- grn %>%
  pull(source) %>%
  unique()

# CellOracle recommends top 2000-3000 variable genes
cells <- FindVariableFeatures(cells, nfeatures = 3000)
genes <- VariableFeatures(cells)
genes <- c(genes, tf_list)

cells_diet <- DietSeurat(
  cells,
  layers = c("counts", "data"),
  features = genes,
  dimreducs = c("pca", "umap")
)
cells_diet <- UpdateSeuratObject(cells_diet)

dir_path <- dirname(output_file)
file_name <- basename(output_file)
as.anndata(
  x = cells_diet,
  file_path = dir_path,
  file_name = file_name
)
