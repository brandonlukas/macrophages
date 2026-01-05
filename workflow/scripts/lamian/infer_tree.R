library(tidyverse)
library(Seurat)
library(Lamian)

input_cells <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
seed_1 <- snakemake@params[["seed_1"]]
seed_2 <- snakemake@params[["seed_2"]]

cells <- input_cells %>%
  readRDS()

pca <- cells@reductions$pca@cell.embeddings %>%
  as.matrix()

cellanno <- cells@meta.data %>%
  mutate(
    cell = rownames(.),
    sample = orig.ident,
    timepoint = case_when(
      grepl("D3", orig.ident) ~ "D3",
      grepl("D6", orig.ident) ~ "D6",
      grepl("D10", orig.ident) ~ "D10"
    )
  ) %>%
  select(cell, sample, timepoint)

# Can update to use seurat clusters
# clusters <- as.numeric(cells$seurat_clusters) # nolint: commented_code_linter.
# names(clusters) <- names(cells$seurat_clusters) # nolint:

set.seed(seed_1)
res <- infer_tree_structure(
  pca = pca,
  cellanno = cellanno,
  origin.celltype = "D3",
  kmeans.seed = seed_2
)

saveRDS(res, output_file)
