library(tidyverse)
library(Seurat)
library(BITFAM2)

input_cells <- snakemake@input[["cells"]]
input_network <- snakemake@input[["network"]]
output_file <- snakemake@output[[1]]
min_targets <- snakemake@params[["min_targets"]]
seed <- snakemake@wildcards[["run"]]

# Main script ====
cells <- input_cells %>%
  readRDS()

genes <- VariableFeatures(cells)
data <- GetAssayData(cells)[genes, ]

# Filtering for TFs with fewer than 10 qualified targets removes
# H1fl1 (2), Cdk9, Gps2, Kdm1a, Lmo2 (3), and Kdm6b (7)

network <- input_network %>%
  read_csv() %>%
  filter(target %in% genes) %>%
  add_count(source) %>%
  filter(n >= min_targets) %>%
  bind_rows(
    tibble(source = "all_genes", target = genes)
  ) %>%
  mutate(value = 1) %>%
  pivot_wider(id_cols = target, names_from = source, values_fill = 0) %>%
  select(-all_genes) %>%
  column_to_rownames("target")

# NEED rownames(data) and rownames(network) to be EQUAL and ALIGNED
data <- data[rownames(network), ]

res <- BITFAM(data, network, seed = as.integer(seed))
saveRDS(res, output_file)
