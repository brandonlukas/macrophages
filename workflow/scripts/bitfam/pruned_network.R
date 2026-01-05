library(tidyverse)
library(Seurat)

input_cells <- snakemake@input[["cells"]]
input_network <- snakemake@input[["network"]]
output_network <- snakemake@output[[1]]

# Main script ====
cells <- input_cells %>%
  readRDS()

genes <- Features(cells)
expressed_genes <- rowSums(cells[["RNA"]]$data > 0) %>%
  enframe() %>%
  filter(value > 0)

network <- input_network %>%
  read_csv() %>%
  filter(source %in% genes, target %in% genes)

expressed_genes %>%
  filter(name %in% network$source) %>%
  arrange(value)

# Here, we note that 100 cells is a decent cutoff
# This removes Grip1 (1), Gfi1 (10), and Tal1 (82), but keeps Zbtb46 (170)
genes2 <- expressed_genes %>%
  filter(value > 100) %>%
  pull(name)

network <- network %>%
  filter(source %in% genes2)

write_csv(network, output_network)
