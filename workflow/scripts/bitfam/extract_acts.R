library(tidyverse)
library(Seurat)
library(BITFAM2)

input_cells <- snakemake@input[["cells"]]
input_res <- snakemake@input[["res"]]
output_file <- snakemake@output[[1]]

# Main script ====
cells <- input_cells %>%
  readRDS()

res <- input_res %>%
  readRDS()

acts <- BITFAM_activities(res)
cells[["Z"]] <- CreateAssayObject(data = t(acts))

saveRDS(cells, output_file)
