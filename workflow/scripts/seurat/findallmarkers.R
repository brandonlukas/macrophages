library(tidyverse)
library(Seurat)

input_cells <- snakemake@input[[1]]
assay <- snakemake@params[["assay"]]
output_file <- snakemake@output[[1]]

# Main script ====
cells <-
  input_cells %>%
  readRDS()

df <- cells %>%
  FindAllMarkers(assay = assay)

write_csv(df, output_file)
