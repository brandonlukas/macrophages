library(tidyverse)
library(Seurat)

input_cells <- snakemake@input[[1]]
assay <- snakemake@params[["assay"]]
output_file <- snakemake@output[[1]]

# Main script ====
cells <-
  input_cells %>%
  readRDS()

DefaultAssay(cells) <- assay
df <-
  lapply(levels(cells$clusterid), function(clusterid) {
    cells %>%
      FindMarkers(
        ident.1 = "db",
        group.by = "condition",
        subset.ident = clusterid
      ) %>%
      mutate(cluster = clusterid) %>%
      rownames_to_column("gene")
  }) %>%
  bind_rows()

write_csv(df, output_file)
