library(tidyverse)
library(Seurat)

input_cells <- snakemake@input[["cells"]]
input_res <- snakemake@input[["res"]]
input_acts <- snakemake@input[["aggregate_acts"]]
output_file <- snakemake@output[[1]]

# Main script ====
cells <-
  input_cells %>%
  readRDS()

res <-
  input_res %>%
  readRDS()

df <-
  res$pseudotime %>%
  enframe(value = "pseudotime") %>%
  distinct() %>%
  left_join(
    res$clusterid %>%
      enframe(value = "clusterid") %>%
      mutate(clusterid = fct_inseq(as.character(clusterid))),
    by = "name"
  )

metadata <-
  cells@meta.data %>%
  merge(df, by.x = 0, by.y = "name") %>%
  column_to_rownames("Row.names")

cells <-
  cells %>%
  AddMetaData(metadata)

Idents(cells) <- "clusterid"

# -------------------------------------------------------------
# add in aggregate ACT scores
acts <- input_acts %>%
  read_csv() %>%
  group_by(gene, name) %>%
  summarise_if(is.numeric, median)

mat <- acts %>%
  pivot_wider(names_from = name, values_from = value) %>% # nolint
  column_to_rownames("gene") %>%
  as.matrix()

colnames(mat) <- gsub("\\.1", "-1", colnames(mat))
mat <- mat[, colnames(cells)]
cells[["Z_avg"]] <- CreateAssayObject(data = mat)
# -------------------------------------------------------------

saveRDS(cells, output_file)
