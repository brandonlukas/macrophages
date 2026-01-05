library(tidyverse)
library(Seurat)
library(Lamian)

input_cells <- snakemake@input[["cells"]]
input_res <- snakemake@input[["res"]]
output_file <- snakemake@output[[1]]
assay <- snakemake@params[["assay"]]
test_type <- snakemake@params[["test_type"]]
nfeatures <- snakemake@params[["nfeatures"]]
ncores <- snakemake@threads

# Main script ====
cells <-
  input_cells %>%
  readRDS()

res <-
  input_res %>%
  readRDS()

design <- cells@meta.data %>%
  distinct(orig.ident, condition) %>%
  mutate(
    intercept = 1,
    test = ifelse(condition == "db", 1, 0)
  ) %>%
  select(-condition) %>%
  remove_rownames() %>%
  column_to_rownames("orig.ident")

branch_names <- names(res$order)
res_list <- lapply(branch_names, function(branch_name) {
  cell_names <- res$order[[branch_name]]
  cell_subset <- subset(cells, cells = cell_names)

  expr <- cell_subset[[assay]]$data
  if (assay == "RNA") {
    top_genes <- cell_subset %>%
      FindVariableFeatures(nfeatures = nfeatures) %>%
      VariableFeatures()
    expr <- expr[top_genes, ]
  }
  expr <- expr[Matrix::rowSums(expr) != 0, ]

  cellanno <-
    cell_subset@meta.data %>%
    mutate(
      cell = rownames(.),
      sample = orig.ident
    ) %>%
    select(cell, sample)

  tryCatch(
    {
      lamian_test(
        expr = as.matrix(expr),
        cellanno = cellanno,
        pseudotime = cell_subset$pseudotime,
        design = design,
        test.type = test_type,
        ncores = ncores,
      )
    },
    error = function(e) {
      message(
        "Error in lamian_test for branch: ", branch_name, "\n",
        "Error: ", e$message
      )
      NULL
    }
  )
})

names(res_list) <- branch_names
saveRDS(res_list, output_file)
