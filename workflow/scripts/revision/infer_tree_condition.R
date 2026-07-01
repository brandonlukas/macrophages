# Revision task T1 / R2-2 (editor-flagged) — per-condition trajectory inference.
#
# Re-runs Lamian's infer_tree_structure() on a SINGLE condition (wt or db) so we
# can test whether the joint ND+DB trajectory reproduces when each condition is
# inferred independently. Cluster labels here are Lamian's own k-means nodes for
# the subset; they are reconciled back to the joint clusters in the concordance
# step (they are NOT assumed to match by number).
#
# Mirrors workflow/scripts/lamian/infer_tree.R so the only difference from the
# published joint run is the condition subset.

library(tidyverse)
library(Seurat)
library(Lamian)

input_cells <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
condition <- snakemake@wildcards[["cond"]]
seed_1 <- snakemake@params[["seed_1"]]
seed_2 <- snakemake@params[["seed_2"]]

cells <- readRDS(input_cells)

# Subset to the requested condition (wt = non-diabetic, db = diabetic).
keep <- cells$condition == condition
if (!any(keep)) {
  stop(sprintf("No cells with condition == '%s'", condition))
}
cells <- cells[, keep]
message(sprintf("Condition '%s': %d cells", condition, ncol(cells)))

pca <- cells@reductions$pca@cell.embeddings %>%
  as.matrix()

cellanno <- cells@meta.data %>%
  mutate(
    cell = rownames(.),
    sample = orig.ident,
    timepoint = case_when(
      grepl("D10", orig.ident) ~ "D10",
      grepl("D6", orig.ident) ~ "D6",
      grepl("D3", orig.ident) ~ "D3"
    )
  ) %>%
  select(cell, sample, timepoint)

set.seed(seed_1)
res <- infer_tree_structure(
  pca = pca,
  cellanno = cellanno,
  origin.celltype = "D3",
  kmeans.seed = seed_2
)

saveRDS(res, output_file)
