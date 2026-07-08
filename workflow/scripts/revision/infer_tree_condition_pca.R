# Revision task T1 (variant C, "most independent") — per-condition trajectory in a
# CONDITION-SPECIFIC PC space. One draw of the ensemble (see revision_condpca_draw):
# we run this repeatedly over {seed} and report the spread honestly.
#
# Recompute the whole embedding on one condition's cells — HVGs, scaling, PCA —
# then run Lamian natively (elbow picks pcadim, kmeans re-clusters). It answers the
# strongest form of the reviewer's question: does the trajectory exist when NOTHING
# is borrowed from the joint analysis?
#
# Reproducibility caveat (stated to the reviewer, NOT patched here): Lamian's
# infer_tree_structure clustering is not reproducible — its internal mykmeans
# ignores the `kmeans.seed` argument (there is no set.seed in its body) and selects
# the cluster number inside an unseeded parallel fork. We pass seeds exactly as the
# main pipeline does (workflow/scripts/lamian/infer_tree.R) so the intent is on the
# record, but we deliberately do NOT modify Lamian's source. Consequently each draw
# is an independent sample; robustness is reported across the ensemble of draws
# (per-draw metrics + edge-level cell-Jaccard reproducibility), not from one tree.
#
# Mirrors workflow/scripts/lamian/infer_tree.R; the only additions are the
# condition subset, the per-condition re-embedding, and saving the full PCA.

library(tidyverse)
library(Seurat)
library(Lamian)

input_cells <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
condition <- snakemake@wildcards[["cond"]]
seed_1 <- snakemake@params[["seed_1"]]
seed_2 <- snakemake@params[["seed_2"]]
npcs <- snakemake@params[["npcs"]]

cells <- readRDS(input_cells)

sub <- cells[, cells$condition == condition]
message(sprintf("Condition '%s': %d cells — recomputing condition-specific PCA", condition, ncol(sub)))

# Condition-specific embedding: HVGs / scaling / PCA all fit on this condition
# alone (this is what makes the PC space condition-specific, per the plan).
DefaultAssay(sub) <- "RNA"
sub <- sub %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = npcs, verbose = FALSE)

pca <- sub@reductions$pca@cell.embeddings %>%
  as.matrix()

cellanno <- sub@meta.data %>%
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

# Seeds passed exactly as the main pipeline does. NB: kmeans.seed is a dead Lamian
# parameter (see caveat above); {seed} still varies the draw via ambient RNG.
set.seed(seed_1)
res <- infer_tree_structure(
  pca = pca,
  cellanno = cellanno,
  origin.celltype = "D3",
  kmeans.seed = seed_2
)

# infer_tree_structure stores only the elbow-truncated PCA in res$pca (ncol ==
# pcadim). Keep the full npcs-dim condition-specific embedding too (deterministic;
# used e.g. for PGD on the condition's PC space).
res$pca_full <- pca

saveRDS(res, output_file)
