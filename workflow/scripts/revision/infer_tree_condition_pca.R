# Revision task T1 (variant C, "most independent") — per-condition trajectory in a
# CONDITION-SPECIFIC PC space.
#
# This recomputes the whole embedding on one condition's cells — HVGs, scaling,
# PCA — then runs Lamian natively (elbow picks pcadim, kmeans reclusters). It
# answers the strongest form of the reviewer's question: does the trajectory
# exist when NOTHING is borrowed from the joint analysis?
#
# Consequence (acknowledged limitation): the PC space, cluster labels and
# pseudotime are condition-specific and NOT directly comparable across conditions
# or to the joint — comparison is only via shared-cell pseudotime rank
# concordance and branch topology (src/revision/trajectory_concordance.R).
#
# Mirrors workflow/scripts/lamian/infer_tree.R; the only additions are the
# condition subset and the per-condition re-embedding.

library(tidyverse)
library(Seurat)
library(Lamian)

# Make Lamian's internal k-means reproducible (its shipped mykmeans ignores its
# seed and forks an unseeded cluster-number search). Without this the trajectory
# is not reproducible run-to-run. See the helper for the full rationale.
source("workflow/scripts/revision/reproducible_kmeans.R")

input_cells <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
condition <- snakemake@wildcards[["cond"]]
seed_1 <- snakemake@params[["seed_1"]]
seed_2 <- snakemake@params[["seed_2"]]
npcs <- snakemake@params[["npcs"]]
nstart <- snakemake@params[["nstart"]]

install_reproducible_kmeans(nstart)

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

set.seed(seed_1)
res <- infer_tree_structure(
  pca = pca,
  cellanno = cellanno,
  origin.celltype = "D3",
  kmeans.seed = seed_2
)

# infer_tree_structure truncates to the elbow-selected dims and stores only those
# in res$pca (ncol == pcadim). Keep the full npcs-dim input embedding too, so the
# complete condition-specific PC space is recoverable downstream.
res$pca_full <- pca

saveRDS(res, output_file)
