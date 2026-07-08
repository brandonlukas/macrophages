# Revision task T1 / R2-2 — per-condition trajectory STABILITY sweep.
#
# The condition-specific inference is reproducible once mykmeans is seeded
# (reproducible_kmeans.R), but reproducible != identifiable. This sweep runs the
# same fixed, well-optimized inference across many k-means seeds and records, per
# condition and seed: the cluster count, the concordance with the joint
# pseudotime (Spearman on shared cells), and the total within-cluster SS (the
# clustering objective). It quantifies how much the recovered trajectory depends
# on which near-optimal clustering is found.
#
# Finding (see docs/revision-trajectory-robustness.md): ND (wt) converges to a
# stable solution (rho ~= 0.76); DB (db) is non-identifiable — near-equal-WSS
# clusterings give rho anywhere from ~0.16 to ~0.74, so DB does not support a
# stable independent trajectory.

library(tidyverse)
library(Seurat)
library(Lamian)

source("workflow/scripts/revision/reproducible_kmeans.R")

input_cells <- snakemake@input[["cells"]]
joint_file <- snakemake@input[["joint"]]
output_file <- snakemake@output[[1]]
seed_1 <- snakemake@params[["seed_1"]]
npcs <- snakemake@params[["npcs"]]
nstart <- snakemake@params[["nstart"]]
n_seed <- snakemake@params[["n_seed"]]

install_reproducible_kmeans(nstart)

cells <- readRDS(input_cells)
joint <- readRDS(joint_file)
joint_pt <- enframe(joint$pseudotime, name = "cell", value = "pt_joint") %>%
  distinct(cell, .keep_all = TRUE)

# Total within-cluster SS of a partition in the space kmeans clustered on.
tot_wss <- function(X, cl) {
  sum(vapply(split(as.data.frame(X), cl), function(g) {
    m <- colMeans(g)
    sum(rowSums(sweep(as.matrix(g), 2, m)^2))
  }, numeric(1)))
}

embed <- function(cond) {
  sub <- cells[, cells$condition == cond]
  DefaultAssay(sub) <- "RNA"
  sub <- sub %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  list(
    pca = as.matrix(sub@reductions$pca@cell.embeddings),
    cellanno = sub@meta.data %>%
      mutate(cell = rownames(.), sample = orig.ident,
             timepoint = case_when(
               grepl("D10", orig.ident) ~ "D10",
               grepl("D6", orig.ident) ~ "D6",
               grepl("D3", orig.ident) ~ "D3")) %>%
      select(cell, sample, timepoint))
}

rows <- list()
for (cond in c("wt", "db")) {
  e <- embed(cond)
  for (km_seed in seq_len(n_seed)) {
    set.seed(seed_1)
    res <- infer_tree_structure(pca = e$pca, cellanno = e$cellanno,
                                origin.celltype = "D3", kmeans.seed = km_seed)
    cid <- res$clusterid
    cp <- enframe(res$pseudotime, name = "cell", value = "pt_cond") %>%
      distinct(cell, .keep_all = TRUE)
    tab <- inner_join(cp, joint_pt, by = "cell")
    rows[[length(rows) + 1]] <- tibble(
      condition = cond, kmeans_seed = km_seed,
      n_clusters = length(unique(cid)), pcadim = ncol(res$pca),
      n_shared = nrow(tab),
      rho_joint = cor(tab$pt_cond, tab$pt_joint, method = "spearman"),
      tot_wss = tot_wss(res$pca[names(cid), , drop = FALSE], cid))
    message(sprintf("%s seed=%d: k=%d rho=%.3f wss=%.1f",
                    cond, km_seed, length(unique(cid)),
                    tail(rows, 1)[[1]]$rho_joint, tail(rows, 1)[[1]]$tot_wss))
  }
}

out <- bind_rows(rows)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(out, output_file)

# Console summary (per condition): the honest headline numbers.
out %>%
  group_by(condition) %>%
  summarise(n = n(), rho_median = median(rho_joint),
            rho_min = min(rho_joint), rho_max = max(rho_joint),
            k_min = min(n_clusters), k_max = max(n_clusters), .groups = "drop") %>%
  print()
