# Revision task T1 / R2-2 — edge-level reproducibility of the joint trajectory in
# the condition-specific ensemble.
#
# The condition-specific draws each have their own clusters (arbitrary labels,
# different counts), so edges cannot be matched by label. Instead we match by the
# CELLS on an edge — exactly the logic Lamian's evaluate_uncertainty uses to score
# branch reproducibility across bootstraps.
#
# For each joint MST edge (a-b), restricted to one condition's cells:
#   Sj = cells whose JOINT cluster is a or b.
# For each draw and each draw edge (x-y):
#   Sd = cells whose DRAW cluster is x or y.
#   Jaccard(Sj, Sd) = |Sj n Sd| / |Sj u Sd|;  overlap = |Sj n Sd| / min(|Sj|,|Sd|).
# The joint edge's score in a draw is the best (max Jaccard) matching draw edge.
# Detection = best Jaccard >= a per-edge null cutoff (99th pct of Jaccard against
# random cell sets of size |Sj|, the same null Lamian uses). Detection RATE is the
# fraction of the 10 draws in which the edge is detected.
#
# Outputs:
#   manifest — per (condition, seed): cluster count, joint pseudotime concordance.
#   edges    — per (condition, joint edge): mean/median/min best Jaccard, mean
#              overlap, null cutoff, detection rate. APC edges (6-7, 5-7) flagged.

library(tidyverse)
library(igraph)

joint_file <- snakemake@input[["joint"]]
draw_files <- unlist(snakemake@input[["draws"]])
manifest_out <- snakemake@output[["manifest"]]
edges_out <- snakemake@output[["edges"]]
n_null <- snakemake@params[["n_null"]]
seed <- snakemake@params[["seed"]]

APC_EDGES <- c("6-7", "5-7")

joint <- readRDS(joint_file)
joint_cl_all <- joint$clusterid
joint_pt <- enframe(joint$pseudotime, name = "cell", value = "pt_joint") %>%
  distinct(cell, .keep_all = TRUE)
jedges <- igraph::as_edgelist(joint$MSTtree)

edge_key <- function(a, b) paste(sort(as.integer(c(a, b))), collapse = "-")
edge_cells <- function(cl, a, b) names(cl)[cl %in% c(a, b)]
jaccard <- function(A, B) length(intersect(A, B)) / length(union(A, B))
overlap <- function(A, B) length(intersect(A, B)) / min(length(A), length(B))
tot_wss <- function(X, cl) sum(vapply(split(as.data.frame(X), cl), function(g) {
  m <- colMeans(g); sum(rowSums(sweep(as.matrix(g), 2, m)^2))
}, numeric(1)))

parse_meta <- function(p) {
  m <- str_match(p, "lamian_condpca_draws/([^/]+)/seed(\\d+)/")
  list(cond = m[, 2], seed = as.integer(m[, 3]))
}
by_cond <- split(draw_files, sapply(draw_files, function(p) parse_meta(p)$cond))

set.seed(seed)
manifest <- list()
edge_rows <- list()

for (cond in names(by_cond)) {
  files <- by_cond[[cond]]
  draws <- lapply(files, readRDS)
  seeds <- sapply(files, function(p) parse_meta(p)$seed)
  cells_c <- names(draws[[1]]$clusterid)
  jcl <- joint_cl_all[cells_c] # joint labels on this condition's cells

  # per-draw manifest (cluster count + joint pseudotime concordance)
  for (i in seq_along(draws)) {
    dr <- draws[[i]]
    cp <- enframe(dr$pseudotime, name = "cell", value = "pt") %>% distinct(cell, .keep_all = TRUE)
    tab <- inner_join(cp, joint_pt, by = "cell")
    manifest[[length(manifest) + 1]] <- tibble(
      condition = cond, seed = seeds[i],
      n_clusters = length(unique(dr$clusterid)), pcadim = ncol(dr$pca),
      rho_joint = cor(tab$pt, tab$pt_joint, method = "spearman"),
      tot_wss = tot_wss(dr$pca[names(dr$clusterid), , drop = FALSE], dr$clusterid))
  }

  # precompute each draw's edge cell-sets
  draw_edge_sets <- lapply(draws, function(dr) {
    de <- igraph::as_edgelist(dr$MSTtree)
    if (nrow(de) == 0) return(list())
    lapply(seq_len(nrow(de)), function(r) edge_cells(dr$clusterid, de[r, 1], de[r, 2]))
  })

  for (i in seq_len(nrow(jedges))) {
    a <- jedges[i, 1]; b <- jedges[i, 2]
    Sj <- edge_cells(jcl, a, b)
    if (length(Sj) == 0) next
    # Lamian-style null: Jaccard of Sj vs random |Sj| cells; cutoff = 99th pct.
    null_j <- replicate(n_null, jaccard(Sj, sample(cells_c, length(Sj))))
    cutoff <- as.numeric(quantile(null_j, 0.99))
    # best match per draw
    best_j <- sapply(draw_edge_sets, function(es) if (length(es)) max(vapply(es, jaccard, numeric(1), B = Sj)) else 0)
    best_o <- sapply(draw_edge_sets, function(es) if (length(es)) max(vapply(es, overlap, numeric(1), B = Sj)) else 0)
    edge_rows[[length(edge_rows) + 1]] <- tibble(
      condition = cond, edge = edge_key(a, b), is_APC = edge_key(a, b) %in% APC_EDGES,
      n_joint_cells = length(Sj), null_cutoff = cutoff,
      jacc_mean = mean(best_j), jacc_median = median(best_j), jacc_min = min(best_j),
      overlap_mean = mean(best_o),
      n_detected = sum(best_j >= cutoff), n_draws = length(draws),
      detection_rate = mean(best_j >= cutoff))
  }
}

manifest_tbl <- bind_rows(manifest) %>% arrange(condition, seed)
edge_tbl <- bind_rows(edge_rows) %>% arrange(condition, desc(is_APC), desc(detection_rate))
write_csv(manifest_tbl, manifest_out)
write_csv(edge_tbl, edges_out)

cat("=== per-draw manifest ===\n"); print(manifest_tbl, n = Inf, width = Inf)
cat("\n=== edge reproducibility (APC edges first) ===\n")
edge_tbl %>% mutate(across(c(jacc_mean, jacc_median, jacc_min, overlap_mean, null_cutoff), ~round(., 3))) %>%
  print(n = Inf, width = Inf)
