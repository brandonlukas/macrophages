# Revision task T1 / R2-2 — consensus across the condition-PCA seed ensemble.
#
# A single reproducible seed can still be a lucky/unlucky draw, so we infer the
# per-condition trajectory across several k-means seeds and summarise the
# consensus. Outputs:
#   manifest  — one row per (condition, seed): cluster count, Spearman vs joint
#               pseudotime, total within-cluster SS.
#   summary   — one row per condition: rho distribution, modal k, how much the
#               seeds agree (mean pairwise adjusted Rand index), the most central
#               "representative" seed and the best-concordance seed.
#   labels    — co-association consensus clustering (per cell): the fraction of
#               seeds in which each cell pair co-clusters, cut to the modal k.
#
# ND is expected to be a tight consensus; DB is multi-modal (see the doc), so its
# low mean ARI is the point — it quantifies that DB has no single trajectory.

library(tidyverse)

joint_file <- snakemake@input[["joint"]]
seed_files <- unlist(snakemake@input[["seeds"]])
manifest_out <- snakemake@output[["manifest"]]
summary_out <- snakemake@output[["summary"]]
labels_out <- snakemake@output[["labels"]]

joint <- readRDS(joint_file)
joint_pt <- enframe(joint$pseudotime, name = "cell", value = "pt_joint") %>%
  distinct(cell, .keep_all = TRUE)

parse_cond <- function(p) str_match(p, "lamian_condpca_seeds/([^/]+)/seed(\\d+)/")[, 2:3]

# Adjusted Rand index between two hard clusterings (no external dependency).
adj_rand <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  sum_ij <- sum(choose(tab, 2))
  sa <- sum(choose(rowSums(tab), 2))
  sb <- sum(choose(colSums(tab), 2))
  expected <- sa * sb / choose(n, 2)
  maxidx <- 0.5 * (sa + sb)
  if (maxidx == expected) return(1)
  (sum_ij - expected) / (maxidx - expected)
}

tot_wss <- function(X, cl) {
  sum(vapply(split(as.data.frame(X), cl), function(g) {
    m <- colMeans(g)
    sum(rowSums(sweep(as.matrix(g), 2, m)^2))
  }, numeric(1)))
}

recs <- map(seed_files, function(f) {
  cs <- parse_cond(f)
  r <- readRDS(f)
  cid <- r$clusterid
  cp <- enframe(r$pseudotime, name = "cell", value = "pt") %>%
    distinct(cell, .keep_all = TRUE)
  tab <- inner_join(cp, joint_pt, by = "cell")
  list(cond = cs[1], seed = as.integer(cs[2]), cid = cid,
       k = length(unique(cid)),
       rho = cor(tab$pt, tab$pt_joint, method = "spearman"),
       wss = tot_wss(r$pca[names(cid), , drop = FALSE], cid))
})

manifest <- map_dfr(recs, ~ tibble(condition = .x$cond, seed = .x$seed,
                                   n_clusters = .x$k, rho_joint = .x$rho,
                                   tot_wss = .x$wss)) %>%
  arrange(condition, seed)

summ <- list()
labels <- list()
for (cond in unique(manifest$condition)) {
  rc <- keep(recs, ~ .x$cond == cond)
  seeds <- map_int(rc, "seed")
  rhos <- map_dbl(rc, "rho")
  ks <- map_int(rc, "k")
  cells <- names(rc[[1]]$cid)
  cids <- map(rc, ~ .x$cid[cells])

  # Pairwise ARI -> agreement + most central seed.
  m <- length(rc)
  A <- matrix(1, m, m)
  for (i in seq_len(m)) for (j in seq_len(m)) if (i < j) {
    A[i, j] <- A[j, i] <- adj_rand(cids[[i]], cids[[j]])
  }
  mean_ari <- (rowSums(A) - 1) / (m - 1)

  # Co-association consensus clustering at the modal k.
  modal_k <- as.integer(names(sort(table(ks), decreasing = TRUE))[1])
  co <- matrix(0, length(cells), length(cells))
  for (cl in cids) for (g in split(seq_along(cells), cl)) co[g, g] <- co[g, g] + 1
  co <- co / m
  cons <- cutree(hclust(as.dist(1 - co), method = "average"), k = modal_k)

  labels[[cond]] <- tibble(condition = cond, cell = cells, consensus_cluster = cons)
  summ[[cond]] <- tibble(
    condition = cond, n_seeds = m,
    modal_k = modal_k, k_min = min(ks), k_max = max(ks),
    rho_median = median(rhos), rho_min = min(rhos), rho_max = max(rhos),
    mean_pairwise_ari = mean(A[upper.tri(A)]),
    representative_seed = seeds[which.max(mean_ari)],
    representative_rho = rhos[which.max(mean_ari)],
    best_rho_seed = seeds[which.max(rhos)],
    best_rho = max(rhos))
}

summary_tbl <- bind_rows(summ)
write_csv(manifest, manifest_out)
write_csv(summary_tbl, summary_out)
write_csv(bind_rows(labels), labels_out)
print(summary_tbl, width = Inf)
