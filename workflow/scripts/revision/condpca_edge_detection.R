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
# For each draw edge (x-y):  Sd = cells whose DRAW cluster is x or y.
#   Jaccard(Sj,Sd) = |Sj n Sd|/|Sj u Sd|;  overlap(Sj,Sd) = |Sj n Sd|/min(|Sj|,|Sd|).
#
# Detection follows Lamian::evaluate_uncertainty exactly:
#   - per-edge null cutoffs js.cut / oc.cut = 99th percentile of Jaccard / overlap
#     against random cell sets of size |Sj| (Lamian's js.null/oc.null construction);
#   - within each draw, threshold the [draw-edge x joint-edge] Jaccard and overlap
#     matrices and resolve to a ONE-TO-ONE matching with Lamian:::get_binary;
#   - a joint edge is "detected" in a draw if it gets a match; js.perc / oc.perc are
#     the fraction of draws detected under each metric;
#   - detection.rate = (js.perc + oc.perc) / 2   (Lamian's average of the two).
# We also report the best-match Jaccard/overlap MAGNITUDE across draws — the more
# informative signal here, since detection is saturated (edges beat the null in
# every draw). APC edges (6-7, 5-7) are flagged.
#
# Pure post-processing of the saved draws; does NOT re-infer any trajectory.
#
# Outputs:
#   manifest — per (condition, seed): cluster count, joint pseudotime concordance.
#   edges    — per (condition, joint edge): best Jaccard/overlap magnitude, null
#              cutoffs, js/oc/averaged detection rate. APC edges flagged.

library(tidyverse)
library(igraph)
library(Lamian) # for get_binary (exact Lamian one-to-one matching)
library(reshape2)

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

# Detected joint edges (columns) in one draw, for a given metric matrix + cutoffs.
# mat is [draw-edge x joint-edge]; get_binary thresholds each column (joint edge)
# by cut[col] and enforces a one-to-one match. Returns logical over joint edges.
detected_edges <- function(mat, cut) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(rep(FALSE, ncol(mat)))
  bin <- Lamian:::get_binary(mat, cut)
  if (is.null(dim(bin))) bin <- matrix(bin, nrow = nrow(mat)) # 1-col safety
  colSums(bin) >= 1
}

set.seed(seed)
manifest <- list()
edge_rows <- list()

for (cond in names(by_cond)) {
  files <- by_cond[[cond]]
  draws <- lapply(files, readRDS)
  seeds <- sapply(files, function(p) parse_meta(p)$seed)
  cells_c <- names(draws[[1]]$clusterid)
  jcl <- joint_cl_all[cells_c]

  # per-draw manifest
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

  # joint edge cell-sets (columns, fixed order) + per-edge null cutoffs
  je <- lapply(seq_len(nrow(jedges)), function(i) edge_cells(jcl, jedges[i, 1], jedges[i, 2]))
  keep <- lengths(je) > 0
  je <- je[keep]
  jkeys <- vapply(which(keep), function(i) edge_key(jedges[i, 1], jedges[i, 2]), character(1))
  js_cut <- vapply(je, function(Sj) quantile(replicate(n_null, jaccard(Sj, sample(cells_c, length(Sj)))), 0.99), numeric(1))
  oc_cut <- vapply(je, function(Sj) quantile(replicate(n_null, overlap(Sj, sample(cells_c, length(Sj)))), 0.99), numeric(1))

  # draw edge cell-sets
  draw_edges <- lapply(draws, function(dr) {
    de <- igraph::as_edgelist(dr$MSTtree)
    if (nrow(de) == 0) return(list())
    lapply(seq_len(nrow(de)), function(r) edge_cells(dr$clusterid, de[r, 1], de[r, 2]))
  })

  # per-draw [draw-edge x joint-edge] Jaccard & overlap matrices
  js_detect <- oc_detect <- matrix(FALSE, length(draws), length(je))
  best_j <- best_o <- matrix(0, length(draws), length(je))
  for (d in seq_along(draws)) {
    es <- draw_edges[[d]]
    if (!length(es)) next
    jm <- sapply(je, function(Sj) vapply(es, jaccard, numeric(1), B = Sj)) # [draw-edge x joint-edge]
    om <- sapply(je, function(Sj) vapply(es, overlap, numeric(1), B = Sj))
    if (is.null(dim(jm))) { jm <- matrix(jm, nrow = length(es)); om <- matrix(om, nrow = length(es)) }
    best_j[d, ] <- apply(jm, 2, max); best_o[d, ] <- apply(om, 2, max)
    js_detect[d, ] <- detected_edges(jm, js_cut)
    oc_detect[d, ] <- detected_edges(om, oc_cut)
  }

  js_perc <- colMeans(js_detect); oc_perc <- colMeans(oc_detect)
  edge_rows[[length(edge_rows) + 1]] <- tibble(
    condition = cond, edge = jkeys, is_APC = jkeys %in% APC_EDGES,
    n_joint_cells = lengths(je),
    jacc_mean = colMeans(best_j), jacc_median = apply(best_j, 2, median), jacc_min = apply(best_j, 2, min),
    overlap_mean = colMeans(best_o),
    js_cut = js_cut, oc_cut = oc_cut,
    js_detect_rate = js_perc, oc_detect_rate = oc_perc,
    detection_rate = pmin((js_perc + oc_perc) / 2, 1),
    n_draws = length(draws))
}

manifest_tbl <- bind_rows(manifest) %>% arrange(condition, seed)
edge_tbl <- bind_rows(edge_rows) %>% arrange(condition, desc(is_APC), desc(jacc_mean))
write_csv(manifest_tbl, manifest_out)
write_csv(edge_tbl, edges_out)

cat("=== per-draw manifest ===\n"); print(manifest_tbl, n = Inf, width = Inf)
cat("\n=== edge reproducibility (APC first) ===\n")
edge_tbl %>% mutate(across(c(jacc_mean, jacc_median, jacc_min, overlap_mean, js_cut, oc_cut), ~round(., 3))) %>%
  print(n = Inf, width = Inf)
