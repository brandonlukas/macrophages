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
# The best-match cell-Jaccard for a joint edge in one draw is the max over that
# draw's edges — a label-agnostic measure of cell-level edge recovery.
#
# TWO scorings are reported:
#
# (A) Lamian-style detection (unchanged): per-edge single-set null cutoffs
#     js.cut / oc.cut (99th pct of Jaccard / overlap vs random size-|Sj| sets),
#     one-to-one matching via Lamian:::get_binary, detection.rate = (js+oc)/2.
#
# (B) Calibrated significance of the best-match magnitude (NEW). The observed
#     statistic is a MAX over the draw's k edges, so a single-set null is
#     anti-conservative. We use a BEST-OF-k null: the max of k independent
#     size-|Sj| random-set Jaccards. Because a size-matched random set's
#     intersection with Sj is hypergeometric, the single-set Jaccard null is
#     analytic (no sampling), and the best-of-k 99th-pct cutoff is the
#     (0.99^(1/k))-quantile of that single-set null. With this cutoff each draw
#     has a 1% false-positive rate under H0, so the number of the n_draws that
#     recover an edge is Binomial(n_draws, 0.01) under H0. That gives a per-edge
#     combined p-value (p_binom), Benjamini-Hochberg adjusted across all
#     edge x condition tests (p_bh). k is the per-condition median edge count.
#     Caveat: the draws share the embedding and differ only by k-means
#     initialization, so they are not fully independent replicates.
#
# Pure post-processing of the saved draws; does NOT re-infer any trajectory.
#
# Outputs:
#   manifest — per (condition, seed): cluster count, joint pseudotime concordance.
#   edges    — per (condition, joint edge): best Jaccard/overlap magnitude, single
#              -set + best-of-k null cutoffs, js/oc/averaged detection rate,
#              best-of-k detection count, binomial p and BH-adjusted p. APC flagged.

library(tidyverse)
library(igraph)
library(Lamian) # for get_binary (exact Lamian one-to-one matching)
library(reshape2)

APC_EDGES <- c("6-7", "5-7")
NULL_LEVEL <- 0.99 # per-edge / per-draw false-positive rate

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

# Detected joint edges (columns) in one draw for a metric matrix + cutoffs, via
# Lamian's one-to-one matching. mat is [draw-edge x joint-edge].
detected_edges <- function(mat, cut) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(rep(FALSE, ncol(mat)))
  bin <- Lamian:::get_binary(mat, cut)
  if (is.null(dim(bin))) bin <- matrix(bin, nrow = nrow(mat))
  colSums(bin) >= 1
}

# Analytic best-of-k Jaccard cutoff. A size-m random set's intersection I with a
# size-m target is Hypergeometric(white = m, black = N - m, draws = m); Jaccard =
# I / (2m - I) (both sets size m). The 99th pct of max-over-k such Jaccards is the
# (level^(1/k))-quantile of the single-set null (max of k iid: F(t)^k = level).
bok_jaccard_cut <- function(m, N, k, level = NULL_LEVEL) {
  q <- level^(1 / k)
  I <- qhyper(q, m = m, n = N - m, k = m)
  I / (2 * m - I)
}

compute_edge_detection <- function(joint_file, draw_files, n_null = 1000, seed = 42) {
  joint <- readRDS(joint_file)
  joint_cl_all <- joint$clusterid
  joint_pt <- enframe(joint$pseudotime, name = "cell", value = "pt_joint") %>%
    distinct(cell, .keep_all = TRUE)
  jedges <- igraph::as_edgelist(joint$MSTtree)
  by_cond <- split(draw_files, sapply(draw_files, function(p) parse_meta(p)$cond))

  set.seed(seed)
  manifest <- list()
  edge_rows <- list()

  for (cond in names(by_cond)) {
    files <- by_cond[[cond]]
    draws <- lapply(files, readRDS)
    seeds <- sapply(files, function(p) parse_meta(p)$seed)
    n_draws <- length(draws)
    cells_c <- names(draws[[1]]$clusterid)
    N <- length(cells_c)
    jcl <- joint_cl_all[cells_c]

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

    # joint edge cell-sets (columns) + single-set null cutoffs (Lamian-style)
    je <- lapply(seq_len(nrow(jedges)), function(i) edge_cells(jcl, jedges[i, 1], jedges[i, 2]))
    keep <- lengths(je) > 0
    je <- je[keep]
    jkeys <- vapply(which(keep), function(i) edge_key(jedges[i, 1], jedges[i, 2]), character(1))
    js_cut <- vapply(je, function(Sj) quantile(replicate(n_null, jaccard(Sj, sample(cells_c, length(Sj)))), NULL_LEVEL), numeric(1))
    oc_cut <- vapply(je, function(Sj) quantile(replicate(n_null, overlap(Sj, sample(cells_c, length(Sj)))), NULL_LEVEL), numeric(1))

    # draw edge cell-sets; representative k = median edge count across draws
    draw_edges <- lapply(draws, function(dr) {
      de <- igraph::as_edgelist(dr$MSTtree)
      if (nrow(de) == 0) return(list())
      lapply(seq_len(nrow(de)), function(r) edge_cells(dr$clusterid, de[r, 1], de[r, 2]))
    })
    k_rep <- as.integer(round(median(lengths(draw_edges))))
    js_bok_cut <- vapply(je, function(Sj) bok_jaccard_cut(length(Sj), N, k_rep), numeric(1))

    # per-draw best-match magnitudes + Lamian detection
    js_detect <- oc_detect <- matrix(FALSE, n_draws, length(je))
    best_j <- best_o <- matrix(0, n_draws, length(je))
    for (d in seq_along(draws)) {
      es <- draw_edges[[d]]
      if (!length(es)) next
      jm <- sapply(je, function(Sj) vapply(es, jaccard, numeric(1), B = Sj))
      om <- sapply(je, function(Sj) vapply(es, overlap, numeric(1), B = Sj))
      if (is.null(dim(jm))) { jm <- matrix(jm, nrow = length(es)); om <- matrix(om, nrow = length(es)) }
      best_j[d, ] <- apply(jm, 2, max); best_o[d, ] <- apply(om, 2, max)
      js_detect[d, ] <- detected_edges(jm, js_cut)
      oc_detect[d, ] <- detected_edges(om, oc_cut)
    }

    # best-of-k calibrated detection + binomial combined p across draws
    n_detected_bok <- colSums(sweep(best_j, 2, js_bok_cut, `>=`))
    p_binom <- pbinom(n_detected_bok - 1, n_draws, 1 - NULL_LEVEL, lower.tail = FALSE)

    js_perc <- colMeans(js_detect); oc_perc <- colMeans(oc_detect)
    edge_rows[[length(edge_rows) + 1]] <- tibble(
      condition = cond, edge = jkeys, is_APC = jkeys %in% APC_EDGES,
      n_joint_cells = lengths(je),
      jacc_mean = colMeans(best_j), jacc_median = apply(best_j, 2, median), jacc_min = apply(best_j, 2, min),
      overlap_mean = colMeans(best_o),
      js_cut = js_cut, oc_cut = oc_cut,
      js_detect_rate = js_perc, oc_detect_rate = oc_perc,
      detection_rate = pmin((js_perc + oc_perc) / 2, 1),
      k_rep = k_rep, js_bok_cut = js_bok_cut,
      n_detected_bok = n_detected_bok, p_binom = p_binom,
      n_draws = n_draws)
  }

  manifest_tbl <- bind_rows(manifest) %>% arrange(condition, seed)
  # BH-FDR across all edge x condition tests
  edge_tbl <- bind_rows(edge_rows) %>%
    mutate(p_bh = p.adjust(p_binom, "BH")) %>%
    arrange(condition, desc(is_APC), desc(jacc_mean))
  list(manifest = manifest_tbl, edges = edge_tbl)
}

if (exists("snakemake")) {
  res <- compute_edge_detection(
    snakemake@input[["joint"]], unlist(snakemake@input[["draws"]]),
    n_null = snakemake@params[["n_null"]], seed = snakemake@params[["seed"]]
  )
  write_csv(res$manifest, snakemake@output[["manifest"]])
  write_csv(res$edges, snakemake@output[["edges"]])

  cat("=== per-draw manifest ===\n"); print(res$manifest, n = Inf, width = Inf)
  cat("\n=== edge reproducibility (APC first) ===\n")
  res$edges %>%
    mutate(across(c(jacc_mean, jacc_median, jacc_min, overlap_mean, js_cut, oc_cut,
                    js_bok_cut, p_binom, p_bh), ~signif(., 3))) %>%
    print(n = Inf, width = Inf)
}
