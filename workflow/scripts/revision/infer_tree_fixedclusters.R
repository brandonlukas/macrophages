# Revision task T1 (variant A, PI/collaborator request) — per-condition trajectory
# with the JOINT clusters held FIXED.
#
# Instead of reclustering each condition (which relabels nodes and made the joint
# vs condition comparison hard), we reuse the joint tree's cluster labels, subset
# to one condition's cells, and refit only the cluster-center minimum spanning
# tree. Because the cluster identities are shared, the resulting trees are
# directly comparable to the joint and to each other (same 1..12 node labels).
#
# The MST is forced to span all clusters, so the point estimate alone is not
# evidence. The test is the bootstrap DETECTION RATE: resample the condition's
# cells, recompute cluster centers with the SAME fixed labels, rebuild the MST,
# and measure how often each edge recurs. A stable edge recurs; a forced/spurious
# edge does not. We report detection for the condition's own MST edges and,
# crucially, for the JOINT edges (does this condition reconstruct them?).
#
# NOTE this deliberately differs from Lamian::evaluate_uncertainty, which
# re-kmeans each bootstrap; here clusters are fixed by design.

library(tidyverse)
library(Seurat)
library(Lamian)
library(igraph)

input_cells <- snakemake@input[["cells"]]
input_joint <- snakemake@input[["joint"]]
out_tree <- snakemake@output[["tree"]]
out_edges <- snakemake@output[["edges"]]
condition <- snakemake@wildcards[["cond"]]
n_permute <- snakemake@params[["n_permute"]]
seed <- snakemake@params[["seed"]]

cells <- readRDS(input_cells)
joint <- readRDS(input_joint)

# Joint cluster labels (manuscript 1..12), named by cell, and the joint MST edges.
joint_clu <- joint$clusterid
# Hold the PC space identical to the joint tree: the SAME first `pcadim` PCs.
# (Do NOT let Lamian's elbow re-pick the dimension on the subset — on a subset it
# collapses to ~2 PCs and the MST comparison becomes a dimensionality artifact,
# not a condition effect. Variant A = vary the cells only.)
pcadim <- ncol(joint$pca) # 10

canon_edges <- function(mst) {
  el <- igraph::as_edgelist(mst)
  # store as ordered "a|b" so a-b and b-a match; labels are cluster ids.
  sort(apply(el, 1, function(r) paste(sort(as.integer(r)), collapse = "|")))
}
joint_edges <- canon_edges(joint$MSTtree)

# Fixed-cluster cluster-center MST in the joint PC space, exactly as
# TSCAN/Lamian builds it (verified: this reproduces joint$MSTtree on joint cells).
build_edges <- function(prm, clum) {
  mcl <- TSCAN::exprmclust(t(prm), cluster = clum, reduce = FALSE)
  canon_edges(mcl$MSTtree)
}

# ---- subset to condition ----------------------------------------------------
# cond == "joint" runs the SAME fixed-cluster bootstrap on ALL ND+DB cells, to
# give a matched baseline for the per-condition detection rates (the joint tree's
# published detection rates come from Lamian's re-kmeans procedure, which is not
# comparable to this fixed-cluster resampling).
keep <- if (condition == "joint") {
  colnames(cells)
} else {
  colnames(cells)[cells$condition == condition]
}
pr <- cells@reductions$pca@cell.embeddings[keep, seq_len(pcadim), drop = FALSE]
clu <- joint_clu[keep]
stopifnot(!anyNA(clu), all(rownames(pr) == names(clu)))

cellanno <- cells@meta.data[keep, ] %>%
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

# ---- point estimate: fixed-cluster tree in the joint 10-PC space -------------
mcl <- TSCAN::exprmclust(t(pr), cluster = clu, reduce = FALSE)

# origin cluster = the cluster most enriched for D3 (matches infer_tree_structure)
tab <- table(cluster = mcl$clusterid, tp = cellanno$timepoint[match(names(mcl$clusterid), cellanno$cell)])
tab <- tab / rowSums(tab)
origin.cluster <- as.integer(names(which.max(tab[, "D3"])))

ord <- TSCAN::TSCANorder(mcl, listbranch = TRUE, orderonly = TRUE)
pt <- unlist(sapply(sapply(ord, length), function(i) seq_len(i)))
names(pt) <- unname(unlist(ord))
newbranch <- Lamian:::findbranch(mst = mcl$MSTtree, order = ord, origin = origin.cluster)

res <- mcl
res$pseudotime <- pt
res$order <- ord
res$branch <- newbranch
res$pca <- pr
res$allsample <- setNames(cellanno$sample, cellanno$cell)
res$origin.cluster <- origin.cluster
saveRDS(res, out_tree)

pr_used <- pr
clu_used <- clu
cond_edges <- canon_edges(mcl$MSTtree)

# ---- bootstrap edge detection (fixed clusters, no re-kmeans) ----------------
set.seed(seed)
rng <- sample.int(.Machine$integer.max, n_permute)
counts <- integer(0)
for (b in seq_len(n_permute)) {
  set.seed(rng[b])
  idx <- sample(seq_len(nrow(pr_used)), nrow(pr_used), replace = TRUE)
  # a cluster can vanish from a bootstrap draw only if very rare; guard it.
  eb <- tryCatch(build_edges(pr_used[idx, , drop = FALSE], clu_used[idx]),
                 error = function(e) character(0))
  for (e in unique(eb)) counts[e] <- (if (is.na(counts[e])) 0L else counts[e]) + 1L
}
detection <- tibble(edge = names(counts), n = as.integer(counts)) %>%
  mutate(detection_rate = n / n_permute)

# ---- assemble per-edge report (union of joint + condition edges) ------------
all_edges <- union(joint_edges, cond_edges)
edge_tbl <- tibble(edge = all_edges) %>%
  separate(edge, into = c("source", "target"), sep = "\\|", remove = FALSE,
           convert = TRUE) %>%
  mutate(
    condition = condition,
    in_joint_mst = edge %in% joint_edges,
    in_condition_mst = edge %in% cond_edges
  ) %>%
  left_join(select(detection, edge, detection_rate), by = "edge") %>%
  mutate(detection_rate = replace_na(detection_rate, 0)) %>%
  arrange(desc(in_joint_mst), desc(detection_rate)) %>%
  relocate(condition, source, target, in_joint_mst, in_condition_mst,
           detection_rate)

write_csv(edge_tbl, out_edges)

message(sprintf("[%s] joint edges recovered in condition MST: %d / %d",
                condition, sum(edge_tbl$in_joint_mst & edge_tbl$in_condition_mst),
                length(joint_edges)))
message("Joint edges by bootstrap detection rate:")
print(edge_tbl %>% filter(in_joint_mst) %>% select(source, target, in_condition_mst, detection_rate))
