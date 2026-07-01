# Revision task T1 — per-condition branch uncertainty (detection.rate), needed
# so the moma figure can draw per-condition backbones with the house-style
# load_edge_tbl(). Parallel reimplementation of Lamian::evaluate_uncertainty:
# per-permutation work is independent (only post-loop aggregation is serial), and
# per-iter RNG seeds are pre-drawn so results are reproducible regardless of
# worker scheduling. Output list matches Lamian's (detection.rate,
# sample.cellcomp.mean/sd), so it is a drop-in for the joint evaluate_uncertainty.

library(Lamian)
library(parallel)
library(reshape2)

evaluate_uncertainty_parallel <- function(inferobj, n.permute, mc.cores = 1L,
                                          subset.cell = NULL, design = NULL,
                                          return.ctcomp = FALSE, seed = 12345) {
  pr <- if (is.null(subset.cell)) inferobj$pca else inferobj$pca[subset.cell, ]
  newbranch <- inferobj$branch
  js.cut <- inferobj$js.cut
  oc.cut <- inferobj$oc.cut
  pt <- inferobj$pseudotime
  ord <- inferobj$order
  alls <- inferobj$allsample
  nclust <- max(inferobj$clusterid)
  set.seed(seed)
  rng.seeds <- sample.int(.Machine$integer.max, n.permute)

  one_perm <- function(pmid) {
    set.seed(rng.seeds[pmid])
    bstid <- unique(sample(seq_len(nrow(pr)), nrow(pr), replace = TRUE))
    pr.pm <- pr[bstid, ]
    invisible(capture.output(
      clu <- Lamian:::mykmeans(pr.pm, number.cluster = nclust)$cluster
    ))
    mcl.pm <- TSCAN::exprmclust(t(pr.pm), cluster = clu, reduce = FALSE)
    pt.pm.mean <- tapply(pt[names(mcl.pm[["clusterid"]])],
                         list(mcl.pm[["clusterid"]]), mean)
    start.cluster <- names(which.min(pt.pm.mean))
    ord.pm <- TSCAN::TSCANorder(mcl.pm, startcluster = start.cluster,
                                listbranch = TRUE, orderonly = TRUE)
    pt.pm <- unlist(sapply(sapply(ord.pm, length), function(i) seq(1, i)))
    names(pt.pm) <- unname(unlist(ord.pm))
    newbranch.pm <- Lamian:::findbranch(mst = mcl.pm$MSTtree, order = ord.pm,
                                        origin = start.cluster)
    js <- sapply(seq_along(newbranch), function(i) {
      id <- which(sapply(paste0(names(ord), ","), function(k)
        grepl(paste0(paste0(newbranch[[i]], collapse = ","), ","), k)))[1]
      cells <- ord[[id]]
      b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k)
        names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
      sapply(seq_along(newbranch.pm), function(j) {
        id <- which(sapply(paste0(names(ord.pm), ","), function(k)
          grepl(paste0(paste0(newbranch.pm[[j]], collapse = ","), ","), k)))[1]
        cells <- ord.pm[[id]]
        b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k)
          names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
        length(intersect(b.pm, b.ori)) / length(union(b.pm, b.ori))
      })
    })
    oc <- sapply(seq_along(newbranch), function(i) {
      id <- which(sapply(paste0(names(ord), ","), function(k)
        grepl(paste0(paste0(newbranch[[i]], collapse = ","), ","), k)))[1]
      cells <- ord[[id]]
      b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k)
        names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
      sapply(seq_along(newbranch.pm), function(j) {
        id <- which(sapply(paste0(names(ord.pm), ","), function(k)
          grepl(paste0(paste0(newbranch.pm[[j]], collapse = ","), ","), k)))[1]
        cells <- ord.pm[[id]]
        b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k)
          names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
        length(intersect(b.pm, b.ori)) / min(length(b.pm), length(b.ori))
      })
    })
    corr <- sapply(seq_along(newbranch), function(i) {
      id <- which(sapply(paste0(names(ord), ","), function(k)
        grepl(paste0(paste0(newbranch[[i]], collapse = ","), ","), k)))[1]
      cells <- ord[[id]]
      b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k)
        names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
      sapply(seq_along(newbranch.pm), function(j) {
        id <- which(sapply(paste0(names(ord.pm), ","), function(k)
          grepl(paste0(paste0(newbranch.pm[[j]], collapse = ","), ","), k)))[1]
        cells <- ord.pm[[id]]
        b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k)
          names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
        ov <- intersect(b.ori, b.pm)
        cor(pt[ov], pt.pm[ov])
      })
    })
    corr[is.na(corr)] <- 0
    colnames(corr) <- colnames(oc) <- colnames(js) <-
      paste0("original", seq_along(newbranch))
    js.binary <- Lamian:::get_binary(js, js.cut)
    corr.score <- corr * js.binary
    js.melt <- melt(js.binary)
    js.melt <- js.melt[js.melt[, 3] != 0, ]
    colnames(js.melt) <- c("permutation.branch", "original.branch", "matched")
    reproduce.js <- as.character(js.melt[, 2])
    oc.binary <- Lamian:::get_binary(oc, oc.cut)
    oc.melt <- melt(oc.binary)
    oc.melt <- oc.melt[oc.melt[, 3] != 0, ]
    reproduce.oc <- as.character(oc.melt[, 2])
    ctcomp.new.logit <- ctcomp.new <- matrix(0,
      nrow = length(unique(alls)), ncol = length(newbranch))
    colnames(ctcomp.new.logit) <- colnames(ctcomp.new) <-
      paste0("origin", seq_along(newbranch))
    rownames(ctcomp.new.logit) <- rownames(ctcomp.new) <- unique(alls)
    if (nrow(js.melt) > 0) {
      ctcomp <- sapply(seq_len(nrow(js.melt)), function(i) {
        c <- names(clu)[clu %in% newbranch.pm[[js.melt[i, 1]]]]
        ctcomp <- rep(0, length(unique(alls)))
        names(ctcomp) <- unique(alls)
        ctcomp[names(table(alls[c]))] <- table(alls[c])
        ctcomp
      })
      colnames(ctcomp) <- paste0("origin", js.melt[, 2])
      ctcomp.logit.tmp <- (ctcomp + 1) / (rowSums(ctcomp) + 1)
      ctcomp.logit <- log(ctcomp.logit.tmp / (1 - ctcomp.logit.tmp))
      ctcomp.new.logit[rownames(ctcomp.logit), colnames(ctcomp.logit)] <-
        ctcomp.logit
      ctcomp <- ctcomp / rowSums(ctcomp)
      ctcomp.new[rownames(ctcomp), colnames(ctcomp)] <- ctcomp
    }
    list(
      corr.score = corr.score,
      reproduce.js = reproduce.js,
      reproduce.oc = reproduce.oc,
      ctcomp = t(ctcomp.new),
      ctcomp.logit = t(ctcomp.new.logit)
    )
  }

  results <- mclapply(seq_len(n.permute), one_perm,
                      mc.cores = mc.cores, mc.preschedule = FALSE)
  fail <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(fail)) {
    first <- results[[which(fail)[1]]]
    stop(sprintf("%d permutation(s) failed; first error: %s",
                 sum(fail), attr(first, "condition")$message))
  }
  corr.score   <- lapply(results, `[[`, "corr.score")
  reproduce.js <- unlist(lapply(results, `[[`, "reproduce.js"))
  reproduce.oc <- unlist(lapply(results, `[[`, "reproduce.oc"))
  ctcomplist       <- lapply(results, `[[`, "ctcomp")
  ctcomplist.logit <- lapply(results, `[[`, "ctcomp.logit")
  js.perc <- rep(0, length(newbranch))
  if (length(reproduce.js))
    js.perc[as.numeric(names(table(reproduce.js)))] <-
      table(reproduce.js) / n.permute
  names(js.perc) <- newbranch
  oc.perc <- rep(0, length(newbranch))
  if (length(reproduce.oc))
    oc.perc[as.numeric(names(table(reproduce.oc)))] <-
      table(reproduce.oc) / n.permute
  names(oc.perc) <- newbranch
  corr.score.m <- do.call(rbind, corr.score)
  corr.score.v <- colSums(corr.score.m) / n.permute
  names(corr.score.v) <- newbranch
  detection.rate <- data.frame(
    detection.rate = (js.perc + oc.perc[names(js.perc)]) / 2,
    stringsAsFactors = FALSE)
  detection.rate[detection.rate[, 1] > 1, 1] <- 1
  sample.cellcomp.mean <- apply(simplify2array(ctcomplist), seq_len(2), mean)
  sample.cellcomp.sd   <- apply(simplify2array(ctcomplist), seq_len(2), sd)
  rownames(sample.cellcomp.mean) <-
    newbranch[as.numeric(sub("origin", "", rownames(sample.cellcomp.mean)))]
  rownames(sample.cellcomp.sd) <-
    newbranch[as.numeric(sub("origin", "", rownames(sample.cellcomp.sd)))]
  if (length(newbranch[[length(newbranch)]]) == 2) {
    nm <- paste0("c(", newbranch[[length(newbranch)]][1], ",",
                 newbranch[[length(newbranch)]][2], ")")
    rownames(detection.rate)[nrow(detection.rate)] <-
      rownames(sample.cellcomp.mean)[nrow(sample.cellcomp.mean)] <-
      rownames(sample.cellcomp.sd)[nrow(sample.cellcomp.sd)] <- nm
  }
  result <- list(detection.rate = detection.rate,
                 sample.cellcomp.mean = sample.cellcomp.mean,
                 sample.cellcomp.sd = sample.cellcomp.sd)
  if (return.ctcomp) {
    result$branchProp <- ctcomplist
    result$branchProp.logit <- ctcomplist.logit
  }
  result
}

# ---- snakemake driver -------------------------------------------------------
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
n_permute <- snakemake@params[["n_permute"]]
mc_cores <- min(snakemake@threads, n_permute)

res <- readRDS(input_file)

cat(sprintf("[%s] evaluate_uncertainty n.permute=%d mc.cores=%d\n",
            format(Sys.time()), n_permute, mc_cores))
t0 <- Sys.time()
tree <- evaluate_uncertainty_parallel(
  res, n.permute = n_permute, mc.cores = mc_cores, seed = 12345
)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("[%s] finished in %.1f s (%.1f min)\n",
            format(Sys.time()), elapsed, elapsed / 60))

saveRDS(tree, output_file)
