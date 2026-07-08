# Reproducible replacement for Lamian:::mykmeans (revision T1 / R2-2).
#
# The shipped Lamian::mykmeans is nondeterministic, for two reasons:
#   1. Its `seed` argument is never used — there is no set.seed() in the body,
#      so `infer_tree_structure(kmeans.seed=...)` has no effect.
#   2. The cluster-number selection runs inside an unseeded mclapply(mc.cores=30)
#      fork, and every kmeans() uses the default single random start.
# Consequence: repeated runs on identical input / PCA / seeds return different
# cluster counts and trajectories (observed: 12 vs 14 clusters, and joint
# concordance ranging 0.16-0.74, for the diabetic condition).
#
# This override (a) honors the seed, (b) uses `nstart` restarts so kmeans reaches
# a near-global optimum, and (c) drops the fork (lapply). It makes the inference
# reproducible; it does NOT make the diabetic trajectory identifiable — DB has
# several near-equal-WSS clusterings with very different topology (see
# docs/revision-trajectory-robustness.md and the condpca stability sweep).
#
# Install by overwriting the function in Lamian's namespace so the internal
# infer_tree_structure() call picks it up:
#   source("workflow/scripts/revision/reproducible_kmeans.R")
#   install_reproducible_kmeans(nstart)
install_reproducible_kmeans <- function(nstart = 100L) {
  mykmeans_det <- function(matrix, number.cluster = NA, maxclunum = 30, seed = 12345) {
    set.seed(seed)
    if (is.na(number.cluster)) {
      rss <- unlist(lapply(seq_len(maxclunum), function(clunum) {
        km <- kmeans(matrix, clunum, iter.max = 1000, nstart = nstart)
        km$betweenss / km$totss
      }))
      x <- 2:maxclunum
      number.cluster <- x[which.min(sapply(seq_len(length(x)), function(i) {
        x2 <- pmax(0, x - i)
        sum(lm(rss[-1] ~ x + x2)$residuals^2)
      }))]
    }
    set.seed(seed) # re-seed so the final fit is independent of the elbow RNG path
    kmeans(matrix, number.cluster, iter.max = 1000, nstart = nstart)
  }
  assignInNamespace("mykmeans", mykmeans_det, ns = "Lamian")
  invisible(mykmeans_det)
}
