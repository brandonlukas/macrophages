# npj revision — heavy trajectory compute only (T1 / R2-2, editor-flagged).
#
# Deliberately minimal: this repo runs only the compute-intensive steps that
# must go through the pipeline and produce results/ artifacts (synced to Box).
# All lightweight post-processing of these outputs — occupancy tables, branch /
# pseudotime extraction, joint-vs-condition concordance stats, and the figures —
# lives in the moma repo (~/code/timkoh/moma), which reads box/results/revision/.
#
# {cond} in {wt, db}: wt = non-diabetic, db = diabetic. Test of whether the joint
# ND+DB trajectory reproduces per condition by re-inferring everything in a
# condition-specific PC space (the "most independent" test).
# See docs/revision-trajectory-robustness.md.

# Seeds for the condition-PCA ensemble (below). A single seed is reproducible but
# can be a lucky/unlucky draw (esp. DB, which is non-identifiable), so we save the
# full tree for each and look at the consensus.
CONDPCA_SEEDS = list(range(1, 11))


# Variant C ("most independent"): per-condition trajectory in a CONDITION-SPECIFIC
# PC space (re-embed HVGs/scaling/PCA on each condition, then native Lamian).
rule revision_infer_tree_condpca:
    input:
        config["inputs"]["cells"],
    output:
        "results/revision/lamian_condpca/{cond}/infer_tree.rds",
    params:
        seed_1=42,
        seed_2=12345,
        npcs=50,
        # nstart: kmeans restarts. Lamian's mykmeans is patched to be reproducible
        # (it otherwise ignores its seed and forks an unseeded cluster search);
        # restarts reach a near-global optimum. See reproducible_kmeans.R.
        nstart=100,
    wildcard_constraints:
        cond="wt|db",
    script:
        "../scripts/revision/infer_tree_condition_pca.R"


rule revision_evaluate_uncertainty_condpca:
    input:
        rules.revision_infer_tree_condpca.output,
    output:
        "results/revision/lamian_condpca/{cond}/evaluate_uncertainty.rds",
    params:
        n_permute=1000,
    threads: 24
    wildcard_constraints:
        cond="wt|db",
    script:
        "../scripts/revision/evaluate_uncertainty_parallel.R"


# Stability sweep: quantifies how much the per-condition trajectory depends on
# which near-optimal k-means solution is found (many seeds, fixed nstart). ND is
# stable (rho ~= 0.76); DB is non-identifiable (rho ranges ~0.16-0.74). This is
# the honest robustness characterization behind the R2-2 response.
rule revision_condpca_stability:
    input:
        cells=config["inputs"]["cells"],
        joint=rules.infer_tree.output,
    output:
        "results/revision/condpca_stability.csv",
    params:
        seed_1=42,
        npcs=50,
        nstart=100,
        n_seed=20,
    script:
        "../scripts/revision/condpca_stability.R"


# Ensemble: save the FULL per-condition tree (with pca_full) for each of a handful
# of k-means seeds, so the consensus across seeds can be inspected instead of
# trusting one draw. Reuses the canonical infer script; the only difference is the
# k-means seed comes from the {seed} wildcard and the output path is per-seed.
rule revision_condpca_seed:
    input:
        config["inputs"]["cells"],
    output:
        "results/revision/lamian_condpca_seeds/{cond}/seed{seed}/infer_tree.rds",
    params:
        seed_1=42,
        seed_2=lambda wc: int(wc.seed),
        npcs=50,
        nstart=100,
    wildcard_constraints:
        cond="wt|db",
        seed=r"\d+",
    script:
        "../scripts/revision/infer_tree_condition_pca.R"


# Consensus across the ensemble: per-seed manifest (k, rho vs joint, WSS), how much
# the seeds agree (pairwise adjusted Rand index), a representative seed (most
# central) and the best-concordance seed, plus a co-association consensus
# clustering. This is the "don't trust one seed" deliverable.
rule revision_condpca_consensus:
    input:
        joint=rules.infer_tree.output,
        seeds=expand(
            "results/revision/lamian_condpca_seeds/{cond}/seed{seed}/infer_tree.rds",
            cond=["wt", "db"],
            seed=CONDPCA_SEEDS,
        ),
    output:
        manifest="results/revision/condpca_consensus_manifest.csv",
        summary="results/revision/condpca_consensus_summary.csv",
        labels="results/revision/condpca_consensus_labels.csv",
    script:
        "../scripts/revision/condpca_consensus.R"
