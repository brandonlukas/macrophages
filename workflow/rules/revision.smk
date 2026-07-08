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
#
# Variant C is reported as an ENSEMBLE of independent draws, not one tree: Lamian's
# clustering is not reproducible (its mykmeans ignores kmeans.seed and forks an
# unseeded cluster-number search — a Lamian bug we document but do not patch). We
# therefore run the inference N times and report robustness across the draws:
#   - per-draw metrics (cluster count, joint pseudotime concordance), and
#   - edge-level cell-Jaccard reproducibility of the joint trajectory's edges
#     (Lamian-style overlap of the cells on an edge), which is label-agnostic and
#     answers "is the APC branch 6->7->5 recovered per condition?".
# See docs/revision-trajectory-robustness.md.

# Draws in the condition-PCA ensemble.
CONDPCA_SEEDS = list(range(1, 11))


# Variant C, one draw. Re-embed HVGs/scaling/PCA on each condition, then native
# Lamian. Seeds are passed exactly as the main pipeline (kmeans.seed is a dead
# Lamian parameter — see the script header); each {seed} is an independent draw.
rule revision_condpca_draw:
    input:
        config["inputs"]["cells"],
    output:
        "results/revision/lamian_condpca_draws/{cond}/seed{seed}/infer_tree.rds",
    params:
        # Vary the RNG seed per draw. set.seed(seed_1) is the seed that actually
        # perturbs the result (via the ambient RNG the fork/kmeans consume);
        # kmeans.seed=seed_2 is passed as the main pipeline does but is a dead
        # Lamian parameter. Distinct seed_1 per draw => 10 independent draws.
        seed_1=lambda wc: int(wc.seed),
        seed_2=12345,
        npcs=50,
    wildcard_constraints:
        cond="wt|db",
        seed=r"\d+",
    script:
        "../scripts/revision/infer_tree_condition_pca.R"


# Ensemble reporting: per-draw manifest (cluster count, joint pseudotime
# concordance) + edge-level reproducibility. For each joint MST edge, restrict the
# joint trajectory to the condition's cells, then measure the best cell-Jaccard /
# overlap-coefficient match among each draw's edges, and a detection rate across
# draws vs a per-edge null. Highlights the APC edges (6-7, 5-7).
rule revision_condpca_edge_detection:
    input:
        joint=rules.infer_tree.output,
        draws=expand(
            "results/revision/lamian_condpca_draws/{cond}/seed{seed}/infer_tree.rds",
            cond=["wt", "db"],
            seed=CONDPCA_SEEDS,
        ),
    output:
        manifest="results/revision/condpca_draw_manifest.csv",
        edges="results/revision/condpca_edge_detection.csv",
    params:
        n_null=1000,
        seed=42,
    script:
        "../scripts/revision/condpca_edge_detection.R"
