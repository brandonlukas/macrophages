# npj revision — heavy trajectory compute only (T1 / R2-2, editor-flagged).
#
# Deliberately minimal: this repo runs only the compute-intensive steps that
# must go through the pipeline and produce results/ artifacts (synced to Box).
# All lightweight post-processing of these outputs — occupancy tables, branch /
# pseudotime extraction, joint-vs-condition concordance stats, and the figures —
# lives in the moma repo (~/code/timkoh/moma), which reads box/results/revision/.
#
# {cond} in {wt, db}: wt = non-diabetic, db = diabetic. Two complementary tests of
# whether the joint ND+DB trajectory reproduces per condition:
#   A = hold the joint clusters + joint PC space FIXED (collaborator request),
#   C = re-infer everything in a condition-specific PC space (most independent).
# See docs/revision-trajectory-robustness.md.


# Variant A (PI/collaborator request): per-condition trajectory holding the JOINT
# clusters FIXED, then bootstrapping the cluster-center MST to get per-edge
# detection rates (does each condition reconstruct the joint edges, and how
# reliably). Directly comparable across joint/wt/db because cluster labels match.
rule revision_infer_tree_fixedclusters:
    input:
        cells=config["inputs"]["cells"],
        joint=rules.infer_tree.output,
    output:
        tree="results/revision/lamian_fixedclusters/{cond}/tree.rds",
        edges="results/revision/lamian_fixedclusters/{cond}/edge_detection.csv",
    params:
        n_permute=1000,
        seed=42,
    wildcard_constraints:
        cond="wt|db",
    script:
        "../scripts/revision/infer_tree_fixedclusters.R"


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
