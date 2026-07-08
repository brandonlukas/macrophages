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
