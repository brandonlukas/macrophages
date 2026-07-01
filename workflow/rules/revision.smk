# npj revision — heavy trajectory compute only (T1 / R2-2, editor-flagged).
#
# Deliberately minimal: this repo runs only the compute-intensive steps that
# must go through the pipeline and produce results/ artifacts (synced to Box).
# All lightweight post-processing of these outputs — occupancy tables, branch /
# pseudotime extraction, joint-vs-condition concordance stats, and the figures —
# lives in the moma repo (~/code/timkoh/moma), which reads box/results/revision/.
#
# {cond} in {wt, db}: wt = non-diabetic, db = diabetic. Re-infers the trajectory
# on each condition alone to test whether the joint ND+DB trajectory reproduces.


# Per-condition trajectory. Same params as the published joint run
# (workflow/rules/lamian.smk::infer_tree); only the condition subset differs.
rule revision_infer_tree_condition:
    input:
        config["inputs"]["cells"],
    output:
        "results/revision/lamian/{cond}/infer_tree.rds",
    params:
        seed_1=42,
        seed_2=12345,
    wildcard_constraints:
        cond="wt|db",
    script:
        "../scripts/revision/infer_tree_condition.R"


# Per-condition branch uncertainty (detection.rate), so moma can draw
# per-condition backbones with load_edge_tbl(). Parallel reimplementation of
# Lamian::evaluate_uncertainty (drop-in output).
rule revision_evaluate_uncertainty_condition:
    input:
        rules.revision_infer_tree_condition.output,
    output:
        "results/revision/lamian/{cond}/evaluate_uncertainty.rds",
    params:
        n_permute=1000,
    threads: 24
    wildcard_constraints:
        cond="wt|db",
    script:
        "../scripts/revision/evaluate_uncertainty_parallel.R"
