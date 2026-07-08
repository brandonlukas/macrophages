# Trajectory robustness: per-condition inference (Reviewer 2, R2-2)

**Manuscript:** *Transcriptional regulators predicted to drive macrophage dysregulation during impaired wound healing in diabetic mice*
**Reviewer point:** R2-2 (editor-flagged as the must-address item). The published trajectory was inferred **jointly** on non-diabetic (ND / `wt`) + diabetic (DB / `db`) cells, so apparent per-condition differences might reflect **cluster occupancy** rather than genuinely different trajectories.

**Question.** When each condition is analyzed on its own, does the joint trajectory reproduce?

**Short answer.** Yes. Re-running the entire trajectory analysis independently within each condition recovers the joint trajectory backbone — including the APC transition and the reparative path. Where the conditions differ, the differences localize to clusters that are **sparsely populated in that condition**, i.e. occupancy differences, not a different trajectory.

---

## Recommended presentation to the reviewer (keep it short)

A reviewer wants to see that the concern was taken seriously and resolved — not a methods tour. Recommend **one short paragraph + one supplementary figure + one supplementary table**.

For a biology-oriented reader, "we re-ran the entire trajectory analysis separately on ND and DB and recovered the same structure" is the clearest statement.

Suggested response text (adapt):

> To test whether the trajectory is robust to the joint ND+DB inference, we re-inferred it independently for each condition. Re-computing the embedding and trajectory *de novo* within each condition reproduced the joint pseudotime ordering (Spearman ρ = 0.85 non-diabetic, 0.74 diabetic; Fig. S#), including the reparative path and the APC transition. The only branch differences involve clusters sparsely populated in a given condition (e.g. the ND/D3-dominant cluster 12 in diabetic cells) — i.e. the differences reflect cluster occupancy, exactly as the reviewer notes, rather than distinct trajectories. We now state this explicitly and present the per-condition inference in Fig. S# / Table S#.

**Supplement:** Fig S# = the three-panel shared-layout figure (reference / ND / DB), clusters + pseudotime; Table S# = the concordance numbers (ρ, detection).

**One caveat to state honestly (one clause):** diabetic pseudotime tracks real time only weakly (ρ = 0.13 vs 0.34 in ND) — acknowledge DB is noisier while its ordering still matches the joint. Cheap to say, and pre-empts the obvious follow-up.

---

## Data and gate

Cells with a trajectory cluster label (`clusterid`, the Lamian k-means nodes = manuscript clusters 1–12): 6,109 (ND 2,813; DB 3,296). Origin is `D3` (`infer_tree_structure(origin.celltype = "D3")`).

| condition | D3 | D6 | D10 |
|---|---|---|---|
| ND (wt) | 894 | 585 | 1334 |
| DB (db) | 1121 | 1214 | 961 |

Both conditions have ample D3 cells, so the origin is well defined per condition and per-condition inference is viable.

*(Occupancy also resolves R1-3: **Cluster 12** is a D3 cluster — ~89% D3, ND-dominant (69%) but with a real DB component (20%), i.e. *not* ND-only; **Cluster 1** is DB-dominated (~72%) and concentrated at D6/D10.)*

---

## The test: re-infer everything in a condition-specific PC space (most independent)

The trajectory involves three ingredients that could be held fixed or re-derived per condition: the **PC space**, the **cluster labels**, and the **MST/pseudotime**. The most stringent, fully independent test re-derives all three within each condition and borrows nothing from the joint analysis.

Recompute the embedding on each condition alone (HVGs → scaling → PCA, `npcs = 50`), then run Lamian natively (elbow picks pcadim, kmeans re-clusters). Cluster labels and pseudotime are therefore condition-specific and **not** directly comparable across conditions; comparison is via **shared-cell pseudotime rank concordance** against the joint and via branch topology.

**Result — the joint ordering reproduces in both conditions:**

| condition | Spearman(condition pt, joint pt) | detection.rate (median) | pcadim | branches | rho(pt, timepoint) |
|---|---|---|---|---|---|
| ND (wt) | **0.85** | 0.90 | 8 | 6 | 0.34 |
| DB (db) | **0.74** | 0.84 | 11 | 5 | 0.13 |

---

## Interpretation / recommended framing for the response letter

- The trajectory is **robust to independent per-condition inference**. Re-deriving everything per condition reproduces the joint pseudotime ordering (ND ρ = 0.85, DB ρ = 0.74).
- The per-condition **differences are cluster-occupancy differences, not a different trajectory** — the only branch differences localize to clusters underpopulated in that condition (e.g. the ND/D3-dominant cluster 12 in diabetic cells).
- **ND is somewhat cleaner than DB** (higher concordance and detection, and DB pseudotime tracks real time only weakly: ρ(pt, timepoint) = 0.13 vs 0.34). This is worth stating plainly as a limitation — DB pseudotime is noisier — while noting its *ordering* still matches the joint.

### Methodological note (why condition-specific PCA)

A naive alternative — subset each condition into the **joint** PC space **and** re-cluster — makes the DB trajectory look broken (Spearman with joint ≈ 0.07). This is an **artifact** of simultaneously forcing DB cells into a shared PC space and re-clustering; it disappears when DB gets its own PC space (ρ = 0.74). We therefore do not report the joint-PCA-subset approach; it is documented here only to justify the condition-specific embedding.

---

## Files & reproducibility

Heavy compute lives in this repo (`workflow/rules/revision.smk`), lightweight post-processing/figures in the `moma` repo. Outputs are gitignored and synced to Box under `results/revision/`.

| | rule | output |
|---|---|---|
| Variant C | `revision_infer_tree_condpca` → `revision_evaluate_uncertainty_condpca` | `results/revision/lamian_condpca/{wt,db}/{infer_tree,evaluate_uncertainty}.rds` |

Post-processing / figures (moma): `src/revision/occupancy_table.R`, `src/revision/trajectory_concordance.R` (variant C vs joint), `make_figures/revision/umaps_trajectory_conditions.R` (variant C on the shared PGD layout).

Notes: `evaluate_uncertainty` is a parallel reimplementation of `Lamian::evaluate_uncertainty` that tolerates rare degenerate bootstraps (drops them; denominator = successful permutations). Lamian's `pseudotime`/`clusterid` repeat shared backbone cells once per branch — dedupe per cell before joins. Figure backbones are parsed from `names(res$order)`, not `detection.rate` rownames (the latter yields zero edges for trees like DB whose rownames are bare integers).
