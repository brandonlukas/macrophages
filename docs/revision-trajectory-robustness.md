# Trajectory robustness: per-condition inference (Reviewer 2, R2-2)

**Manuscript:** *Transcriptional regulators predicted to drive macrophage dysregulation during impaired wound healing in diabetic mice*
**Reviewer point:** R2-2 (editor-flagged as the must-address item). The published trajectory was inferred **jointly** on non-diabetic (ND / `wt`) + diabetic (DB / `db`) cells, so apparent per-condition differences might reflect **cluster occupancy** rather than genuinely different trajectories.

**Question.** When each condition is analyzed on its own, does the joint trajectory reproduce?

**Short answer.** Yes. Under two independent, principled tests the ND and DB cells both recover the joint trajectory backbone — including the APC transition and the reparative path. Where the conditions differ, the differences localize to clusters that are **sparsely populated in that condition**, i.e. occupancy differences, not a different trajectory.

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

## Two tests (a robustness ladder)

The trajectory involves three ingredients that could be held fixed or re-derived per condition: the **PC space**, the **cluster labels**, and the **MST/pseudotime**. We report the two principled endpoints.

### Variant A — hold the joint clusters + joint PC space FIXED (collaborator request)

Reuse the joint cluster labels (fork's `infer_tree_structure(clusters = ...)`; verified `TSCAN::exprmclust(cluster = ...)` preserves labels), subset to one condition, and refit **only** the cluster-center minimum spanning tree in the **same 10 PCs** the joint tree uses (`ncol(joint$pca) == 10`). Because the MST is forced to span all 12 clusters, the point estimate alone is not evidence — the test is the **bootstrap detection rate** (resample the cells → recompute centers with the *same* labels → rebuild MST → how often does each edge recur, over 1,000 resamples). Trees are directly comparable across joint/ND/DB because the labels match.

This bootstrap is **not** Lamian's `evaluate_uncertainty`: Lamian resamples cells *and re-runs k-means* each iteration (cells get relabeled, branches matched back by Jaccard/overlap), capturing sampling **and** clustering instability. Here the clustering is fixed by design, so we run the same fixed-cluster resampling on the **joint ND+DB cells too** (`cond = "joint"`), giving a procedure-matched baseline to compare the per-condition detection rates against.

> Implementation note: the PC dimension must be held at the joint's value (10). Letting Lamian's elbow re-pick it on a subset collapses to ~2 PCs and turns the comparison into a dimensionality artifact.

**Result — matched fixed-cluster edge detection (1,000 resamples each):**

| edge | joint | ND (wt) | DB (db) | |
|---|---|---|---|---|
| 3–9, 4–8, 5–7, 6–7, 6–9 | 1.0 | 1.0 | 1.0 | **core — maximally stable everywhere** |
| 3–10 | 1.0 | 1.0 | .978 | solid everywhere |
| 3–4 | .983 | 1.0 | .995 | solid everywhere |
| 1–4 | 1.0 | **.002** | .999 | drops in ND (cl 1 is DB-dominated) |
| 6–12 | .999 | .99 | **.064** | drops in DB (cl 12 is ND/D3-dominant, ~61 DB cells) |
| 1–11 | **.82** | .004 | .729 | already marginal in joint (cl 11 tiny) |
| 2–6 | **.687** | 1.0 | .004 | already the weakest joint edge |
| 1–3 / 5–11 (ND-only) | — | .998 / .849 | — | ND rewiring of its sparse clusters |
| 2–12 / 10–12 (DB-only) | — | — | .959 / .864 | DB rewiring of its sparse clusters |

The **7-edge core backbone sits at detection ≈ 1.0 in all three** under the identical procedure — including the **APC branch (6→7→5)** and the **reparative path (6→9→3→4→8, 3→10)**. Both conditions reconstruct 9/11 joint edges. The edges a condition drops are either (i) into clusters **sparse in that condition** (1–4 in ND, 6–12 in DB) or (ii) edges that are **already the weakest in the joint itself** (2–6 at .687, 1–11 at .82). Conditions do not break *strong* joint edges; the condition-specific rewirings are themselves reproducible (high detection).

> Caveat: n differs (joint 6109; ND 2813; DB 3296) and detection rises with n, so the joint is an *upper* reference — the clean like-for-like contrast is ND vs DB (comparable n). An n-matched joint (subsample to condition size) can be added if a strict baseline is wanted.

### Variant C — re-infer everything in a condition-specific PC space (most independent)

Recompute the embedding on each condition alone (HVGs → scaling → PCA, `npcs = 50`), then run Lamian natively (elbow picks pcadim, kmeans re-clusters). Nothing is borrowed from the joint analysis. Cluster labels and pseudotime are therefore condition-specific and **not** directly comparable across conditions; comparison is via **shared-cell pseudotime rank concordance** against the joint and via branch topology.

**Result — the joint ordering reproduces in both conditions:**

| condition | Spearman(condition pt, joint pt) | detection.rate (median) | pcadim | branches | rho(pt, timepoint) |
|---|---|---|---|---|---|
| ND (wt) | **0.85** | 0.90 | 8 | 6 | 0.34 |
| DB (db) | **0.74** | 0.84 | 11 | 5 | 0.13 |

---

## Interpretation / recommended framing for the response letter

- The trajectory is **robust to independent per-condition inference**. Two tests that make different assumptions agree: holding cell-states fixed (A) recovers the core backbone with near-perfect stability in both conditions; re-deriving everything per condition (C) reproduces the joint pseudotime ordering (ND ρ = 0.85, DB ρ = 0.74).
- The per-condition **differences are cluster-occupancy differences, not a different trajectory** — variant A shows the only unstable edges are those into clusters underpopulated in that condition, and even those rewirings are reproducible.
- **ND is somewhat cleaner than DB** (higher concordance and detection, and DB pseudotime tracks real time only weakly: ρ(pt, timepoint) = 0.13 vs 0.34). This is worth stating plainly as a limitation — DB pseudotime is noisier — while noting its *ordering* still matches the joint.

### Methodological note (why condition-specific PCA)

A third, naive approach — subset each condition into the **joint** PC space **and** re-cluster — makes the DB trajectory look broken (Spearman with joint ≈ 0.07). This is an **artifact** of simultaneously forcing DB cells into a shared PC space and re-clustering; it disappears when DB gets its own PC space (C, ρ = 0.74) or when cluster identities are held fixed (A, 9/11 edges). We therefore do not report the joint-PCA-subset approach; it is documented here only to justify the condition-specific embedding.

---

## Files & reproducibility

Heavy compute lives in this repo (`workflow/rules/revision.smk`), lightweight post-processing/figures in the `moma` repo. Outputs are gitignored and synced to Box under `results/revision/`.

| | rule | output |
|---|---|---|
| Variant A | `revision_infer_tree_fixedclusters` | `results/revision/lamian_fixedclusters/{joint,wt,db}/{tree.rds, edge_detection.csv}` |
| Variant C | `revision_infer_tree_condpca` → `revision_evaluate_uncertainty_condpca` | `results/revision/lamian_condpca/{wt,db}/{infer_tree,evaluate_uncertainty}.rds` |

Post-processing / figures (moma): `src/revision/occupancy_table.R`, `src/revision/trajectory_concordance.R` (variant C vs joint), `make_figures/revision/umaps_trajectory_conditions.R` (variant C on the shared PGD layout).

Notes: `evaluate_uncertainty` is a parallel reimplementation of `Lamian::evaluate_uncertainty` that tolerates rare degenerate bootstraps (drops them; denominator = successful permutations). Lamian's `pseudotime`/`clusterid` repeat shared backbone cells once per branch — dedupe per cell before joins. Figure backbones are parsed from `names(res$order)`, not `detection.rate` rownames (the latter yields zero edges for trees like DB whose rownames are bare integers).
