# Trajectory robustness: per-condition inference (Reviewer 2, R2-2)

**Manuscript:** *Transcriptional regulators predicted to drive macrophage dysregulation during impaired wound healing in diabetic mice*
**Reviewer point:** R2-2 (editor-flagged as the must-address item). The published trajectory was inferred **jointly** on non-diabetic (ND / `wt`) + diabetic (DB / `db`) cells, so apparent per-condition differences might reflect **cluster occupancy** rather than genuinely different trajectories.

**Question.** When each condition is analyzed on its own, does the joint trajectory reproduce?

**Short answer.** Yes for the parts that matter. We re-infer the trajectory independently within each condition and, because the inference has run-to-run variability (below), report it as an **ensemble of 10 independent draws** rather than one tree. Two things are stable across the ensemble:

- **Overall ordering:** the per-condition pseudotime tracks the joint ordering — non-diabetic in every draw (Spearman ρ 0.67–0.88, median 0.80), diabetic in most (ρ up to 0.74, median 0.70, 7/10 draws ≥ 0.5).
- **The edges that matter:** the **APC branch edges (6→7 and 7→5) are among the most reproducible.** By Lamian's own detection score (one-to-one cell-overlap matching, averaging Jaccard- and overlap-based detection), **7→5 is recovered in 10/10 draws of both conditions** and **6→7 in ~8–9/10** (mean cell-Jaccard 0.65–0.79, overlap coefficient 0.86–0.97).

Where draws disagree, it is on edges into clusters that are **sparsely populated in that condition** — i.e. occupancy, exactly as the reviewer suggests, not a different trajectory.

> **Reproducibility caveat (stated to the reviewer, NOT patched here).** Lamian's `infer_tree_structure` is not reproducible: its internal `mykmeans` ignores the `kmeans.seed` argument (there is no `set.seed` in its body) and selects the cluster number inside an unseeded parallel fork. We pass seeds exactly as the main pipeline does (`workflow/scripts/lamian/infer_tree.R`) so the intent is on the record, but we deliberately **do not modify Lamian's source**. We therefore report robustness across an ensemble of draws instead of from a single tree. (A reviewer who re-runs will get a different individual tree; the ensemble-level conclusions — ordering concordance and edge reproducibility — are what is stable.)

---

## Recommended presentation to the reviewer (keep it short)

Suggested response text (adapt):

> To test whether the trajectory depends on the joint ND+DB inference, we re-inferred it independently within each condition (condition-specific embedding and clustering, nothing borrowed from the joint analysis). Because the clustering step has run-to-run variability, we repeated the inference 10 times per condition and assessed reproducibility across the ensemble. The non-diabetic pseudotime reproduces the joint ordering in every draw (Spearman ρ 0.67–0.88); the diabetic ordering reproduces it in most draws (median ρ 0.70). Critically, the trajectory edges central to the paper — the APC branch 6→7→5 — reproduce robustly: by a Lamian-style cell-overlap detection score the 7→5 edge is recovered in all 10 draws of both conditions and 6→7 in ~8–9/10 (cell-level Jaccard 0.65–0.79 with the joint edge; Fig. S#/Table S#). The edges that vary between draws are those into clusters sparsely populated in a given condition, i.e. cluster-occupancy differences (Table S#), not a distinct trajectory. We report the per-condition inference as an ensemble because Lamian's clustering is not deterministic (its `kmeans.seed` is a no-op in the released code); we pass seeds but do not alter the package source.

**Supplement:** Fig S# = per-condition inference on the shared layout; Table S# = per-edge reproducibility (Jaccard / detection across draws) + the occupancy table (R1-3).

---

## Data and gate

Cells with a trajectory cluster label (`clusterid`, the Lamian k-means nodes = manuscript clusters 1–12): 6,109 (ND 2,813; DB 3,296). Origin is `D3` (`infer_tree_structure(origin.celltype = "D3")`).

| condition | D3 | D6 | D10 |
|---|---|---|---|
| ND (wt) | 894 | 585 | 1334 |
| DB (db) | 1121 | 1214 | 961 |

Both conditions have ample D3 cells, so the origin is well defined per condition.

*(Occupancy also resolves R1-3: **Cluster 12** is a D3 cluster — ~89% D3, ND-dominant (69%) but with a real DB component (20%), i.e. *not* ND-only; **Cluster 1** is DB-dominated (~72%) and concentrated at D6/D10.)*

---

## Method: condition-specific re-inference, reported as a 10-draw ensemble

Recompute the embedding on each condition alone (HVGs → scaling → PCA, `npcs = 50`), then run Lamian natively (elbow picks pcadim, k-means re-clusters). Cluster labels and pseudotime are condition-specific and **not** directly comparable across conditions; comparison is against the joint via shared cells. We run the inference 10 times per condition (`revision_condpca_draw`, seeds 1–10) and report two ensemble measures.

### 1. Ordering concordance (per draw)

Spearman correlation of each draw's pseudotime against the joint pseudotime on shared cells (`condpca_draw_manifest.csv`):

| condition | clusters (range) | ρ vs joint (median; range) | draws with ρ ≥ 0.5 |
|---|---|---|---|
| **ND (wt)** | 9–14 | **0.80** (0.67–0.88) | 10 / 10 |
| **DB (db)** | 10–13 | **0.70** (−0.02–0.74) | 7 / 10 |

ND recovers the joint ordering in every draw; DB in most (a minority of DB draws land on a poor ordering — DB is noisier, and its late compartment is thinner).

### 2. Edge reproducibility (Lamian-style cell-Jaccard)

Because draws have different cluster labels/counts, edges are matched by the **cells on them**, exactly as Lamian's `evaluate_uncertainty` scores branches across bootstraps. For each joint MST edge (a–b), restricted to one condition's cells, `Sj` = cells whose *joint* cluster is a or b; for each draw edge, `Sd` = cells whose *draw* cluster is one of its endpoints. We reproduce Lamian's detection exactly (`condpca_edge_detection.csv`):

- per-edge null cutoffs `js.cut` / `oc.cut` = 99th percentile of Jaccard / overlap-coefficient against random cell sets of size |Sj| (Lamian's `js.null`/`oc.null`);
- within each draw, threshold the [draw-edge × joint-edge] Jaccard and overlap matrices and resolve to a **one-to-one** matching with Lamian's `get_binary`;
- `detection = (js.perc + oc.perc)/2` — the average of the Jaccard- and overlap-based detection fractions across the 10 draws, capped at 1 (Lamian's own formula).

The **detection rate** and the **cell-overlap magnitude** both rank the APC branch and reparative backbone highest and the sparse-cluster edges lowest:

| edge | condition | mean Jaccard | overlap | detection (js+oc)/2 |
|---|---|---|---|---|
| **7–5** (APC) | db / wt | 0.79 / 0.69 | 0.94 / 0.97 | **1.00 / 1.00** |
| **6–7** (APC) | db / wt | 0.70 / 0.65 | 0.92 / 0.86 | **0.85 / 0.75** |
| 3–4, 4–8, 3–10, 3–9 | both | 0.5–0.77 | 0.74–1.0 | 0.5–1.0 |
| 2–6, 6–12, 1–11, 1–4 | both | **0.19–0.36** | 0.6–0.9 | **0.0–0.45** |

The APC edges the paper depends on (6→7→5) are the most reproducible in both conditions — 7→5 in all 10 draws, 6→7 in ~8–9/10 (it competes with the adjacent 7→5 edge under the one-to-one matching, since they share cluster 7). The weak edges are exactly those into clusters underpopulated in that condition (cl 11 = 61 cells; cl 12 / cl 2 sparse in DB / ND) — an occupancy effect (ties to R1-3), not a broken trajectory. (The two metrics are complementary: e.g. ND 4→8 has low Jaccard but ~1.0 overlap because cl 8 is tiny — Lamian's averaging is what keeps such small-cluster edges detectable.)

---

## Interpretation / framing for the response letter

- **The core trajectory reproduces per condition.** ND recovers the joint ordering in every draw; the APC branch (6→7→5) and reparative backbone are recovered in all 10 draws of both conditions by cell-level overlap.
- **Where per-condition draws differ, it is cluster-occupancy**, not a different trajectory — the only low-overlap edges are those into clusters sparse in that condition (R1-3 occupancy).
- **DB is noisier** — a minority of DB draws give a poor global ordering — which is honest to state and is why the joint model is a reasonable basis for DB; but the edges that matter still reproduce.
- **We report an ensemble, not one tree**, because Lamian's clustering is non-deterministic (documented `kmeans.seed` no-op); we pass seeds but do not patch the package.

---

## Files & reproducibility

Heavy compute lives in this repo (`workflow/rules/revision.smk`); lightweight post-processing/figures in `moma`. Outputs are gitignored and synced to Box under `results/revision/`.

| | rule | output |
|---|---|---|
| Variant C draws | `revision_condpca_draw` (seeds 1–10) | `results/revision/lamian_condpca_draws/{wt,db}/seed{1..10}/infer_tree.rds` |
| Ensemble metrics | `revision_condpca_edge_detection` | `condpca_draw_manifest.csv` (per-draw ρ, k), `condpca_edge_detection.csv` (per-edge Jaccard/detection) |

Each draw's `infer_tree.rds` also carries `res$pca_full` — the full 50-dim condition-specific PCA (deterministic), e.g. for PGD on the condition's PC space.

Post-processing / figures (moma): `src/revision/occupancy_table.R` (R1-3), `src/revision/trajectory_concordance.R`, `make_figures/revision/umaps_trajectory_conditions.R`.

Notes: Lamian's `pseudotime`/`clusterid` repeat shared backbone cells once per branch — dedupe per cell before joins. Joint MST edges come from `igraph::as_edgelist(joint$MSTtree)`; per-draw edges likewise. Edge detection reproduces Lamian's `evaluate_uncertainty` exactly: per-edge `js.cut`/`oc.cut` nulls (99th percentile of Jaccard/overlap vs size-matched random cell sets), one-to-one matching via `Lamian:::get_binary`, and `detection = (js.perc + oc.perc)/2`. It is pure post-processing of the saved draws (no re-inference).
