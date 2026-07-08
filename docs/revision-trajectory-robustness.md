# Trajectory robustness: per-condition inference (Reviewer 2, R2-2)

**Manuscript:** *Transcriptional regulators predicted to drive macrophage dysregulation during impaired wound healing in diabetic mice*
**Reviewer point:** R2-2 (editor-flagged as the must-address item). The published trajectory was inferred **jointly** on non-diabetic (ND / `wt`) + diabetic (DB / `db`) cells, so apparent per-condition differences might reflect **cluster occupancy** rather than genuinely different trajectories.

**Question.** When each condition is analyzed on its own, does the joint trajectory reproduce?

**Short answer.** For **non-diabetic** cells, yes: an independent re-inference robustly reproduces the joint pseudotime ordering (Spearman ρ ≈ 0.76, stable). For **diabetic** cells, **no stable independent trajectory exists** — the diabetic data admit many near-equally-good clusterings whose trajectories disagree wildly (ρ ranges ≈ 0.16–0.74), so DB is not identifiable on its own. This is itself the answer to the reviewer: DB is too noisy/underpowered to re-infer independently, which is exactly why the **joint** inference is the appropriate analysis for DB, and the per-condition differences that appear are cluster-occupancy differences rather than a genuinely different DB trajectory.

> **Reproducibility note (important).** Lamian's clustering is nondeterministic as shipped — `mykmeans` ignores its `seed` argument and selects the cluster number inside an unseeded parallel fork, so repeated runs on identical input give different trajectories. We patch it to be seeded + well-optimized (`workflow/scripts/revision/reproducible_kmeans.R`). The ND/DB contrast above only becomes visible *after* this fix; earlier single-run numbers (e.g. ND 0.85 / DB 0.74) were individual non-reproducible draws and should not be reported.

---

## Recommended presentation to the reviewer (keep it short)

Lead with the honest, defensible statement — it resolves the concern without over-claiming:

Suggested response text (adapt):

> To test whether the trajectory is robust to the joint ND+DB inference, we re-inferred it independently within each condition (condition-specific embedding and clustering, nothing borrowed from the joint analysis). In non-diabetic cells the independent inference reproduces the joint pseudotime ordering (Spearman ρ ≈ 0.76), confirming the trajectory backbone — including the reparative path and the APC transition — is not an artifact of pooling. Diabetic cells alone do not support a stable independent trajectory: because the diabetic compartment is noisier and several clusters are sparsely populated, independent re-clustering is unstable and does not converge on a single ordering. This is consistent with the reviewer's concern — the apparent per-condition differences reflect **cluster occupancy** (Table S#), not a distinct diabetic trajectory — and is why the joint model is the appropriate basis for the diabetic analysis. We now state this explicitly and present the per-condition inference in Fig. S# / Table S#.

**Supplement:** Fig S# = per-condition inference on the shared layout (reference / ND / DB); Table S# = the occupancy table (R1-3) + the concordance/stability numbers.

**State the DB limitation plainly** — it is the crux, not a footnote: the diabetic trajectory is not independently identifiable; diabetic pseudotime also tracks real time only weakly. Both point to the same thing (DB is noisier), and both support using the joint model for DB.

---

## Data and gate

Cells with a trajectory cluster label (`clusterid`, the Lamian k-means nodes = manuscript clusters 1–12): 6,109 (ND 2,813; DB 3,296). Origin is `D3` (`infer_tree_structure(origin.celltype = "D3")`).

| condition | D3 | D6 | D10 |
|---|---|---|---|
| ND (wt) | 894 | 585 | 1334 |
| DB (db) | 1121 | 1214 | 961 |

Both conditions have ample D3 cells, so the origin is well defined per condition. Note DB's cells skew earlier (more D3/D6, fewer D10) — fewer late cells is part of why the DB late trajectory is unstable.

*(Occupancy also resolves R1-3: **Cluster 12** is a D3 cluster — ~89% D3, ND-dominant (69%) but with a real DB component (20%), i.e. *not* ND-only; **Cluster 1** is DB-dominated (~72%) and concentrated at D6/D10.)*

---

## The test: re-infer everything in a condition-specific PC space (most independent)

The trajectory involves three ingredients that could be held fixed or re-derived per condition: the **PC space**, the **cluster labels**, and the **MST/pseudotime**. The most stringent test re-derives all three within each condition, borrowing nothing from the joint analysis: recompute the embedding on each condition alone (HVGs → scaling → PCA, `npcs = 50`), then run Lamian natively (elbow picks pcadim, k-means re-clusters). Cluster labels and pseudotime are condition-specific and **not** directly comparable across conditions; comparison is via **shared-cell pseudotime rank concordance** against the joint.

### The clustering must be made reproducible first

As shipped, `Lamian:::mykmeans` is nondeterministic (it never uses its `seed`, and its cluster-number search runs in an unseeded `mclapply` fork). Two consequences we verified directly:

- Back-to-back runs on identical data/PCA/seeds returned **12 vs 14 clusters** for DB.
- Concordance with the joint swung across draws from ρ ≈ **−0.12 to 0.83**.

We install a seeded, restart-based (`nstart`) replacement (`reproducible_kmeans.R`). This makes a single run reproducible; it does **not** make DB identifiable (below).

### Result — ND is stable, DB is not

Stability sweep (`revision_condpca_stability`, `results/revision/condpca_stability.csv`): the fixed, well-optimized inference run across 20 k-means seeds each.

| condition | clusters | ρ vs joint (across seeds) | verdict |
|---|---|---|---|
| **ND (wt)** | 9 (all 20 seeds) | **0.762 – 0.767**, median 0.763 | robustly recovers the joint ordering |
| **DB (db)** | 12–13 (varies) | **0.163 – 0.736**, median 0.62 | **non-identifiable** |

The decisive detail for DB: clustering **quality and concordance are decoupled**. Across the sweep the total within-cluster SS spans 306,900–356,900 and ρ spans 0.16–0.74, but the two are *unrelated* — the *lowest*-WSS (best-fitting) clusterings, which pick 13 clusters, give ρ ≈ 0.74, while many 12-cluster clusterings of only slightly higher WSS give ρ ≈ 0.16, and a still-worse-fitting one also gives ρ ≈ 0.74. You cannot predict the diabetic trajectory from the clustering objective. Even with heavy restarts (`nstart = 200`) DB does not converge on a single ordering. There is therefore no defensible single DB concordance value; DB partially and inconsistently recovers the joint ordering (median ρ ≈ 0.62, range 0.16–0.74).

### Consensus across a seed ensemble — the instability is in the ordering, not the clustering

To avoid trusting one seed we save the full tree for 10 seeds per condition and summarise the consensus (`revision_condpca_seed` → `revision_condpca_consensus`; `condpca_consensus_{manifest,summary,labels}.csv` + `lamian_condpca_seeds/{cond}/seed{1..10}/`).

| condition | modal k | ρ (median; range) | mean pairwise ARI | representative seed | best-ρ seed |
|---|---|---|---|---|---|
| **ND (wt)** | 9 | 0.763 (0.762–0.767) | **0.99** | 4 (ρ 0.76) | 9 (ρ 0.77) |
| **DB (db)** | 12 | 0.62 (0.163–0.736) | **0.85** | 2 (ρ 0.62) | 10 (ρ 0.74) |

The DB **adjusted Rand index of 0.85 is the crux**: the diabetic *clustering* is largely consistent across seeds (cells fall into essentially the same groups), so DB instability is **not** a clustering artifact. What is unstable is the **MST / pseudotime built on those clusters** — small partition differences and the 12↔13 cluster-count flip re-wire the branch topology and reorder pseudotime, which is what swings ρ. In other words: the diabetic cell *states* are reproducible, but a single diabetic *ordering* among them is not identifiable.

Practical guidance: for ND any seed is fine (use the representative, seed 4). For DB there is no canonical trajectory — report the ensemble (median ρ ≈ 0.62, range 0.16–0.74). If a single DB tree is needed for a figure, seed 2 is the most-central/consensus tree (ρ ≈ 0.62); seed 10 is the strongest draw (ρ ≈ 0.74, the "July-1-like" trajectory) but is a favorable minority topology (2/10 seeds) and should be labelled as such, not presented as the diabetic trajectory.

---

## Interpretation / recommended framing for the response letter

- **ND independently reproduces the joint trajectory** (ρ ≈ 0.76, stable across clusterings) — the backbone, including the APC transition and reparative path, is not a pooling artifact.
- **DB has reproducible cell states but no identifiable independent trajectory.** The diabetic clustering is consistent across seeds (ARI 0.85), but the branch/pseudotime ordering built on it is not — equally-good clusterings give trajectories ranging from strong to no concordance. The right conclusion is not "DB reproduces the joint" nor "DB contradicts it," but that **DB alone is underpowered to re-infer a stable ordering** — so the joint model is the appropriate basis for the diabetic analysis.
- **The per-condition differences are cluster-occupancy differences**, not a distinct diabetic trajectory (occupancy table, R1-3). This is exactly the reviewer's own framing.
- **DB is noisier** on every axis: unstable independent clustering, weak pseudotime-vs-real-time coupling. Worth stating plainly as a limitation.

### Methodological note (why condition-specific PCA)

A naive alternative — subset each condition into the **joint** PC space **and** re-cluster — makes the DB trajectory look broken. That is an artifact of forcing DB cells into a shared PC space and re-clustering; the condition-specific embedding removes it. (For DB the distinction is moot, since even the condition-specific embedding is non-identifiable.)

---

## Files & reproducibility

Heavy compute lives in this repo (`workflow/rules/revision.smk`), lightweight post-processing/figures in the `moma` repo. Outputs are gitignored and synced to Box under `results/revision/`.

| | rule | output |
|---|---|---|
| Per-condition tree | `revision_infer_tree_condpca` → `revision_evaluate_uncertainty_condpca` | `results/revision/lamian_condpca/{wt,db}/{infer_tree,evaluate_uncertainty}.rds` |
| Stability sweep | `revision_condpca_stability` | `results/revision/condpca_stability.csv` |
| Seed ensemble | `revision_condpca_seed` (seeds 1–10) | `results/revision/lamian_condpca_seeds/{wt,db}/seed{1..10}/infer_tree.rds` |
| Consensus | `revision_condpca_consensus` | `results/revision/condpca_consensus_{manifest,summary,labels}.csv` |
| Reproducible k-means | `workflow/scripts/revision/reproducible_kmeans.R` | (sourced by infer/stability/seed rules) |

The stored `lamian_condpca/{wt,db}/infer_tree.rds` are now reproducible (fixed seed + `nstart`). The DB tree is one arbitrary point in the non-identifiable range — it is illustrative, not canonical; any DB figure must be read alongside the stability sweep, not as "the" diabetic trajectory.

Post-processing / figures (moma): `src/revision/occupancy_table.R` (R1-3), `src/revision/trajectory_concordance.R` (per-condition vs joint), `make_figures/revision/umaps_trajectory_conditions.R` (per-condition on the shared PGD layout).

Notes: `evaluate_uncertainty` is a parallel reimplementation of `Lamian::evaluate_uncertainty` that pre-draws per-iteration seeds (reproducible) and tolerates rare degenerate bootstraps (drops them; denominator = successful permutations). Lamian's `pseudotime`/`clusterid` repeat shared backbone cells once per branch — dedupe per cell before joins. Figure backbones are parsed from `names(res$order)`, not `detection.rate` rownames.
