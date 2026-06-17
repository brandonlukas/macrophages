# Reviewer Revision — Computational/Bioinformatics Action Items

**Manuscript:** *Transcriptional regulators predicted to drive macrophage dysregulation during impaired wound healing in diabetic mice*
**Decision:** Revise & resubmit (Reviewer 2: major revision). Submission ID `da864ce4-898c-411d-80fb-43abccbaa655`.
**Scope of this doc:** Only the analysis/computation items. Pure writing/tone items (e.g. citing Wicks et al. 2015, softening "drive"/"validated" language, Lyve1 balanced-discussion framing) are tracked separately — they are *response-letter* tasks, not analyses, though several analyses below feed those rewrites.

Editor (M. Craig) explicitly flagged **trajectory inference robustness (R2-2)** as the must-address point.

---

## Priority 1 — Editor-flagged / central conclusions

### T1. Trajectory robustness: separate ND vs DB inference (R2 — trajectory)
The current trajectory was inferred jointly on ND+DB cells, so differences may reflect cluster occupancy, not real per-condition trajectories.
- [ ] Re-run trajectory inference (Lamian) **separately on ND-only and DB-only** cells; compare branch structure and pseudotime ordering against the joint trajectory.
- [ ] OR/AND benchmark against an **alternative TI method** (e.g. Slingshot, PAGA, Monocle3) to show the branch topology is method-robust.
- [ ] Explicitly report which conclusions hold per-condition vs only in the joint model; reframe any that don't as cluster-occupancy differences.

### T2. ND-only vs DB-only differential TR activity / expression (R1-5)
Reviewer 1 wants direct evidence that specific TRs are actually altered in diabetic cells, not just KO predictions on pooled cells.
- [ ] Test whether highlighted TRs (Cebpa, Irf8, Irf4, Pparg, Nr1h3, Nr3c1, Ogt, Brd2, Brd3, Zmynd11, Zbtb46) are **significantly under-/over-expressed or differential in BITFAM activity** between ND and DB at the critical timepoints (D3/D6/D10).
- [ ] Present as a per-timepoint ND-vs-DB comparison (expression + inferred activity) so prediction is anchored to observed dysregulation.

---

## Priority 2 — Method validation & statistical rigor

### T3. CellOracle perturbation statistical assessment (R2 — CellOracle)
Perturbation scores / Markov simulations are currently descriptive.
- [ ] Add a **formal robustness/significance assessment** of perturbation scores (e.g. permutation/null distribution, bootstrap over cells, ranking stability across reruns).
- [ ] Reframe Results to present ranked regulators as a **prioritization framework**, not established effects.

### T4. GRN dependence & circularity (R2 — network)
BITFAM and CellOracle both lean on the same ChIP-Atlas-derived network → potential circularity.
- [ ] **Sensitivity to alternative regulatory priors** — re-run with a different base GRN and report overlap of top hits.
- [ ] Quantify the **impact of the manual Zbtb46 target inclusion** (results with vs without it).
- [ ] Add network-robustness discussion supported by the above.

### T5. PGD benchmarking (R2 — PGD) — *you flagged as possibly out of scope; doing anyway*
Reviewer asks for validation that PGD adds value beyond standard visualization. Three parts:
- [ ] **Smoothness:** quantify pseudotime-gradient smoothness in 2D, PGD vs UMAP (Moran's I of pseudotime over the embedding; local pseudotime SD among 2D kNN; rank corr of 2D distance vs |Δpseudotime|).
- [ ] **Structure preservation** (defeats the "over-smoothing" critique): trustworthiness & continuity, kNN overlap, cluster silhouette — PGD vs UMAP relative to high-dim/PCA space. Headline: *smoother gradient without worse neighborhood preservation.*
- [ ] **Parameter sensitivity:** vary diffusion time t / k / bandwidth, show layout + smoothness metric stable.
- [ ] **Reframe PGD as visualization-only:** state no biological conclusion depends on it (claims come from Lamian/BITFAM/CellOracle); reproduce key figures on plain UMAP in supplement.

---

## Priority 3 — Targeted data checks (small, fast)

### T6. Verify & fix contradictory cluster compositions (R1-3)
Conflicting statements about whether **Cluster 12** is ND-D3 vs ND+DB-D3, and similar for **Cluster 1**.
- [ ] Recompute per-cluster composition by condition × timepoint (occupancy table/fractions); establish the correct statement and reconcile all text against it.

### T7. Cell-state annotation support (R2 — annotations)
- [ ] **Cluster 5 (APC):** marker analysis distinguishing monocyte-derived DC vs activated antigen-presenting macrophage (Ciita, MHC-II, Cd209a + DC-specific vs macrophage markers) — matters because of downstream Zbtb46 emphasis.
- [ ] **Cluster 1 (foam-cell-like):** check expression of additional lipid-associated macrophage markers (e.g. Trem2, Lpl, Cd36, Gpnmb, Plin2) to strengthen the call.

### T8. Supplemental video annotations (R1-4)
- [ ] Clarify what each CellOracle simulation video shows: which KO (e.g. Cebpa), what the left/right branched timelines represent, whether Irf8 is included; regenerate/annotate as needed so the figure legend matches the text.

---

## Not analyses (response-letter only — listed for completeness)
- R1-1: cite Wicks et al., *Diabetes* 2015 (Cebpa support).
- R1-2: clarify wounds-per-animal / biological replicate structure (Methods text).
- R2-1: soften causal "drive" language (Abstract/Results/Discussion).
- R2-2-lang: distinguish computational prediction vs literature-consistency vs experimental validation; stop calling predictions "validated."
- R2-annot: more balanced Lyve1+ discussion (homeostasis/repair contexts).
