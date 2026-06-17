# MoMa Diabetic Wound Healing — npj Revision Meeting

**Date:** 2026-06-17
**Attendees:** Timothy Koh (Dr. Koh), Yang Dai (Dr. Dai), Brandon Lukas, Jingbo Pang (King Bo)
**Topic:** Planning responses to reviewer comments for the npj (Systems Biology and Applications) revision of the macrophage diabetic wound-healing manuscript.

---

## Summary

The team walked through the reviewer comments one by one, deciding which require new analysis vs. wording/discussion changes. The reviewer appears to be a biology-oriented reader (possibly the same reviewer from the earlier Frontiers paper, given they singled out *ZBTB46* and mis-stated the journal name). The reviews were generally regarded as fair and reasonable. Most comments are addressable with clarifications and modest additional analysis. A second, recently-published IRF4 paper and a soon-to-publish "Apc-like/Light" macrophage paper can be leaned on as supporting literature.

**Submission deadline:** July 30 (an extension is available if needed).
**Drafting plan:** Brandon takes the first crack at the response letter (especially the computational/analysis comments); Dr. Koh handles the biology-heavy comments and will share his in-progress Frontiers-paper response letter as a template.

---

## Comment-by-comment outcomes

### Comment 1 — In-silico knockouts / transcription factor activity at "time points"
- Reviewer (a biology reader) is likely conflating **time points with pseudotime**, asking whether Cell Oracle in-silico knockout predictions align with observed/experimental data.
- Example discussed: knocking out **IRF4** inhibits the cluster-7 → APC-like transition; diabetics have a reduced APC-like population. If observed RNA data shows IRF4 downregulation, that supports the Cell Oracle result.
- Key point to make in the response: **TF activity ≠ expression level** (BITFAM logic), and similar overall transcriptome in a cluster doesn't imply identical individual-gene expression.
- Plan: Use the **combined-cell analysis**, then retrieve the **simulated knockout profiles** for DB vs. ND cells and compare; frame as "simulated profile closer to transcriptome at a certain pseudotime point."
- Verdict: **Easy to address.** Reviewer just wants "something."

### Comment 2 — "Validated" / mechanistic claims wording
- Reviewer objects to language implying mechanistic/validated results without experimental validation.
- Plan: Replace **"validated" → "consistent with"** (e.g., "consistent with the literature finding"); avoid "literature-validated" and "experimental validation." Change "driving transitions" language to **"candidate regulators / candidate drivers"** and emphasize "predicted."
- Verdict: **No new analysis — wording changes only.** Lean on the now-published IRF4 paper as support.

### Comment 3 — Trajectory inference robustness
- Reviewer wants additional analysis demonstrating robustness of trajectories (separate ND/DB analysis and/or comparison with an alternative trajectory-inference method).
- Context: Earlier non-diabetic-mice (Frontiers) paper had a fairly similar trajectory structure (early→late reparative, early→APC). Because the macrophages are well quality-controlled and not too heterogeneous, trajectories are not expected to "go too crazy."
- Plan: **Run the separate/alternative trajectory analysis to satisfy the reviewer.** If results are clean and interpretable → put in a **supplementary figure**; if messy → make an argument against the approach. Either outcome gives a defensible response.
- Verdict: **Requires new analysis** (low risk).

### Comment 4 — Justification for the cell-data / population (e.g., ZBTB46, "Light/Apc-like" macrophages)
- This population is gaining interest (Dr. Koh got comments at the recent meeting; another group publishing soon).
- Plan: **No new analysis.** Add discussion, collate more marker genes describing the phenotype, and lean on the newly published paper.
- Verdict: **Discussion/wording only.**

### Comment 5 — PGD (Brandon)
- Reviewer asks whether PGD is only for visualization and requests **benchmarking** of PGD.
- Brandon felt it's slightly out of scope (PGD is just an application here) but plans to do this benchmarking anyway, so he doesn't mind addressing it.
- Plan: Emphasize Cell Oracle results are a **starting point for future work**; expand the **methods** (technical detail on Cell Oracle), possibly in the **supplement**, since Cell Oracle itself doesn't provide formal statistical significance assessment.
- Verdict: **Addressable** — benchmarking + methods expansion.

### Comment 6 — Gene regulatory network / possible circularity
- Reviewer worried the GRN reasoning may be circular.
- Clarification reached in the meeting: BITFAM (TF activity) and Cell Oracle both ultimately trace to the **same seed/base GRN derived from ChIP-Atlas** as prior, but the analyses are **independent and not circular** — Cell Oracle is **not** built from BITFAM results.
- Note of concern (Dr. Dai): using different GRNs for reporting TF activity (BITFAM) vs. the Cell Oracle analysis could invite criticism; the BITFAM results focus only on **TF activity, not targets**.
- Plan: **Clarify** the relationship concisely without over-explaining. Dai/Brandon to draft this response.
- Verdict: **Clarification only.**

---

## Action Items

| # | Action | Owner | Notes / Deadline |
|---|--------|-------|------------------|
| 1 | Draft the full reviewer response letter (first pass) | **Brandon** | Lead on computational/analysis comments |
| 2 | Share in-progress Frontiers-paper response letter as a template/example | **Dr. Koh** | For Brandon to match framing |
| 3 | Review & write biology-heavy comment responses; dig into literature (incl. the ~10-yr-old CEBPA paper) | **Dr. Koh** | Comments 1, 4 |
| 4 | **Comment 1:** Run combined-cell Cell Oracle, retrieve simulated knockout profiles for DB vs. ND, compare; check IRF4 expression in cluster 7 (and earlier) | **Brandon** | Frame vs. pseudotime; note TF activity ≠ expression |
| 5 | **Comment 2:** Replace "validated" with "consistent with"; change "driving transitions" → "candidate regulators/drivers"; emphasize "predicted" | **Brandon / Dr. Koh** | Wording only; cite published IRF4 paper |
| 6 | **Comment 3:** Run separate ND/DB and/or alternative trajectory-inference analysis; decide supplementary figure vs. counter-argument | **Brandon** (w/ Jingbo input) | New analysis |
| 7 | **Comment 4:** Add discussion + collate marker genes for the Light/Apc-like population; lean on newly published paper | **Dr. Koh** | Discussion only |
| 8 | **Comment 5:** Benchmark PGD; expand Cell Oracle methods (possibly supplement); frame results as starting point for future work | **Brandon** | Benchmarking + methods |
| 9 | **Comment 6:** Write clarification that BITFAM and Cell Oracle share the ChIP-Atlas-derived seed GRN but are non-circular/independent; note BITFAM focuses on TF activity, not targets | **Brandon / Dr. Dai** | Clarification only |
| 10 | Set up recurring revision meeting — **every other week** initially, moving to weekly as deadline nears | **Brandon** | Same time slot |
| 11 | Track submission deadline | All | **July 30** (extension available if needed) |

---

## Other notes
- Congrats to **Brandon** on a graduate research award (~$2,000) usable toward his next conference.
- Brandon's late-breaking abstract submitted to **ISMB** (likely poster).
- The **PGD** paper accepted for proceedings at **IEEE EMBC**; the bioRxiv preprint has ~300+ downloads.
- Reviewers/journal regarded as fair; review took a long time but comments are reasonable.
