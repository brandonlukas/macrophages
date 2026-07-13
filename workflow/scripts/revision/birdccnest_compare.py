#!/usr/bin/env python
"""Compare the BIRDccNEST OMST to the Lamian MST over the SAME clusters.

Run after `birdccnest.py --cluster-col clusterid` (which builds the flow network +
OMST over the Lamian clusters) and after exporting the Lamian MST edges to
`lamian_mst_edges.csv` in the same outdir. Produces a 2-panel figure:

  (A) UMAP overlay — cluster centroids as nodes; edges colored by agreement
      (shared / BIRDccNEST-only / Lamian-only); BIRDccNEST edges carry direction
      arrows; the APC branch (6-7, 5-7) is highlighted.
  (B) BIRDccNEST OMST as a tidy tree on a real pseudotime x-axis, edges colored by
      the same agreement categories — shows the inferred trajectory + orientation.

Design: agreement is encoded by BOTH color AND line style (shared = solid dark,
method-specific = dashed), so it survives CVD/grayscale; a legend is always shown;
the APC branch gets a red halo as a third, redundant channel.

Usage:
  python birdccnest_compare.py --h5ad results/celloracle/cells.h5ad \
      --outdir results/revision/birdccnest/joint_lamianclusters
"""

import argparse

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SHARED = "#2c3e50"       # dark slate — agreement
BIRD_ONLY = "#e67e22"    # orange — BIRDccNEST-specific
LAMIAN_ONLY = "#7d3c98"  # purple — Lamian-specific
APC = "#e74c3c"          # red halo — the APC branch
APC_EDGES = {frozenset((6, 7)), frozenset((5, 7))}


def edge_category(u, v, shared):
    return SHARED if frozenset((u, v)) in shared else BIRD_ONLY


def draw_edge(ax, p0, p1, color, directed, style="-", lw=2.4, z=3, halo=False):
    if halo:  # red underlay marking the APC branch
        ax.annotate("", xy=p1, xytext=p0, zorder=z - 1,
                    arrowprops=dict(arrowstyle="-", color=APC, lw=lw + 4, alpha=0.55,
                                    shrinkA=16, shrinkB=16))
    arrow = "-|>" if directed else "-"
    ax.annotate("", xy=p1, xytext=p0, zorder=z,
                arrowprops=dict(arrowstyle=arrow, color=color, lw=lw, linestyle=style,
                                shrinkA=16, shrinkB=16, mutation_scale=18))


def tidy_tree_y(children, root):
    """Standard tidy-tree y: leaves get successive slots, parents center on children."""
    y, counter = {}, [0]
    def rec(n):
        ch = children.get(n, [])
        if not ch:
            y[n] = counter[0]; counter[0] += 1
        else:
            for c in ch:
                rec(c)
            y[n] = float(np.mean([y[c] for c in ch]))
    rec(root)
    return y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label-col", default="clusterid")
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    umap = np.asarray(adata.obsm["X_umap"])
    lab = adata.obs[args.label_col].astype(int).values

    bird = pd.read_csv(f"{args.outdir}/omst_rooted_edges.csv")          # directed
    lam = pd.read_csv(f"{args.outdir}/lamian_mst_edges.csv")            # undirected
    comm_pt = pd.read_csv(f"{args.outdir}/community_mean_pseudotime.csv", index_col=0).iloc[:, 0]
    comm_pt.index = comm_pt.index.astype(int)

    B_dir = [(int(r.u), int(r.v)) for r in bird.itertuples()]
    B_und = {frozenset(e) for e in B_dir}
    L_und = {frozenset((int(r.a), int(r.b))) for r in lam.itertuples()}
    shared, l_only = B_und & L_und, L_und - B_und

    clusters = sorted(np.unique(lab))
    centroid = {c: np.array([np.median(umap[lab == c, 0]), np.median(umap[lab == c, 1])]) for c in clusters}
    size = {c: int((lab == c).sum()) for c in clusters}
    smax = max(size.values())

    fig, (axU, axT) = plt.subplots(1, 2, figsize=(20, 9))

    # ---- Panel A: UMAP overlay ----
    axU.scatter(umap[:, 0], umap[:, 1], s=3, c="0.85", linewidths=0, zorder=0)
    for e in l_only:  # Lamian-only: dashed purple, undirected
        u, v = tuple(e)
        draw_edge(axU, centroid[u], centroid[v], LAMIAN_ONLY, directed=False, style=(0, (4, 3)))
    for u, v in B_dir:  # BIRDccNEST edges: directed, solid; shared=dark, else orange
        draw_edge(axU, centroid[u], centroid[v], edge_category(u, v, shared), directed=True,
                  halo=frozenset((u, v)) in APC_EDGES)
    for c in clusters:
        p = centroid[c]
        axU.scatter([p[0]], [p[1]], s=260 + 900 * size[c] / smax, color="white",
                    edgecolors="black", linewidths=1.3, zorder=4)
        axU.text(p[0], p[1], str(c), ha="center", va="center", fontsize=11, fontweight="bold", zorder=5)
    axU.set_title("Trajectory over the 12 Lamian clusters, on the UMAP\n"
                  f"BIRDccNEST OMST vs Lamian MST — {len(shared)}/{len(L_und)} edges shared", fontsize=12)
    axU.set_xticks([]); axU.set_yticks([])
    for s in axU.spines.values():
        s.set_visible(False)

    # ---- Panel B: BIRDccNEST OMST on a pseudotime axis ----
    root = int(comm_pt.idxmin())
    children = {}
    for u, v in B_dir:
        children.setdefault(u, []).append(v)
    ty = tidy_tree_y(children, root)
    posT = {c: np.array([comm_pt.get(c, np.nan), ty.get(c, 0)]) for c in clusters}
    for u, v in B_dir:
        draw_edge(axT, posT[u], posT[v], edge_category(u, v, shared), directed=True,
                  halo=frozenset((u, v)) in APC_EDGES)
    for c in clusters:
        if not np.isfinite(posT[c][0]):
            continue
        axT.scatter([posT[c][0]], [posT[c][1]], s=260 + 900 * size[c] / smax, color="white",
                    edgecolors="black", linewidths=1.3, zorder=4)
        axT.text(posT[c][0], posT[c][1], str(c), ha="center", va="center", fontsize=11, fontweight="bold", zorder=5)
    axT.set_title("BIRDccNEST OMST (directed), laid out by cluster mean pseudotime", fontsize=12)
    axT.set_xlabel("cluster mean pseudotime  →", fontsize=10)
    axT.set_yticks([])
    axT.margins(x=0.10, y=0.14)
    for s in ["top", "right", "left"]:
        axT.spines[s].set_visible(False)

    legend = [
        Line2D([0], [0], color=SHARED, lw=2.6, label=f"shared ({len(shared)})"),
        Line2D([0], [0], color=BIRD_ONLY, lw=2.6, label=f"BIRDccNEST-only ({len(B_und - L_und)})"),
        Line2D([0], [0], color=LAMIAN_ONLY, lw=2.6, ls=(0, (4, 3)), label=f"Lamian-only ({len(l_only)})"),
        Line2D([0], [0], color=APC, lw=6, alpha=0.55, label="APC branch (6→7→5)"),
    ]
    axU.legend(handles=legend, loc="upper left", fontsize=10, frameon=True)

    fig.suptitle("BIRDccNEST vs Lamian trajectory over shared clusters (joint)", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = f"{args.outdir}/birdccnest_vs_lamian_trajectory.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"[compare] wrote {out}  (shared={len(shared)}/{len(L_und)}, "
          f"APC 6-7 & 5-7 recovered: {APC_EDGES <= B_und})")


if __name__ == "__main__":
    main()
