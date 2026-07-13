#!/usr/bin/env python
"""Visualize a BIRDccNEST run against the Lamian reference (feature branch).

Reads the CSV outputs of birdccnest.py (no RegDiffusion re-run) plus the shared
embedding from the h5ad, and renders one 4-panel figure:

  (A) shared UMAP colored by Lamian clusterid    (B) same UMAP colored by BIRDccNEST community
  (C) overlap heatmap community x clusterid       (D) cluster flow network with the OMST
      (row-normalized, co-ordered to diagonalize)     (root-anchored, laid out by trajectory depth)

Design choices (following the dataviz principles that transfer to a static figure):
  - categorical hues in a FIXED order, and every category is ALSO text-labeled
    (direct labels on the UMAPs / heatmap ticks / network nodes) so identity never
    rests on color alone — necessary because 12-13 categories cannot all be made
    CVD-distinct.
  - the heatmap is a single-hue sequential scale (magnitude), not a rainbow.
  - community colors are shared between panel B and panel D (same entity, same hue).

Usage:
  python birdccnest_viz.py --h5ad results/celloracle/cells.h5ad \
      --outdir results/revision/birdccnest/joint --title "BIRDccNEST joint (res=1.0)"
"""

import argparse

import numpy as np
import pandas as pd
import networkx as nx
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def categorical_colors(ids):
    """Fixed-order categorical hues keyed by id. tab20 gives 20 distinct hues; for
    >20 we wrap (labels disambiguate). Returns {id: rgba}."""
    ids = list(ids)
    cmap = plt.get_cmap("tab20")
    return {k: cmap(i % 20) for i, k in enumerate(ids)}


def scatter_by_label(ax, xy, labels, color_map, title):
    """UMAP scatter colored by a categorical label, with each category direct-labeled
    at its median position (secondary encoding so color isn't the only channel)."""
    labels = np.asarray(labels).astype(str)
    for lab in sorted(pd.unique(labels), key=lambda s: (len(s), s)):
        m = labels == lab
        ax.scatter(xy[m, 0], xy[m, 1], s=3, alpha=0.55, linewidths=0,
                   color=color_map.get(lab, (0.6, 0.6, 0.6, 1.0)))
        cx, cy = np.median(xy[m, 0]), np.median(xy[m, 1])
        ax.text(cx, cy, lab, fontsize=8, fontweight="bold", ha="center", va="center",
                color="black", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def ordered_overlap(grid, comm_pt):
    """Row-normalize the community x label grid, order rows by community mean
    pseudotime, and order columns to push mass onto the diagonal."""
    grid = grid.copy()
    grid.index = grid.index.astype(int)
    row_order = comm_pt.reindex(grid.index).sort_values().index
    grid = grid.loc[row_order]
    rownorm = grid.div(grid.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    # column centroid = weighted mean row-position; sort columns ascending
    pos = np.arange(len(rownorm))
    centroid = {c: (pos * rownorm[c].values).sum() / max(rownorm[c].values.sum(), 1e-9)
                for c in rownorm.columns}
    col_order = sorted(rownorm.columns, key=lambda c: centroid[c])
    return rownorm[col_order]


def heatmap(ax, rownorm, title, label_col):
    im = ax.imshow(rownorm.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(rownorm.shape[1])); ax.set_xticklabels(rownorm.columns, fontsize=8)
    ax.set_yticks(range(rownorm.shape[0])); ax.set_yticklabels(rownorm.index, fontsize=8)
    ax.set_xlabel(f"Lamian {label_col}", fontsize=9)
    ax.set_ylabel("BIRDccNEST community\n(ordered by pseudotime)", fontsize=9)
    ax.set_title(title, fontsize=10)
    # annotate dominant cells so the alignment is legible without reading color alone
    for i in range(rownorm.shape[0]):
        for j in range(rownorm.shape[1]):
            v = rownorm.values[i, j]
            if v >= 0.15:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if v > 0.5 else "black")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="fraction of community")


def layered_positions(tree, root):
    """x = trajectory depth (shortest-path from root over the rooted tree); y spread
    within depth. Shows the trajectory flowing outward from the origin community."""
    depth = nx.shortest_path_length(tree.to_undirected(), root)
    by_depth = {}
    for n, d in depth.items():
        by_depth.setdefault(d, []).append(n)
    pos = {}
    for d, nodes in by_depth.items():
        nodes = sorted(nodes)
        for i, n in enumerate(nodes):
            pos[n] = (d, i - (len(nodes) - 1) / 2.0)
    return pos


def network(ax, outdir, comm_color, dom_cluster, comm_size, comm_pt, title):
    flow = pd.read_csv(f"{outdir}/cluster_flow_edges.csv")
    try:
        omst = pd.read_csv(f"{outdir}/omst_rooted_edges.csv")
    except FileNotFoundError:
        omst = pd.read_csv(f"{outdir}/omst_edges.csv")
    root = int(comm_pt.idxmin())
    tree = nx.DiGraph()
    tree.add_edges_from([(int(r.u), int(r.v)) for r in omst.itertuples()])
    for n in comm_size.index:
        tree.add_node(int(n))
    pos = layered_positions(tree, root)

    # faint flow edges underneath (all inter-community edges), width ~ weight
    wmax = flow.weight.max()
    for r in flow.itertuples():
        u, v = int(r.u), int(r.v)
        if u in pos and v in pos:
            ax.annotate("", xy=pos[v], xytext=pos[u],
                        arrowprops=dict(arrowstyle="-", color="0.85",
                                        lw=0.3 + 2.0 * r.weight / wmax, alpha=0.5))
    # bold OMST edges on top, with arrows (the inferred trajectory)
    for r in omst.itertuples():
        u, v = int(r.u), int(r.v)
        ax.annotate("", xy=pos[v], xytext=pos[u],
                    arrowprops=dict(arrowstyle="-|>", color="0.25", lw=2.0,
                                    shrinkA=12, shrinkB=12))
    # nodes: color = community (shared with panel B), size ~ #cells, label = dominant cluster
    smax = comm_size.max()
    for n in tree.nodes():
        x, y = pos[n]
        ax.scatter([x], [y], s=120 + 700 * comm_size.get(n, 0) / smax,
                   color=comm_color.get(str(n), (0.6, 0.6, 0.6, 1.0)),
                   edgecolors="black", linewidths=0.6, zorder=3)
        ax.text(x, y, f"{n}\ncl{dom_cluster.get(n, '?')}", ha="center", va="center",
                fontsize=6.5, fontweight="bold", zorder=4)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("trajectory depth from origin community", fontsize=9)
    ax.set_yticks([])
    ax.margins(x=0.14, y=0.16)  # keep large nodes/labels off the panel edges
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--title", default="BIRDccNEST")
    ap.add_argument("--label-col", default="clusterid")
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    umap = np.asarray(adata.obsm["X_umap"])
    cc = pd.read_csv(f"{args.outdir}/cell_communities.csv").set_index("cell")
    cc = cc.reindex(adata.obs_names)
    comm_pt = pd.read_csv(f"{args.outdir}/community_mean_pseudotime.csv", index_col=0).iloc[:, 0]
    comm_pt.index = comm_pt.index.astype(int)
    grid = pd.read_csv(f"{args.outdir}/community_overlap_{args.label_col}.csv", index_col=0)

    # shared categorical maps
    comm_ids = sorted(cc["community"].dropna().astype(int).unique())
    comm_color = {str(k): v for k, v in categorical_colors(comm_ids).items()}
    lamian_ids = sorted(adata.obs[args.label_col].astype(str).unique(), key=lambda s: (len(s), s))
    lamian_color = categorical_colors(lamian_ids)

    dom_cluster = (cc.dropna(subset=["community"])
                     .groupby(cc["community"].astype("Int64"))[args.label_col]
                     .agg(lambda s: s.value_counts().index[0]))
    comm_size = cc["community"].dropna().astype(int).value_counts().sort_index()

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    scatter_by_label(axes[0, 0], umap, adata.obs[args.label_col].values, lamian_color,
                     f"UMAP · Lamian {args.label_col}")
    scatter_by_label(axes[0, 1], umap, cc["community"].astype("Int64").astype(str).values,
                     comm_color, "UMAP · BIRDccNEST community")
    heatmap(axes[1, 0], ordered_overlap(grid, comm_pt),
            "Community → cluster overlap (row-normalized)", args.label_col)
    network(axes[1, 1], args.outdir, comm_color, dom_cluster.to_dict(),
            comm_size, comm_pt, "Cluster flow network + OMST trajectory")

    fig.suptitle(args.title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = f"{args.outdir}/birdccnest_figure.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    main()
