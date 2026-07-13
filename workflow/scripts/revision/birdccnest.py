#!/usr/bin/env python
"""BIRDccNEST — BI-directional RegDiffusion Cell-Cell Networks EStimating Trajectories.

Exploratory alternative trajectory-inference method for the npj revision (feature
branch; NOT wired into the main pipeline yet). Reimplements Cicekli et al. (MLCB
2025, "BIRDccNEST: Interpretable single cell characterization with inferred
directed cell networks") from the authors' Colab + paper.

Idea: RegDiffusion normally infers a gene-gene regulatory network. If instead you
run it on the TRANSPOSED expression matrix, the cells become the network nodes and
you get a directed cell-cell similarity network. Community detection over that
network yields cell clusters/states, and a directed "cluster flow network" (edges
between communities weighted by the number of cell-cell edges crossing them) yields
a trajectory as an oriented maximum spanning tree (OMST, Edmonds' arborescence).

Pipeline (paper Section 2):
  1. transpose normalized expression -> genes x cells (cells = nodes)
  2. RegDiffusion -> asymmetric cell x cell adjacency A (n x n)
  3. prune to the top-weight edges, kept high enough to stay strongly connected
     (paper: top quartile for all datasets except mDC, which used top 40%)
  4. Louvain communities on the pruned directed graph
  5. cluster flow network: node per community; edge u->v weighted by the count of
     pruned cell->cell edges from community u to community v
  6. OMST = nx.maximum_spanning_arborescence(flow network) -> directed trajectory

Caveats we handle explicitly (both flagged in the paper):
  - Orientation is fragile: Edmonds roots the arborescence at high-outgoing-weight
    nodes, and community-size imbalance can reverse the whole trajectory. We do NOT
    trust the raw OMST direction; we report direction concordance against the joint
    Lamian pseudotime and identify the biologically-expected (origin/D3-dominant)
    root community.
  - Sample regime: transposed, RegDiffusion's "samples" are genes. Our data has more
    cells (nodes) than genes (samples), unlike the BEELINE datasets. Use as many
    informative genes as available (the --h5ad export controls this).

Usage:
  python birdccnest.py --h5ad results/celloracle/cells.h5ad \
      --outdir results/revision/birdccnest/joint \
      [--top-percent 0.25] [--resolution 1.0] [--n-steps 1000] [--seed 42] [--device cuda]
"""

import argparse
import os
import random

import numpy as np
import pandas as pd
import scipy.sparse as sp
import networkx as nx
import anndata as ad
import regdiffusion as rd
from sklearn.metrics import davies_bouldin_score, silhouette_score


# ----------------------------------------------------------------------------------
# Reproducibility
# ----------------------------------------------------------------------------------
def set_seed(seed):
    """Seed python / numpy / torch RNGs. NB: RegDiffusion on GPU may retain minor
    nondeterminism (cuDNN); we still pin every seed we can and report it."""
    random.seed(seed)
    np.random.seed(seed)
    import torch
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# ----------------------------------------------------------------------------------
# Step 2: RegDiffusion cell-cell network
# ----------------------------------------------------------------------------------
def infer_cell_network(adata, n_steps, device, seed):
    """Run RegDiffusion on the transposed matrix -> (adj_matrix, cell_names).

    adata is cells x genes (normalized). We transpose so RegDiffusion sees samples =
    genes and variables (nodes) = cells. adj[i, j] is a directed cell_i -> cell_j
    weight."""
    adata_t = adata.transpose()  # genes x cells
    X = adata_t.X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    X = X.astype(np.float32)

    set_seed(seed)
    trainer = rd.RegDiffusionTrainer(X, device=device, n_steps=n_steps)
    trainer.train()

    cell_names = np.asarray(adata_t.var_names)
    cnet = trainer.get_grn(cell_names)
    adj = np.asarray(cnet.adj_matrix, dtype=np.float32)
    names = cnet.gene_names  # "genes" here are cells (RegDiffusion's naming)
    if len(names) and isinstance(names[0], bytes):
        names = [x.decode("utf-8") for x in names]
    return adj, np.asarray(names)


# ----------------------------------------------------------------------------------
# Step 3: prune to top-weight edges, kept high enough to stay strongly connected
# ----------------------------------------------------------------------------------
def prune_to_strongly_connected(adj, top_percent, escalate=(0.25, 0.30, 0.35, 0.40, 0.50)):
    """Zero out all but the top `top_percent` of edge weights; if the resulting
    directed graph is not strongly connected, escalate the kept fraction (paper's
    rule: keep the top quartile 'so long as the resulting network is strongly
    connected', top 40% for mDC). Returns (sparse_pruned_adj, kept_fraction)."""
    from scipy.sparse.csgraph import connected_components
    fractions = [top_percent] + [f for f in escalate if f > top_percent]
    n = adj.shape[0]
    off_diag = ~np.eye(n, dtype=bool)
    for frac in fractions:
        cut = np.quantile(adj, 1 - frac)
        keep = (adj >= cut) & off_diag  # drop self-loops; they carry no flow
        sparse = sp.csr_matrix(np.where(keep, adj, 0.0))
        # strong connectivity via scipy (far cheaper than building a networkx graph)
        n_comp, _ = connected_components(sparse, directed=True, connection="strong")
        if n_comp == 1:
            return sparse, frac, True
    # none strongly connected: return the largest kept fraction, flag it
    return sparse, fractions[-1], False


# ----------------------------------------------------------------------------------
# Steps 4-6: communities -> cluster flow network -> OMST
# ----------------------------------------------------------------------------------
def louvain_communities(sparse_adj, names, resolution, seed):
    """Louvain over the pruned cell-cell network. Returns comm_dict {cell_name: idx}
    and the list of community sets."""
    G = nx.from_scipy_sparse_array(sparse_adj, create_using=nx.DiGraph)
    G = nx.relabel_nodes(G, {i: name for i, name in enumerate(names)})
    comms = nx.community.louvain_communities(G, resolution=resolution, seed=seed)
    comm_dict = {}
    for ci, comm in enumerate(comms):
        for node in comm:
            comm_dict[node] = ci
    return G, comm_dict, comms


def make_cluster_flow_network(G, comm_dict):
    """Directed community graph; edge u->v weight = # of cell->cell edges from a cell
    in community u to a cell in community v (u != v). Faithful to the paper/notebook."""
    from collections import defaultdict
    w = defaultdict(int)
    for u, v in G.edges():
        cu, cv = comm_dict[u], comm_dict[v]
        if cu != cv:
            w[(cu, cv)] += 1
    flowG = nx.DiGraph()
    for (cu, cv), weight in w.items():
        flowG.add_edge(cu, cv, weight=weight)
    return flowG


def extract_omst(flowG):
    """Oriented maximum spanning tree = maximum spanning arborescence (Edmonds)."""
    return nx.maximum_spanning_arborescence(flowG)


# ----------------------------------------------------------------------------------
# Evaluation / validation against the existing joint labels
# ----------------------------------------------------------------------------------
def overlap_grid(comm_dict, label_series):
    """community x label count grid (rows = community idx, cols = label values)."""
    df = pd.DataFrame({"community": pd.Series(comm_dict),
                       "label": label_series.reindex(list(comm_dict.keys())).values})
    grid = (df.groupby(["community", "label"]).size().unstack(fill_value=0)
              .sort_index().sort_index(axis=1))
    return grid


def coherence(embedding, labels):
    """DBI (lower better) + silhouette (higher better) on an embedding for a labeling.
    Mirrors the paper's cluster-quality comparison."""
    labs = np.asarray(labels)
    mask = pd.notnull(labs)
    if len(np.unique(labs[mask])) < 2:
        return np.nan, np.nan
    return (davies_bouldin_score(embedding[mask], labs[mask]),
            silhouette_score(embedding[mask], labs[mask]))


def root_anchored_orientation(omst, root):
    """The raw OMST orientation is set by Edmonds' rooting, which the paper notes can
    flip the whole trajectory under community-size imbalance. We keep the OMST's
    UNDIRECTED backbone (its max-weight tree structure) but re-orient every edge to
    point away from a biologically-chosen `root` community (BFS). Returns a DiGraph."""
    undirected = omst.to_undirected()
    oriented = nx.DiGraph()
    oriented.add_nodes_from(undirected.nodes(data=True))
    for u, v in nx.bfs_edges(undirected, root):
        w = undirected[u][v].get("weight")
        oriented.add_edge(u, v, weight=w)
    return oriented


def direction_concordance(omst, comm_dict, pseudotime):
    """For each OMST edge u->v, does mean pseudotime(u) < mean pseudotime(v)? Reports
    the fraction of edges oriented with increasing pseudotime — a check on the OMST's
    orientation independent of Edmonds' rooting."""
    pt = pd.Series(pseudotime)
    cell_comm = pd.Series(comm_dict, name="community")
    comm_pt = pt.reindex(cell_comm.index).groupby(cell_comm).mean()
    rows, agree = [], 0
    for u, v, d in omst.edges(data=True):
        inc = comm_pt.get(u, np.nan) < comm_pt.get(v, np.nan)
        agree += int(bool(inc))
        rows.append({"u": u, "v": v, "weight": d.get("weight"),
                     "pt_u": comm_pt.get(u), "pt_v": comm_pt.get(v),
                     "increasing_pt": bool(inc)})
    frac = agree / max(len(rows), 1)
    return frac, pd.DataFrame(rows), comm_pt


# ----------------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------------
def build_or_load_graph(args):
    """Return (sparse_adj, names, kept_frac, strong, adata). RegDiffusion training +
    pruning is the expensive part; cache it so Louvain/OMST can be re-swept cheaply."""
    adata = ad.read_h5ad(args.h5ad)
    print(f"[BIRDccNEST] {adata.shape[0]} cells x {adata.shape[1]} genes from {args.h5ad}")
    cache_adj = f"{args.outdir}/pruned_adj.npz"
    cache_names = f"{args.outdir}/cell_names.npy"
    if args.from_cache and os.path.exists(cache_adj):
        print(f"[BIRDccNEST] loading cached pruned network from {cache_adj}")
        sparse_adj = sp.load_npz(cache_adj)
        names = np.load(cache_names, allow_pickle=True)
        from scipy.sparse.csgraph import connected_components
        strong = connected_components(sparse_adj, directed=True, connection="strong")[0] == 1
        return sparse_adj, names, args.top_percent, strong, adata

    print(f"[BIRDccNEST] RegDiffusion ({args.n_steps} steps, {args.device}) on transposed matrix...")
    adj, names = infer_cell_network(adata, args.n_steps, args.device, args.seed)
    print(f"[BIRDccNEST] cell-cell adjacency: {adj.shape[0]} x {adj.shape[1]}")
    sparse_adj, kept_frac, strong = prune_to_strongly_connected(adj, args.top_percent)
    print(f"[BIRDccNEST] pruned to top {kept_frac:.0%} edges "
          f"({sparse_adj.nnz} edges; strongly connected: {strong})")
    sp.save_npz(cache_adj, sparse_adj)
    np.save(cache_names, names)
    return sparse_adj, names, kept_frac, strong, adata


def main():
    ap = argparse.ArgumentParser(description="BIRDccNEST trajectory inference (joint).")
    ap.add_argument("--h5ad", required=True, help="normalized AnnData (cells x genes) with obs labels")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-percent", type=float, default=0.25, help="starting top edge fraction kept")
    ap.add_argument("--resolution", type=float, default=1.0, help="Louvain resolution")
    ap.add_argument("--n-steps", type=int, default=1000, help="RegDiffusion training steps")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--label-col", default="clusterid", help="obs column of the reference (Lamian) clusters")
    ap.add_argument("--from-cache", action="store_true", help="reuse a cached pruned network in --outdir")
    ap.add_argument("--resolution-sweep", default=None,
                    help="comma list of resolutions; reports community counts + coherence, then exits")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sparse_adj, names, kept_frac, strong, adata = build_or_load_graph(args)

    # optional: sweep Louvain resolution to find granularity matching the reference
    if args.resolution_sweep:
        emb = np.asarray(adata.obsm["X_pca"]) if "X_pca" in adata.obsm else None
        rows = []
        for res in [float(x) for x in args.resolution_sweep.split(",")]:
            _, cdict, cl = louvain_communities(sparse_adj, names, res, args.seed)
            row = {"resolution": res, "n_communities": len(cl)}
            if emb is not None:
                labs = pd.Series(cdict).reindex(adata.obs_names).values
                row["DBI"], row["silhouette"] = coherence(emb, labs)
            rows.append(row)
            print(f"[sweep] res={res}: {len(cl)} communities"
                  + (f", DBI={row['DBI']:.3f}, sil={row['silhouette']:.3f}" if emb is not None else ""))
        pd.DataFrame(rows).to_csv(f"{args.outdir}/resolution_sweep.csv", index=False)
        print(f"[BIRDccNEST] sweep -> {args.outdir}/resolution_sweep.csv")
        return

    # 4) Louvain communities
    G, comm_dict, comms = louvain_communities(sparse_adj, names, args.resolution, args.seed)
    print(f"[BIRDccNEST] Louvain (res={args.resolution}): {len(comms)} communities")

    # 5) cluster flow network + 6) OMST
    flowG = make_cluster_flow_network(G, comm_dict)
    omst = extract_omst(flowG)
    print(f"[BIRDccNEST] flow network: {flowG.number_of_nodes()} nodes / {flowG.number_of_edges()} edges; "
          f"OMST: {omst.number_of_edges()} edges (DAG: {nx.is_directed_acyclic_graph(omst)})")

    # ---- assemble outputs ----
    obs = adata.obs
    cell_df = pd.DataFrame({"cell": list(comm_dict.keys()),
                            "community": list(comm_dict.values())}).set_index("cell")
    for col in [args.label_col, "condition", "pseudotime"]:
        if col in obs.columns:
            cell_df[col] = obs[col].reindex(cell_df.index)
    cell_df.reset_index().to_csv(f"{args.outdir}/cell_communities.csv", index=False)

    if args.label_col in obs.columns:
        overlap_grid(comm_dict, obs[args.label_col]).to_csv(f"{args.outdir}/community_overlap_{args.label_col}.csv")

    pd.DataFrame([(u, v, d["weight"]) for u, v, d in flowG.edges(data=True)],
                 columns=["u", "v", "weight"]).to_csv(f"{args.outdir}/cluster_flow_edges.csv", index=False)
    pd.DataFrame([(u, v, d["weight"]) for u, v, d in omst.edges(data=True)],
                 columns=["u", "v", "weight"]).to_csv(f"{args.outdir}/omst_edges.csv", index=False)

    # coherence: BIRDccNEST communities vs reference clusters, on the PCA embedding
    report = {"kept_fraction": kept_frac, "strongly_connected": strong,
              "n_communities": len(comms), "resolution": args.resolution,
              "n_steps": args.n_steps, "seed": args.seed}
    if "X_pca" in adata.obsm:
        emb = np.asarray(adata.obsm["X_pca"])
        comm_labels = pd.Series(comm_dict).reindex(adata.obs_names).values
        dbi_c, sil_c = coherence(emb, comm_labels)
        report["community_DBI"], report["community_silhouette"] = dbi_c, sil_c
        if args.label_col in obs.columns:
            dbi_r, sil_r = coherence(emb, obs[args.label_col].values)
            report["reference_DBI"], report["reference_silhouette"] = dbi_r, sil_r

    # orientation validation against joint pseudotime + root-anchored re-orientation
    if "pseudotime" in obs.columns:
        frac, dir_df, comm_pt = direction_concordance(omst, comm_dict, obs["pseudotime"].to_dict())
        dir_df.to_csv(f"{args.outdir}/omst_direction_check.csv", index=False)
        comm_pt.rename("mean_pseudotime").to_csv(f"{args.outdir}/community_mean_pseudotime.csv")
        report["omst_pt_direction_concordance"] = frac
        # biological root = origin (lowest mean pseudotime); re-orient the backbone from it
        root = int(comm_pt.idxmin())
        report["root_community_by_min_pt"] = root
        rooted = root_anchored_orientation(omst, root)
        pd.DataFrame([(u, v, d["weight"]) for u, v, d in rooted.edges(data=True)],
                     columns=["u", "v", "weight"]).to_csv(f"{args.outdir}/omst_rooted_edges.csv", index=False)
        rfrac, _, _ = direction_concordance(rooted, comm_dict, obs["pseudotime"].to_dict())
        report["rooted_pt_direction_concordance"] = rfrac

    pd.Series(report).to_csv(f"{args.outdir}/summary.csv")
    print("[BIRDccNEST] summary:")
    for k, v in report.items():
        print(f"    {k}: {v}")
    print(f"[BIRDccNEST] outputs -> {args.outdir}")


if __name__ == "__main__":
    main()
