"""
Clustering circularity control (Leng EC)
=========================================

The double-positive enrichment in a single k-means cluster could in principle be
circular: if RELN and RORB are themselves in the highly-variable-gene set that
defines the clusters, cells are partly clustered ON the markers being counted.

This re-runs the identical pipeline with RELN and RORB REMOVED from the HVG set,
so the subtype is defined only by its other genes. Double-positive labels are
taken from the raw counts before any feature selection.

Reported in manuscript section 3.2 and Supplementary Fig. S3e.

Inputs : data/leng2021_ec_full.h5ad
Outputs: clustering_circularity.csv (printed to stdout)
Method : scanpy preprocessing + sklearn KMeans + scipy hypergeometric test

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown.
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import DATA, require
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from sklearn.cluster import KMeans
from scipy.stats import hypergeom

N_HVG, N_PC, K, SEED = 2000, 30, 8, 0


def load_excitatory(path=None):
    if path is None:
        path = require("leng2021_ec_full.h5ad")
    a = ad.read_h5ad(path)
    mask = a.obs["clusterCellType"].astype(str).str.contains("Exc", case=False, na=False)
    L = a[mask].copy()
    L.var_names = L.var["feature_name"].astype(str).values
    L.var_names_make_unique()
    return L


def double_positive(L):
    """DP labels from RAW counts -- fixed before any feature selection."""
    def col(sym):
        i = list(L.var_names).index(sym)
        x = L.X[:, i]
        return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()
    return (col("RELN") > 0) & (col("RORB") > 0)


def cluster_and_enrich(L, dp, exclude_markers):
    B = L.copy()
    sc.pp.normalize_total(B, target_sum=1e4)
    sc.pp.log1p(B)
    sc.pp.highly_variable_genes(B, n_top_genes=N_HVG)
    hv = B.var["highly_variable"].copy()
    removed = []
    if exclude_markers:
        for g in ("RELN", "RORB"):
            if g in hv.index and bool(hv[g]):
                hv[g] = False
                removed.append(g)
    B = B[:, hv.values].copy()
    sc.pp.scale(B, max_value=10)
    sc.tl.pca(B, n_comps=N_PC, svd_solver="arpack", random_state=SEED)
    lab = KMeans(n_clusters=K, random_state=SEED, n_init=10).fit(B.obsm["X_pca"]).labels_

    N, tot = len(lab), int(dp.sum())
    best = None
    for c in range(K):
        m = lab == c
        n = int(m.sum())
        if n == 0:
            continue
        k = int(dp[m].sum())
        fold = (k / n) / (tot / N)
        pv = hypergeom.sf(k - 1, N, tot, n)
        if best is None or fold > best[0]:
            best = (fold, pv, n, k)
    return best, removed, int(hv.sum())


def main():
    L = load_excitatory()
    dp = double_positive(L)
    print(f"EC excitatory nuclei: {L.shape[0]}  |  double-positive: {int(dp.sum())}\n")

    rows = []
    for exclude, label in ((False, "markers IN HVG set (as published)"),
                           (True, "RELN/RORB REMOVED from HVG set")):
        (fold, pv, n, k), removed, nhv = cluster_and_enrich(L, dp, exclude)
        rows.append({"configuration": label, "fold_enrichment": round(fold, 3),
                     "hypergeom_p": pv, "cluster_n": n, "cluster_dp": k,
                     "n_hvg": nhv, "removed": ",".join(removed) or "-"})
        print(f"{label:38s} fold={fold:.2f}  p={pv:.2e}  (cluster n={n}, DP={k}, HVGs={nhv})")

    print("\nIf the enrichment is essentially unchanged, the double-positive population is\n"
          "not an artifact of clustering on the markers being counted.")
    pd.DataFrame(rows).to_csv("clustering_circularity.csv", index=False)
    print("\nwrote clustering_circularity.csv")


if __name__ == "__main__":
    main()
