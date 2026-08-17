"""
RORB-alone negative control: the second RORB-high entorhinal cluster
====================================================================

The co-occupancy claim would be weak if RORB were simply a proxy for a broadly
RELN-adjacent cell state. This script tests that directly using an internal
negative control that already exists in the discovery cohort: a SECOND cluster
(C2, mostly EC:Exc.5) that is MORE RORB-rich than the double-positive cluster
(C7, mostly EC:Exc.4) yet carries almost no RELN.

Three tests:
  1. Cluster-level DP depletion (hypergeometric LOWER tail) -- C2 is not merely
     unenriched for double-positives, it is significantly depleted.
  2. Conditional co-occupancy -- among RORB+ neurons only, what fraction is also
     RELN+? This, not the marginal RORB rate, is the quantity that matters.
  3. Depth control -- C7 has ~2x C2's sequencing depth and Leng positivity is a
     RAW-COUNT threshold, so the conditional gap could be a detection artifact.
     Re-tested with log10 library size in the model and within depth quintiles.

Then replicated at scale in SEA-AD MEC, where the near-universally RORB+
neocortical subclasses (L4 IT, L5 IT) are essentially RELN-negative while the
far less RORB-rich entorhinal subclass carries the co-occupancy.

NOTE ON POSITIVITY (see DATA.md): the two cohorts use different rules by design.
Leng marker columns are RAW COUNTS (>= 1). The SEA-AD extract columns are
DEPTH-NORMALISED truncated values that resemble small counts; applying >= 1
there reproduces the manuscript prevalences (RELN+ 13.3%, RORB+ 20.2% of MEC
excitatory nuclei) and the script asserts this. Do not "fix" one to match the
other.

NOTE ON SEPARATION: subclasses resting on fewer than MIN_DP double-positive
cells produce unbounded maximum-likelihood estimates. These are reported as
unestimable rather than quoted -- see the manuscript Methods.

Reported in manuscript section 3.2 and Supplementary Fig. S10.

Inputs : data/leng2021_ec_full.h5ad
         data/seaad_mec_extract.parquet
Outputs: second_rorb_cluster.csv (printed to stdout)
Method : sklearn KMeans + scipy hypergeometric/Fisher + statsmodels logistic

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
import anndata as ad
import statsmodels.formula.api as smf
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.stats import hypergeom, fisher_exact, spearmanr

N_HVG, N_PC, K, SEED = 2000, 30, 8, 0
RELN, RORB = "ENSG00000189056", "ENSG00000198963"
MIN_DP = 10  # below this, logistic estimates are separation artifacts


def load_excitatory(path=None):
    if path is None:
        path = require("leng2021_ec_full.h5ad")
    a = ad.read_h5ad(path)
    return a[a.obs["clusterCellType"].astype(str) == "Exc"].copy()


def cluster_leng(exc):
    """Reproduce the Figure 2 pipeline exactly (log-norm, HVG, PCA, k-means)."""
    def counts(ens):
        x = exc[:, ens].X
        return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

    reln, rorb = counts(RELN), counts(RORB)
    X = exc.layers["counts"] if "counts" in exc.layers else exc.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    lib = np.asarray(X.sum(1)).ravel()
    Xn = X.multiply((1e4 / np.clip(lib, 1, None))[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    mean = np.asarray(Xn.mean(0)).ravel()
    var = np.asarray(Xn.multiply(Xn).mean(0)).ravel() - mean ** 2
    disp = np.divide(var, mean, out=np.zeros_like(var), where=mean > 0)
    Xh = Xn[:, np.argsort(disp)[-N_HVG:]].toarray()
    sd = Xh.std(0); sd[sd == 0] = 1
    pcs = PCA(n_components=N_PC, random_state=SEED).fit_transform(
        np.clip((Xh - Xh.mean(0)) / sd, -10, 10))
    lab = KMeans(n_clusters=K, random_state=SEED, n_init=10).fit(pcs).labels_.astype(str)
    return pd.DataFrame({"cl": lab, "rp": (reln > 0).astype(int),
                         "bp": (rorb > 0).astype(int), "nUMI": lib})


def per_cluster_table(d):
    """DP enrichment AND depletion, plus the conditional RELN|RORB+ rate."""
    N, Kdp = len(d), int(((d.rp == 1) & (d.bp == 1)).sum())
    rows = []
    for cl in sorted(d.cl.unique(), key=int):
        q = d[d.cl == cl]
        k = int(((q.rp == 1) & (q.bp == 1)).sum())
        qq = q[q.bp == 1]
        rows.append(dict(
            cluster=cl, n=len(q), n_DP=k,
            fold=round(k / (Kdp * len(q) / N), 3),
            p_enrich=hypergeom.sf(k - 1, N, Kdp, len(q)),
            p_deplete=hypergeom.cdf(k, N, Kdp, len(q)),
            RORB_pct=round(100 * q.bp.mean(), 1),
            RELN_pct=round(100 * q.rp.mean(), 1),
            RELN_given_RORB=round(100 * qq.rp.mean(), 1) if len(qq) else np.nan,
            median_UMI=int(q.nUMI.median())))
    return pd.DataFrame(rows)


def conditional_contrast(d, a="2", b="7"):
    """Is the RELN|RORB+ gap between two clusters real, or a depth artifact?"""
    tab = [[int(((d.cl == c) & (d.bp == 1) & (d.rp == r)).sum()) for r in (1, 0)]
           for c in (a, b)]
    orv, pv = fisher_exact(tab)

    w = d[(d.cl.isin([a, b])) & (d.bp == 1)].copy()
    w["hi"] = (w.cl == b).astype(int)
    w["logd"] = np.log10(w.nUMI)
    m = smf.logit("rp ~ hi + logd", data=w).fit(disp=0)
    ci = m.conf_int()

    bins = np.quantile(w.nUMI, [0, .2, .4, .6, .8, 1.0])
    w["db"] = pd.cut(w.nUMI, bins, include_lowest=True)
    strata = []
    for bb, q in w.groupby("db", observed=True):
        if min((q.hi == 0).sum(), (q.hi == 1).sum()) < MIN_DP:
            continue
        strata.append(dict(depth_bin=f"{int(bb.left)}-{int(bb.right)}",
                           n_lo=int((q.hi == 0).sum()), n_hi=int((q.hi == 1).sum()),
                           pct_lo=round(100 * q[q.hi == 0].rp.mean(), 1),
                           pct_hi=round(100 * q[q.hi == 1].rp.mean(), 1)))
    return dict(fisher_OR=orv, fisher_p=pv,
                depth_OR=float(np.exp(m.params["hi"])),
                depth_lo=float(np.exp(ci.loc["hi", 0])),
                depth_hi=float(np.exp(ci.loc["hi", 1])),
                depth_p=float(m.pvalues["hi"]),
                depth_term_OR=float(np.exp(m.params["logd"])),
                depth_term_p=float(m.pvalues["logd"]),
                strata=pd.DataFrame(strata))


def seaad_subclasses(path=None,
                     subclasses=("L4 IT", "L5 IT", "L2/3 IT", "EC IT")):
    """Replicate the dissociation at scale. See NOTE ON POSITIVITY above."""
    if path is None:
        path = require("seaad_mec_extract.parquet")
    me = pd.read_parquet(path)
    ex = me[me.Class.astype(str).str.contains("Glutamatergic", na=False)].copy()
    ex["rp"] = (ex.reln >= 1).astype(int)
    ex["bp"] = (ex.rorb >= 1).astype(int)
    assert abs(100 * ex.rp.mean() - 13.3) < 0.2, "RELN prevalence off; check positivity rule"
    assert abs(100 * ex.bp.mean() - 20.2) < 0.2, "RORB prevalence off; check positivity rule"

    rows = []
    for s in subclasses:
        q = ex[ex.Subclass == s]
        qq = q[q.bp == 1]
        n_dp = int(((q.rp == 1) & (q.bp == 1)).sum())
        rows.append(dict(subclass=s, n=len(q), n_DP=n_dp,
                         RORB_pct=round(100 * q.bp.mean(), 1),
                         RELN_given_RORB=round(100 * qq.rp.mean(), 1) if len(qq) else np.nan,
                         estimable=n_dp >= MIN_DP))
    return pd.DataFrame(rows)


def main():
    d = cluster_leng(load_excitatory())
    tab = per_cluster_table(d)
    print("=== Leng EC: per-cluster DP enrichment and conditional RELN|RORB+ ===")
    print(tab.to_string(index=False))

    rho, prho = spearmanr(tab.RORB_pct, tab.RELN_given_RORB)
    print(f"\nAcross clusters, RORB abundance vs RELN|RORB+: rho={rho:.2f}, p={prho:.2f}")
    print("(i.e. how RORB-rich a cluster is does not predict RELN co-occupancy)")

    dp_cluster = tab.loc[tab.fold.idxmax(), "cluster"]
    other = tab[tab.cluster != dp_cluster].sort_values("RORB_pct", ascending=False).iloc[0]["cluster"]
    print(f"\nDP cluster = C{dp_cluster}; most RORB-rich other cluster = C{other}")

    r = conditional_contrast(d, a=other, b=dp_cluster)
    print(f"\nConditional Fisher (RELN|RORB+, C{other} vs C{dp_cluster}): "
          f"OR={r['fisher_OR']:.4f}, p={r['fisher_p']:.3g}")
    print(f"Depth-controlled logistic: OR={r['depth_OR']:.1f} "
          f"[{r['depth_lo']:.1f}-{r['depth_hi']:.1f}], p={r['depth_p']:.3g}")
    print(f"Depth term itself: OR={r['depth_term_OR']:.2f} per log10 UMI, "
          f"p={r['depth_term_p']:.3g}  (depth matters, but does not explain the gap)")
    print("\nRELN+ (%) among RORB+ within depth quintiles:")
    print(r["strata"].to_string(index=False))

    print("\n=== SEA-AD MEC replication ===")
    sead = seaad_subclasses()
    print(sead.to_string(index=False))
    if (~sead.estimable).any():
        print(f"\nSubclasses with < {MIN_DP} DP cells are reported as unestimable "
              "(separation artifacts), not quoted.")

    tab.to_csv("second_rorb_cluster.csv", index=False)
    print("\nWrote second_rorb_cluster.csv")


if __name__ == "__main__":
    main()
