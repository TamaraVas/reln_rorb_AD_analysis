"""
Within-donor regional control: MEC vs hippocampus (MERFISH)
===========================================================

Both SEA-AD spatial sections come from ONE donor (H24.30.005), confirmed by the
data generators. That makes the hippocampal section a WITHIN-DONOR regional
control rather than an independent replication: same brain, same processing,
same segmentation, different region -- so a regional difference cannot be a
donor-level confound.

PROVENANCE NOTE -- the generators also confirmed these sections were profiled to
localize transcriptionally defined cell types in a NEUROTYPICAL context, and that
this donor is unlikely to carry significant pathological burden. MEC/HIP spatial
data across levels of pathological burden have not yet been generated. The spatial
analysis therefore establishes ANATOMY (the double-positive population exists in
intact tissue, is glutamatergic, spatially clustered and entorhinal-specific), not
the disease gradient, which rests entirely on the staged snRNA-seq cohorts.

Also note the spatial and single-nucleus ORs are not on a common measurement
scale: MERFISH records ~16 RELN transcripts in a positive neuron vs ~1 in
dissociated nuclei. Reported as Figure 9F.

Inputs : data/seaad_merfish_mec.h5ad, data/seaad_merfish_hpf.h5ad
Outputs: within_donor_regional.csv (printed to stdout)
Method : statsmodels Logit (depth-controlled) + label-shuffle kNN permutation

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
from sklearn.neighbors import NearestNeighbors

DONOR = "H24.30.005"   # both sections, confirmed by the data generators
K_NEIGHBOURS = 10
N_PERM = 1000
SEED = 0


def load_section(path):
    """Glutamatergic neurons with RELN/RORB counts and spatial coordinates."""
    a = ad.read_h5ad(path)
    cls = a.obs.get("Class_scANVI", a.obs.get("Class"))
    keep = cls.astype(str).str.contains("Glut", case=False, na=False)
    s = a[keep.values].copy()

    def gene(sym):
        i = list(s.var_names).index(sym)
        x = s.X[:, i]
        return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

    total = np.asarray(s.X.sum(axis=1)).ravel() if sp.issparse(s.X) else s.X.sum(axis=1)
    xy = s.obsm["spatial"] if "spatial" in s.obsm else s.obs[["x", "y"]].values
    return pd.DataFrame({
        "reln": gene("RELN"), "rorb": gene("RORB"), "total": total,
        "x": xy[:, 0], "y": xy[:, 1]})


def same_cell_or(d):
    d = d.copy()
    d["reln_pos"] = (d.reln > 0).astype(int)
    d["rorb_pos"] = (d.rorb > 0).astype(int)
    d["logdepth"] = np.log10(np.clip(d.total, 1, None))
    m = smf.logit("rorb_pos ~ reln_pos + logdepth", data=d).fit(disp=0)
    return float(np.exp(m.params["reln_pos"])), float(m.pvalues["reln_pos"])


def dp_spatial_clustering(d, k=K_NEIGHBOURS, n_perm=N_PERM, seed=SEED):
    """Do double-positive cells neighbour each other more than chance?

    The permutation shuffles DP labels while holding positions and the DP count
    fixed, so the null preserves both geometry and marginal prevalence.
    """
    rng = np.random.default_rng(seed)
    dp = ((d.reln > 0) & (d.rorb > 0)).values
    if dp.sum() < 10:
        return np.nan, int(dp.sum())
    xy = d[["x", "y"]].values
    nn = NearestNeighbors(n_neighbors=k + 1).fit(xy)
    idx = nn.kneighbors(xy, return_distance=False)[:, 1:]

    obs = dp[idx][dp].mean()
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(dp)
        null[i] = perm[idx][perm].mean()
    z = (obs - null.mean()) / null.std()
    return float(z), int(dp.sum())


def main():
    rows = []
    mec_path, hpf_path = require("seaad_merfish_mec.h5ad", "seaad_merfish_hpf.h5ad")
    for region, path in (("MEC (entorhinal)", mec_path),
                         ("HPF (hippocampus)", hpf_path)):
        d = load_section(path)
        o, p = same_cell_or(d)
        z, ndp = dp_spatial_clustering(d)
        rows.append({"region": region, "donor": DONOR, "n_glut": len(d),
                     "n_double_positive": ndp,
                     "same_cell_OR_depth_controlled": round(o, 3), "p": p,
                     "dp_clustering_z": round(z, 1)})
        print(f"{region:20s} n={len(d):6d}  DP={ndp:5d}  "
              f"depth-controlled OR={o:.3f} (p={p:.2e})  clustering z={z:.1f}")

    print(f"\nBoth sections are from donor {DONOR}, so the regional difference is a\n"
          "within-donor contrast: it cannot be explained by donor-level confounding.\n"
          "Tissue is neurotypical by design -- this establishes anatomical specificity,\n"
          "not a disease effect.")
    pd.DataFrame(rows).to_csv("within_donor_regional.csv", index=False)
    print("\nwrote within_donor_regional.csv")


if __name__ == "__main__":
    main()
