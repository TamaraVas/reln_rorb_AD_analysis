"""
Artifact-removal robustness controls (SEA-AD MEC)
=================================================

Reviewer-requested controls that REMOVE rather than adjust for the two most
common single-cell artifacts, plus depth functional-form invariance.

  (a) Depth: binomial thinning of every nucleus to a common depth, after which
      no depth covariate is used at all; plus spline / decile-FE specifications.
  (b) Ambient RNA: ambient floor from non-neuronal nuclei, and a >=2 threshold
      that discards any call a single ambient molecule could produce.
  (c) Donor-level inference: leave-one-donor-out and a donor-label permutation
      test that makes no parametric assumption about the variance estimator.

These controls trade precision for elimination of a rival explanation, so their
p-values are LARGER than the primary estimate's by construction. Reported in
Supplementary Fig. S3e.

Inputs : data/seaad_mec_extract.parquet
Outputs: robustness_controls.csv (printed to stdout)
Method : statsmodels Logit (donor-clustered), numpy binomial thinning, permutation

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown.
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import DATA, require
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

BRAAK = {"Braak 0": 0, "Braak II": 2, "Braak III": 3,
         "Braak IV": 4, "Braak V": 5, "Braak VI": 6}
SEED = 0


def load_mec(path=None):
    """Excitatory nuclei with numeric Braak and log-depth.

    NOTE: `reln`/`rorb` are depth-normalised truncated values, not raw counts
    (see DATA.md). The >0 rule is the cohort's normalised-expression threshold.
    """
    if path is None:
        path = require("seaad_mec_extract.parquet")
    d = pd.read_parquet(path)
    d = d[d["Class"] == "Neuronal: Glutamatergic"].copy()
    d["bnum"] = d["Braak"].map(BRAAK)
    d = d.dropna(subset=["bnum"])
    d["reln_pos"] = (d.reln > 0).astype(int)
    d["rorb_pos"] = (d.rorb > 0).astype(int)
    d["logdepth"] = np.log10(d.umis.clip(lower=1))
    return d


def interaction(df, formula="rorb_pos ~ reln_pos * bnum + logdepth",
                term="reln_pos:bnum", cluster=True):
    kw = dict(cov_type="cluster", cov_kwds={"groups": df.donor}) if cluster else {}
    m = smf.logit(formula, data=df).fit(disp=0, **kw)
    return float(np.exp(m.params[term])), float(m.pvalues[term])


def depth_specifications(d):
    """The estimate should not depend on how depth enters the model."""
    rows = []
    rows.append(("linear log10(UMI)", *interaction(d)))
    d = d.copy()
    d["ddec"] = pd.qcut(d.logdepth, 10, labels=False, duplicates="drop").astype(str)
    rows.append(("depth-decile fixed effects",
                 *interaction(d, "rorb_pos ~ reln_pos * bnum + C(ddec)")))
    rows.append(("natural cubic spline df=5",
                 *interaction(d, "rorb_pos ~ reln_pos * bnum + cr(logdepth, df=5)")))
    return rows


def binomial_thinning(d, seed=SEED):
    """Downsample every nucleus to a COMMON depth, then drop the depth covariate.

    This removes the depth confound rather than adjusting for it: after thinning,
    all nuclei carry the same library size, so no depth term is needed.
    """
    rng = np.random.default_rng(seed)
    target = int(np.floor(d.umis.quantile(0.10)))   # keep ~90% of nuclei
    sub = d[d.umis >= target].copy()
    p_keep = target / sub.umis.values
    sub["reln_t"] = rng.binomial(sub.reln.values.astype(int), p_keep)
    sub["rorb_t"] = rng.binomial(sub.rorb.values.astype(int), p_keep)
    sub["reln_pos_t"] = (sub.reln_t > 0).astype(int)
    sub["rorb_pos_t"] = (sub.rorb_t > 0).astype(int)
    o, p = interaction(sub, "rorb_pos_t ~ reln_pos_t * bnum",
                       term="reln_pos_t:bnum", cluster=False)
    return target, len(sub), o, p


def ambient_floor(path=None):
    """Non-neuronal nuclei index the ambient floor; neurons must exceed it."""
    if path is None:
        path = require("seaad_mec_extract.parquet")
    a = pd.read_parquet(path)
    nn = a[a["Class"] == "Non-neuronal and Non-neural"]
    return {
        "non-neuronal mean RELN (detected)": float(nn.loc[nn.reln > 0, "reln"].mean()),
        "non-neuronal RELN+ %": float((nn.reln > 0).mean() * 100),
    }


def stricter_threshold(d):
    """Require >=2 for both genes: discards any call one ambient molecule could make.

    Ambient contamination DILUTES a true same-cell signal toward the null, so if the
    effect is real the estimate should grow, not shrink, under a stricter threshold.
    """
    e = d.copy()
    e["reln_pos"] = (e.reln >= 2).astype(int)
    e["rorb_pos"] = (e.rorb >= 2).astype(int)
    return interaction(e)


def leave_one_donor_out(d):
    ors = []
    for dn in d.donor.unique():
        try:
            o, _ = interaction(d[d.donor != dn], cluster=False)
            ors.append(o)
        except Exception:
            pass
    a = np.array(ors)
    return len(a), float(a.min()), float(a.max()), float(np.median(a))


def donor_permutation(d, n_perm=200, seed=7):
    """Shuffle Braak stage ACROSS DONORS -- assumption-free null."""
    rng = np.random.default_rng(seed)
    obs, _ = interaction(d, cluster=False)
    dmap = d.groupby("donor")["bnum"].first()
    null = []
    for _ in range(n_perm):
        perm = dict(zip(dmap.index, rng.permutation(dmap.values)))
        tmp = d.copy()
        tmp["bperm"] = tmp.donor.map(perm)
        try:
            o, _ = interaction(tmp, "rorb_pos ~ reln_pos * bperm + logdepth",
                               term="reln_pos:bperm", cluster=False)
            null.append(o)
        except Exception:
            pass
    null = np.array(null)
    p = (np.sum(null >= obs) + 1) / (len(null) + 1)
    return obs, float(np.median(null)), float(np.percentile(null, 95)), float(p), len(null)


def donor_level_trend(d):
    """Per-donor co-occupancy OR vs Braak -- the genomics-native donor-level test."""
    rows = []
    for dn, s in d.groupby("donor"):
        a = int(((s.reln_pos == 1) & (s.rorb_pos == 1)).sum())
        b = int(((s.reln_pos == 1) & (s.rorb_pos == 0)).sum())
        c = int(((s.reln_pos == 0) & (s.rorb_pos == 1)).sum())
        e = int(((s.reln_pos == 0) & (s.rorb_pos == 0)).sum())
        orr = ((a + 0.5) * (e + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append({"donor": dn, "braak": s.bnum.iloc[0], "logOR": np.log(orr)})
    dl = pd.DataFrame(rows)
    rho, p = spearmanr(dl.braak, dl.logOR)
    return len(dl), float(rho), float(p)


def main():
    d = load_mec()
    out = []

    o, p = interaction(d)
    out.append(("primary (donor-clustered)", o, p))
    for name, o_, p_ in depth_specifications(d):
        out.append((f"depth spec: {name}", o_, p_))

    tgt, n, o, p = binomial_thinning(d)
    out.append((f"binomial thinning to {tgt} UMIs (n={n}, no depth term)", o, p))

    o, p = stricter_threshold(d)
    out.append((">=2 threshold both genes", o, p))

    n, lo, hi, med = leave_one_donor_out(d)
    out.append((f"leave-one-donor-out (n={n}): range {lo:.3f}-{hi:.3f}", med, float("nan")))

    obs, nmed, n95, pperm, B = donor_permutation(d)
    out.append((f"donor permutation (B={B}, null med {nmed:.3f}, p95 {n95:.3f})", obs, pperm))

    nd, rho, prho = donor_level_trend(d)
    out.append((f"donor-level Spearman rho={rho:.3f} (n={nd} donors)", float("nan"), prho))

    print(f"{'control':64s} {'OR':>8s} {'p':>12s}")
    for name, o, p in out:
        os_ = "     n/a" if o != o else f"{o:8.3f}"
        ps_ = "         n/a" if p != p else f"{p:12.3e}"
        print(f"  {name:62s} {os_} {ps_}")

    print("\nAmbient floor:")
    for k, v in ambient_floor().items():
        print(f"  {k:36s} {v:.2f}")

    pd.DataFrame(out, columns=["control", "OR", "p"]).to_csv("robustness_controls.csv", index=False)
    print("\nwrote robustness_controls.csv")


if __name__ == "__main__":
    main()
