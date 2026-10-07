"""
Continuous pathology axes: does the convergence track tau or amyloid?
=====================================================================

Braak and Thal are ordinal summaries, and in the ordinal model both tau and amyloid
stage remain associated with the convergence, so neither can be excluded. The 2026
release publishes donor-level pseudo-progression scores that model tau and amyloid
burden as SEPARATE CONTINUOUS axes, at brain-wide and region-local scope. That
permits a mutual adjustment the ordinal variables cannot support.

Result, and it cuts both ways -- both directions are reported:
  brain-wide : tau survives adjustment for amyloid (OR 1.44/SD, p = 0.011);
               amyloid does not survive adjustment for tau (p = 0.13)
  MEC-local  : the ordering INVERTS (local amyloid p = 0.020, local tau p = 0.152)

The local tau axis is by far the most range-restricted of the four (SD 0.162 vs
0.246 for local amyloid on a 0-1 scale), which is what would be expected in the
region where tau originates and is near-saturating at cohort entry -- but that is an
explanation, not a demonstration, so the reversal is reported rather than dismissed.
Axes are standardised before mutual adjustment so the two are on a common scale.

Inputs : data/seaad2026_MEC_final.parquet, data/Global_and_Local_CPS.csv  (script 25)
Outputs: cps_tau_vs_amyloid.csv, cps_axis_spread.csv
Method : depth-controlled logistic regression, donor-clustered SEs

Part of the RELN-RORB EC convergence analysis. See README.md.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import require

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

AXES = ["CPS_Global_pTau", "CPS_Global_ABeta", "CPS_Local_pTau", "CPS_Local_ABeta"]
REGION = "MEC"


def main():
    mec_p, cps_p = require("seaad2026_MEC_final.parquet", "Global_and_Local_CPS.csv")
    d = pd.read_parquet(mec_p)
    d = d[d["Class"].astype(str).str.contains("Glutamatergic", na=False)].copy()
    d["rp"] = (d["RELN_norm"] >= 1).astype(int)
    d["bp"] = (d["RORB_norm"] >= 1).astype(int)
    d["logd"] = np.log10(pd.to_numeric(d["Number of UMIs"], errors="coerce").clip(lower=1))

    cps = pd.read_csv(cps_p)
    cps = cps[cps["Brain Region"] == REGION].set_index("Donor ID")
    for c in AXES:
        d[c] = d["Donor ID"].astype(str).map(cps[c])
    q = d.dropna(subset=AXES + ["logd"]).copy()
    print("n = %d nuclei, %d donors" % (len(q), q["Donor ID"].nunique()), flush=True)
    g = q["Donor ID"].astype("category").cat.codes.values

    rows = []
    for c in AXES:
        m = smf.logit("bp ~ rp * %s + logd" % c, data=q).fit(
            disp=0, cov_type="cluster", cov_kwds={"groups": g})
        k, ci = "rp:%s" % c, m.conf_int()
        rows.append(dict(axis=c, model="single axis", OR=round(float(np.exp(m.params[k])), 3),
                         ci="%.2f-%.2f" % (np.exp(ci.loc[k, 0]), np.exp(ci.loc[k, 1])),
                         p=float(m.pvalues[k])))

    # Standardise so the two axes are comparable per SD, then adjust mutually.
    # DONOR-LEVEL moments, not cell-level: CPS is a donor property, so "per SD" must
    # mean per SD of the donor distribution. Standardising on the nucleus-weighted
    # distribution lets high-yield donors dominate the SD and shifts the per-SD OR
    # (global tau 1.438 -> 1.427, global amyloid 1.287 -> 1.303). p-values are
    # unaffected either way, since rescaling a predictor cannot change them.
    dl_z = q.groupby("Donor ID")[AXES].first()
    for c in AXES:
        q[c + "_z"] = (q[c] - dl_z[c].mean()) / dl_z[c].std()
    for scope in ("Global", "Local"):
        tau, amy = "CPS_%s_pTau_z" % scope, "CPS_%s_ABeta_z" % scope
        m = smf.logit("bp ~ rp*%s + rp*%s + logd" % (tau, amy), data=q).fit(
            disp=0, cov_type="cluster", cov_kwds={"groups": g})
        ci = m.conf_int()
        for k in ("rp:" + tau, "rp:" + amy):
            rows.append(dict(axis=k.replace("rp:", "").replace("_z", "") + " (per SD)",
                             model="mutually adjusted", OR=round(float(np.exp(m.params[k])), 3),
                             ci="%.2f-%.2f" % (np.exp(ci.loc[k, 0]), np.exp(ci.loc[k, 1])),
                             p=float(m.pvalues[k])))
    R = pd.DataFrame(rows)
    R["n_nuclei"], R["n_donors"] = len(q), q["Donor ID"].nunique()
    R.to_csv("cps_tau_vs_amyloid.csv", index=False)
    print(R.to_string(index=False), flush=True)

    dl = q.groupby("Donor ID")[AXES].first()
    sp = pd.DataFrame([dict(axis=c, sd=round(float(dl[c].std()), 3),
                            q25=round(float(dl[c].quantile(.25)), 3),
                            median=round(float(dl[c].median()), 3),
                            note="0-1 scale; smaller SD = more range-restricted")
                       for c in AXES])
    sp.to_csv("cps_axis_spread.csv", index=False)
    print("\n=== axis spread across donors (why the local reversal is not interpreted) ===", flush=True)
    print(sp.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
