"""
Regional and subclass localisation of the RELN-RORB convergence (ten regions)
=============================================================================

Three questions the four-region version of this analysis could not answer:

  1. Is the tau-graded convergence entorhinal, or does it appear wherever RORB is
     abundant? Ten regions, run in TWO donor cohorts so the difference is visible
     rather than assumed:
        "all"       - every donor profiled in that region
        "no_severe" - Severely Affected Donor == 'N' (the release's own flag)

  2. Region-level models are composition-sensitive. Two regions return a
     significant interaction under both cohorts (MEC and MTG), so the interaction
     is refit WITHIN subclass. This is the decisive test: MEC's double-positives
     are 97.7% entorhinal-IT (24.19% of RORB+ nuclei are RELN+), MTG's are 84.6%
     L4 IT (0.44%) -- the RORB-high/RELN-negative population this study uses as a
     negative control. A 55-fold difference, and MTG's within-subclass interaction
     is not significant. The MTG region-level signal is a detection-floor drift in
     a very large RORB+ population, not a weak entorhinal convergence.

  3. Do MEC and MTG reflect one donor-level process? MEC, MTG and PFC were profiled
     in the SAME donors, so the per-donor conditional rate can be correlated across
     regions directly. It is not correlated (rho = -0.16, n = 74), and the 95% CI
     excludes rho > 0.07 -- a shared moderate donor-level process is ruled out.

GUARD: a region or region x cohort cell resting on fewer than MIN_DP double-positive
nuclei is reported UNESTIMABLE, not as an odds ratio. Maximum-likelihood separation
inflates such estimates without bound -- hippocampus has 6 DP nuclei in 119,096 and
returns OR 8.3 with a "significant" negative interaction if left unguarded.

Inputs : data/seaad2026_<REGION>_final.parquet  (script 25)
Outputs: regional_two_cohorts.csv, subclass_interaction.csv,
         regional_braak_trajectories.csv, cross_region_donor_correlation.csv
Method : depth-controlled logistic regression, donor-clustered SEs; Spearman
         and Braak-partialled Spearman across donors

Part of the RELN-RORB EC convergence analysis. See README.md.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import DATA, require

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact, norm, spearmanr

REGIONS = ["MEC", "LEC", "HIP", "ITG", "MTG", "STG", "FI", "AnG", "PFC", "V1C"]
# Braak I is absent from this cohort in all ten regions; the map has no entry for it
# so a future object containing stage I fails the assertion below rather than silently
# collapsing into an adjacent stage.
BRAAK = {"Braak 0": 0, "Braak II": 2, "Braak III": 3, "Braak IV": 4, "Braak V": 5, "Braak VI": 6}
MIN_DP = 10        # below this, report UNESTIMABLE (separation)
MIN_RORB = 50      # per-donor RORB+ floor for the donor-level conditional rate
SUBCLASS_MIN_N = 2000


def load(region):
    p = require("seaad2026_%s_final.parquet" % region)
    d = pd.read_parquet(p)
    d = d[d["Class"].astype(str).str.contains("Glutamatergic", na=False)].copy()
    stages = set(d["Braak"].astype(str).unique())
    unknown = stages - set(BRAAK)
    assert not unknown, ("unmapped Braak levels in %s: %s -- if this includes "
                         "'Braak I', revisit every per-stage table" % (region, unknown))
    d["rp"] = (d["RELN_norm"] >= 1).astype(int)     # manuscript SEA-AD positivity rule
    d["bp"] = (d["RORB_norm"] >= 1).astype(int)
    d["bk"] = d["Braak"].astype(str).map(BRAAK)
    d["logd"] = np.log10(pd.to_numeric(d["Number of UMIs"], errors="coerce").clip(lower=1))
    return d.dropna(subset=["bk", "logd"])


def fit(d):
    """Baseline co-occupancy (Fisher) and RELN x Braak interaction (RORB outcome)."""
    dp = int(((d.rp == 1) & (d.bp == 1)).sum())
    row = dict(n=len(d), donors=d["Donor ID"].nunique(), DP_cells=dp,
               RELN_pct=round(100 * d.rp.mean(), 3), RORB_pct=round(100 * d.bp.mean(), 3),
               DP_pct=round(100 * dp / max(len(d), 1), 4), estimable=dp >= MIN_DP)
    if dp < MIN_DP:
        row.update(co_OR="UNESTIMABLE (DP<%d)" % MIN_DP, co_p="", int_OR="UNESTIMABLE (DP<%d)" % MIN_DP, int_p="")
        return row
    tab = [[dp, int(((d.rp == 1) & (d.bp == 0)).sum())],
           [int(((d.rp == 0) & (d.bp == 1)).sum()), int(((d.rp == 0) & (d.bp == 0)).sum())]]
    o, p = fisher_exact(tab)
    g = d["Donor ID"].astype("category").cat.codes.values
    m = smf.logit("bp ~ rp * bk + logd", data=d).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": g})
    row.update(co_OR=round(float(o), 3), co_p="%.3g" % p,
               int_OR=round(float(np.exp(m.params["rp:bk"])), 3), int_p="%.3g" % m.pvalues["rp:bk"])
    return row


def main():
    frames, main_rows, traj = {}, [], []
    for rg in REGIONS:
        d = load(rg)
        frames[rg] = d
        for cohort in ("all", "no_severe"):
            q = d if cohort == "all" else d[d["Severely Affected Donor"].astype(str) == "N"]
            main_rows.append(dict(region=rg, cohort=cohort, **fit(q)))
        for s, q in d.groupby("bk"):
            traj.append(dict(region=rg, braak=int(s), n=len(q), donors=q["Donor ID"].nunique(),
                             DP_pct=round(100 * ((q.rp == 1) & (q.bp == 1)).mean(), 4),
                             n_RORBpos=int(q.bp.sum()),
                             RELN_given_RORB=(round(100 * q[q.bp == 1].rp.mean(), 3)
                                              if q.bp.sum() else np.nan)))
    M = pd.DataFrame(main_rows)
    M.to_csv("regional_two_cohorts.csv", index=False)
    pd.DataFrame(traj).to_csv("regional_braak_trajectories.csv", index=False)
    print("=== region x cohort ===", flush=True)
    print(M[["region", "cohort", "n", "donors", "DP_cells", "co_OR", "int_OR", "int_p", "estimable"]]
          .to_string(index=False), flush=True)

    # ---- subclass localisation, for the regions the region-level model flags ----
    sub_rows = []
    for rg, subs in [("MEC", ["EC IT", "L4 IT"]), ("LEC", ["EC IT"]),
                     ("MTG", ["L4 IT", "L2/3 IT"]), ("V1C", ["L4 IT"])]:
        d = frames[rg]
        for sub in subs:
            q = d[d["Subclass"].astype(str) == sub]
            if len(q) < SUBCLASS_MIN_N:
                continue
            r = dict(region=rg, subclass=sub, n=len(q))
            r.update({k: v for k, v in fit(q).items() if k in ("DP_cells", "int_OR", "int_p", "estimable")})
            r["RELN_given_RORB"] = round(100 * q[q.bp == 1].rp.mean(), 3) if q.bp.sum() else np.nan
            tot = int(((d.rp == 1) & (d.bp == 1)).sum())
            r["share_of_region_DP_pct"] = round(100 * r["DP_cells"] / max(tot, 1), 1)
            sub_rows.append(r)
    S = pd.DataFrame(sub_rows)
    S.to_csv("subclass_interaction.csv", index=False)
    print("\n=== subclass localisation ===", flush=True)
    print(S.to_string(index=False), flush=True)

    # ---- same-donor cross-region correlation (MEC / MTG / PFC only) ----
    def donor_rate(rg):
        d = frames[rg]
        g = (d[d.bp == 1].groupby("Donor ID")
             .agg(frac=("rp", "mean"), nc=("rp", "size"), bk=("bk", "first")))
        return g[g.nc >= MIN_RORB]

    tabs = {rg: donor_rate(rg) for rg in ("MEC", "MTG", "PFC")}
    J = tabs["MEC"][["frac", "bk"]].rename(columns={"frac": "MEC"})
    for rg in ("MTG", "PFC"):
        J = J.join(tabs[rg][["frac"]].rename(columns={"frac": rg}), how="inner")
    rows = []
    for a, b in [("MEC", "MTG"), ("MEC", "PFC"), ("MTG", "PFC")]:
        q = J[[a, b, "bk"]].dropna()
        rho, p = spearmanr(q[a], q[b])
        ra = sm.OLS(q[a].rank(), sm.add_constant(q["bk"])).fit().resid
        rb = sm.OLS(q[b].rank(), sm.add_constant(q["bk"])).fit().resid
        prho, pp = spearmanr(ra, rb)
        n = len(q)
        hi = np.tanh(np.arctanh(rho) + norm.ppf(0.975) / np.sqrt(n - 3))
        rows.append(dict(pair="%s~%s" % (a, b), n=n, rho=round(rho, 3), p="%.4f" % p,
                         partial_rho_braak=round(prho, 3), partial_p="%.4f" % pp,
                         ci_upper=round(hi, 3),
                         detectable_at_80pct_power=round(
                             np.tanh((norm.ppf(0.975) + norm.ppf(0.80)) / np.sqrt(n - 3)), 2)))
    C = pd.DataFrame(rows)
    C.to_csv("cross_region_donor_correlation.csv", index=False)
    print("\n=== same-donor cross-region correlation ===", flush=True)
    print(C.to_string(index=False), flush=True)
    print("\nA shared donor-level process would produce a POSITIVE correlation here.", flush=True)
    print("It does not; MEC and MTG are separate region-specific associations.", flush=True)


if __name__ == "__main__":
    main()
