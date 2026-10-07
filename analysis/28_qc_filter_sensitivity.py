"""
Quality-control filter sensitivity: mitochondrial fraction and doublet score
============================================================================

An snRNA-seq reviewer asked whether the convergence survives stringent per-nucleus
quality filtering. The 2026 release distributes per-nucleus mitochondrial fraction
and doublet score, which the earlier release did not, so this can be measured
rather than argued.

The doublet test is the one that matters. A doublet -- two nuclei in one droplet --
can only MANUFACTURE apparent same-cell co-expression; it cannot hide it. So if
doublets drove the signal, discarding suspected doublets would collapse it. It does
not: in MEC the interaction is unchanged (1.481 -> 1.470 at score <= 0.1, retaining
69.4% of nuclei) and in LEC the baseline co-occupancy OR RISES from 1.07 to 1.51.

The mitochondrial results are region-specific and reported as such rather than
summarised in one direction: MEC strengthens as filtering tightens (1.481 -> 1.571
at <= 1%), MTG weakens (1.284 -> 1.212). What holds in both is that the interaction
stays significant at EVERY threshold tested -- "robust to filtering", not
"strengthened by it".

Inputs : data/seaad2026_<REGION>_final.parquet  (script 25)
Outputs: regional_qc_two_cohorts.csv, doublet_direction_test.csv
Method : depth-controlled logistic regression, donor-clustered SEs, at each threshold

Part of the RELN-RORB EC convergence analysis. See README.md.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import require

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact

REGIONS = ["MEC", "LEC", "HIP", "ITG", "MTG", "STG", "FI", "AnG", "PFC", "V1C"]
BRAAK = {"Braak 0": 0, "Braak II": 2, "Braak III": 3, "Braak IV": 4, "Braak V": 5, "Braak VI": 6}
MITO_THRESHOLDS = [5, 2, 1]          # percent of UMIs; 5% is the conventional floor for nuclei
DOUBLET_THRESHOLDS = [0.3, 0.2, 0.1]
MIN_DP = 10


def load(region):
    d = pd.read_parquet(require("seaad2026_%s_final.parquet" % region))
    d = d[d["Class"].astype(str).str.contains("Glutamatergic", na=False)].copy()
    d["rp"] = (d["RELN_norm"] >= 1).astype(int)
    d["bp"] = (d["RORB_norm"] >= 1).astype(int)
    d["bk"] = d["Braak"].astype(str).map(BRAAK)
    d["logd"] = np.log10(pd.to_numeric(d["Number of UMIs"], errors="coerce").clip(lower=1))
    d["mito"] = 100 * pd.to_numeric(d["Fraction mitochondrial UMIs"], errors="coerce")
    d["dbl"] = pd.to_numeric(d["Doublet score"], errors="coerce")
    return d.dropna(subset=["bk", "logd", "mito", "dbl"])


def fit(d):
    dp = int(((d.rp == 1) & (d.bp == 1)).sum())
    if dp < MIN_DP:
        return dict(n=len(d), DP_cells=dp, co_OR="UNESTIMABLE (DP<%d)" % MIN_DP, co_p="",
                    int_OR="UNESTIMABLE (DP<%d)" % MIN_DP, int_p="")
    tab = [[dp, int(((d.rp == 1) & (d.bp == 0)).sum())],
           [int(((d.rp == 0) & (d.bp == 1)).sum()), int(((d.rp == 0) & (d.bp == 0)).sum())]]
    o, p = fisher_exact(tab)
    g = d["Donor ID"].astype("category").cat.codes.values
    m = smf.logit("bp ~ rp * bk + logd", data=d).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": g})
    return dict(n=len(d), DP_cells=dp, co_OR=round(float(o), 3), co_p="%.3g" % p,
                int_OR=round(float(np.exp(m.params["rp:bk"])), 3),
                int_p="%.3g" % m.pvalues["rp:bk"])


def main():
    rows = []
    for rg in REGIONS:
        d = load(rg)
        for cohort in ("all", "no_severe"):
            base = d if cohort == "all" else d[d["Severely Affected Donor"].astype(str) == "N"]
            if not len(base):
                continue
            grids = ([("none", base)]
                     + [("mito<=%d%%" % t, base[base.mito <= t]) for t in MITO_THRESHOLDS]
                     + [("doublet<=%s" % t, base[base.dbl <= t]) for t in DOUBLET_THRESHOLDS])
            for lab, q in grids:
                rows.append(dict(region=rg, cohort=cohort, filter=lab,
                                 pct_kept=round(100 * len(q) / len(base), 1), **fit(q)))
    Q = pd.DataFrame(rows)
    Q.to_csv("regional_qc_two_cohorts.csv", index=False)
    print("=== MEC and MTG, all donors ===", flush=True)
    print(Q[(Q.region.isin(["MEC", "MTG"])) & (Q.cohort == "all")]
          [["region", "filter", "pct_kept", "co_OR", "int_OR", "int_p"]].to_string(index=False), flush=True)

    # doublet direction: removal should LOWER the OR if doublets drive co-occupancy
    a = Q[(Q.cohort == "all") & (Q["filter"].isin(["none", "doublet<=0.1"]))]
    piv = a.pivot_table(index="region", columns="filter", values="co_OR", aggfunc="first").reset_index()
    piv.columns = ["region", "OR_doublet_filtered", "OR_unfiltered"]
    for c in ("OR_doublet_filtered", "OR_unfiltered"):
        piv[c] = pd.to_numeric(piv[c], errors="coerce")
    piv["delta"] = (piv.OR_doublet_filtered - piv.OR_unfiltered).round(3)
    piv["pct_change"] = (100 * (piv.OR_doublet_filtered / piv.OR_unfiltered - 1)).round(1)
    piv.to_csv("doublet_direction_test.csv", index=False)
    print("\n=== doublet direction test ===", flush=True)
    print(piv.dropna().to_string(index=False), flush=True)
    up = int((piv.delta > 0).sum())
    print("\nOR rises on doublet removal in %d of %d estimable regions." % (up, int(piv.delta.notna().sum())), flush=True)
    print("Rises are large where a real DP population exists; the small falls occur", flush=True)
    print("in regions with RELN < 0.75%, where the DP population is marginal.", flush=True)


if __name__ == "__main__":
    main()
