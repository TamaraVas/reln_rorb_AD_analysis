"""
Reciprocal test in an independent onset cohort (Rodriguez-Rodriguez 2025)
========================================================================

An independent onset-focused entorhinal study (Braak 0 vs Braak II) reported no
significant RELN/RORB co-localization by duplex in-situ hybridization. This applies
OUR same-cell test to THEIR snRNA-seq data to check whether the disagreement is
substantive or a stage-window effect.

DESIGN NOTE -- the published excitatory object is integrated with the Leng cohort,
which is our own primary discovery cohort. A naive analysis of the whole object
would therefore be partly circular. We subset to THEIR nuclei only, using the
`set` metadata column, and additionally honour the authors' own exclusion of the
two samples whose amyloid did not match their pathology group.

Their cohort spans Braak 0-II only. Our own staged data show the convergence
emerging from Braak III onward, so a null here is the predicted result, not a
contradiction. Reported in section 4.2 and Supplementary Fig. S9.

STEP 1 (R, requires Seurat) extracts counts+metadata to parquet; STEP 2 (this
script) runs the statistics. The R snippet is in the EXTRACT_R string below.

Inputs : data/rod_reln_rorb.parquet  (see EXTRACT_R; from GEO GSE287652)
Outputs: rodriguez_reciprocal.csv (printed to stdout)
Method : statsmodels Logit, depth-controlled, donor-clustered

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown.
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import DATA, require
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact

EXTRACT_R = r"""
# STEP 1 -- run once in R. Requires Seurat + arrow.
# Source: GEO GSE287652, file GSE287652_Exc_neu_integrated_with_Leng.rds.gz
library(Seurat); library(arrow)
obj  <- readRDS("data/Exc_neu_integrated_with_Leng.rds")
md   <- obj@meta.data
keep <- rownames(md)[md$set == "Rodriguez"]      # exclude integrated Leng nuclei
rna  <- GetAssayData(obj, assay = "RNA", layer = "counts")
df <- data.frame(
  reln_count = as.numeric(rna["RELN", keep]),
  rorb_count = as.numeric(rna["RORB", keep]),
  umis       = md[keep, "nCount_RNA"],
  sample     = md[keep, "sample"],
  braak      = md[keep, "BraakStage"])
write_parquet(df, "data/rod_reln_rorb.parquet")
"""

# The authors excluded these two samples from their own DE analysis because
# amyloid levels did not match their pathology group; the object encodes that
# exclusion as BraakStage == "other".
CONTROLS = ["C1", "C2", "C3", "C4"]
EARLY_AD = ["AD2b", "AD3b", "AD4", "AD5"]


def load(path=None):
    if path is None:
        path = require("rod_reln_rorb.parquet")
    d = pd.read_parquet(path)
    d["group"] = np.where(d["sample"].isin(CONTROLS), "control",
                  np.where(d["sample"].isin(EARLY_AD), "earlyAD", "excluded"))
    d = d[d.group != "excluded"].copy()
    d["reln_pos"] = (d.reln_count > 0).astype(int)
    d["rorb_pos"] = (d.rorb_count > 0).astype(int)
    d["logdepth"] = np.log10(d.umis.clip(lower=1))
    d["isAD"] = (d.group == "earlyAD").astype(int)
    return d


def same_cell_or(df, cluster_col="sample"):
    m = smf.logit("rorb_pos ~ reln_pos + logdepth", data=df).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": df[cluster_col]})
    ci = m.conf_int().loc["reln_pos"]
    return (float(np.exp(m.params["reln_pos"])), float(np.exp(ci[0])),
            float(np.exp(ci[1])), float(m.pvalues["reln_pos"]))


def per_donor(df):
    rows = []
    for dn, s in df.groupby("sample"):
        a = int(((s.reln_pos == 1) & (s.rorb_pos == 1)).sum())
        b = int(((s.reln_pos == 1) & (s.rorb_pos == 0)).sum())
        c = int(((s.reln_pos == 0) & (s.rorb_pos == 1)).sum())
        e = int(((s.reln_pos == 0) & (s.rorb_pos == 0)).sum())
        orr = fisher_exact([[a, b], [c, e]])[0] if min(b, c) > 0 else np.nan
        rows.append({"donor": dn, "group": s.group.iloc[0], "n": len(s),
                     "RELN+%": round(s.reln_pos.mean() * 100, 2),
                     "RORB+%": round(s.rorb_pos.mean() * 100, 2),
                     "DP%": round(((s.reln_pos & s.rorb_pos).mean()) * 100, 2),
                     "FisherOR": round(orr, 3) if orr == orr else np.nan})
    return pd.DataFrame(rows)


def main():
    d = load()
    print(f"Rodriguez nuclei analysed: {len(d)}  |  donors: {d['sample'].nunique()}")
    print("(integrated Leng nuclei excluded; C5/AD1 excluded per the authors' own criterion)\n")

    rows = []
    for g in ("control", "earlyAD"):
        o, lo, hi, p = same_cell_or(d[d.group == g])
        rows.append({"stratum": g, "OR": round(o, 3),
                     "CI_low": round(lo, 3), "CI_high": round(hi, 3), "p": p})
        print(f"  {g:9s} same-cell OR = {o:.3f}  95% CI [{lo:.3f}, {hi:.3f}]  p = {p:.3f}")

    m = smf.logit("rorb_pos ~ reln_pos * isAD + logdepth", data=d).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": d["sample"]})
    io = float(np.exp(m.params["reln_pos:isAD"]))
    ip = float(m.pvalues["reln_pos:isAD"])
    ci = m.conf_int().loc["reln_pos:isAD"]
    print(f"\n  interaction OR = {io:.3f}  95% CI "
          f"[{np.exp(ci[0]):.3f}, {np.exp(ci[1]):.3f}]  p = {ip:.3f}")
    rows.append({"stratum": "interaction (AD vs control)", "OR": round(io, 3),
                 "CI_low": round(float(np.exp(ci[0])), 3),
                 "CI_high": round(float(np.exp(ci[1])), 3), "p": ip})

    print("\nPer-donor breakdown (watch for single-donor drivers):")
    print(per_donor(d).to_string(index=False))

    print("\nInterpretation: this cohort spans Braak 0-II. Our SEA-AD data give a same-cell\n"
          "OR well below 1 at those stages, crossing above 1 only at Braak V-VI, so a null\n"
          "here is what our own model predicts -- the cohorts agree where they overlap.\n"
          "Power is limited (4 donors per group after exclusions); read as consistent-with,\n"
          "not as proof of, an absence.")
    pd.DataFrame(rows).to_csv("rodriguez_reciprocal.csv", index=False)
    print("\nwrote rodriguez_reciprocal.csv")


if __name__ == "__main__":
    main()
