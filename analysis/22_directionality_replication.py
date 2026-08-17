"""
Directionality asymmetry: why no directional claim is made
==========================================================

In the Leng discovery cohort the two conditional acquisition slopes are
asymmetric: RELN co-expression rises steeply within RORB+ neurons while RORB
changes little within RELN+ neurons. Read naively that suggests the double-
positive state is reached by adding RELN to a stable RORB+ entorhinal identity.

The manuscript does NOT advance that reading, and this script is the reason.
Two independent problems:

  1. IDENTIFIABILITY. Acquisition and differential survival of RELN+RORB- neurons
     predict the SAME cross-sectional asymmetry. Autopsy tissue sampled once per
     donor cannot distinguish them. No test in this script can fix that; it is
     stated here so the limitation travels with the code.

  2. NON-REPLICATION. The DIRECTION of the asymmetry reverses in the larger
     staged cohort. In SEA-AD MEC, RORB acquisition within RELN+ neurons is
     significant while RELN acquisition within RORB+ neurons is not -- the
     opposite assignment to Leng. Within the stable EC-IT subclass the
     RELN-acquisition slope is flat.

Why the reversal is expected rather than contradictory: the two conditional
slopes differ only by the difference in the markers' MARGINAL trends, and those
marginals differ between cohorts (RORB+ prevalence is flat end-to-end in Leng but
rises in SEA-AD). The asymmetry is therefore a description of a cohort, not a
cohort-independent property.

What DOES replicate across both cohorts is the convergence itself -- that
co-occupancy increases with tau stage. This script reports both quantities side
by side so the distinction is auditable.

NOTE ON POSITIVITY (see DATA.md): Leng marker columns are RAW COUNTS (`> 0`).
The SEA-AD extract columns are DEPTH-NORMALISED truncated values for which the
rule is `>= 1`; the script asserts the published MEC prevalences. The two rules
are not interchangeable.

Reported in manuscript section 3.12 and Discussion section 4.1.

Inputs : data/leng2021_ec_full.h5ad
         data/seaad_mec_extract.parquet
Outputs: directionality_replication.csv (printed to stdout)
Method : statsmodels Logit with donor-clustered robust covariance

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown.
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import DATA, require
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.formula.api as smf

RELN_ENS = "ENSG00000189056"
RORB_ENS = "ENSG00000198963"
BRAAK_MAP = {"Braak 0": 0, "Braak I": 1, "Braak II": 2,
             "Braak III": 3, "Braak IV": 4, "Braak V": 5, "Braak VI": 6}


def fit(df, formula, term):
    """Logistic fit with donor-clustered robust errors. Returns OR, CI, p, n."""
    groups = df["donor"].astype("category").cat.codes.values
    m = smf.logit(formula, data=df).fit(disp=0, cov_type="cluster",
                                        cov_kwds={"groups": groups})
    ci = m.conf_int()
    return dict(OR=round(float(np.exp(m.params[term])), 3),
                ci=f"{np.exp(ci.loc[term, 0]):.2f}-{np.exp(ci.loc[term, 1]):.2f}",
                p=float(m.pvalues[term]), n=len(df),
                donors=int(df["donor"].nunique()))


def load_leng():
    """Leng EC excitatory nuclei with RAW-COUNT positivity."""
    import anndata as ad
    A = ad.read_h5ad(require("leng2021_ec_full.h5ad"))
    exc = A[A.obs["clusterCellType"].astype(str) == "Exc"].copy()

    def counts(ens):
        x = exc[:, ens].X
        return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

    X = exc.layers["counts"] if "counts" in exc.layers else exc.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)

    donor_col = next(c for c in exc.obs.columns
                     if c.lower() in ("sampleid", "donor", "patient", "donor_id"))
    braak_col = next(c for c in exc.obs.columns if "braak" in c.lower())

    df = pd.DataFrame({
        "rp": (counts(RELN_ENS) > 0).astype(int),   # RAW COUNT rule, see DATA.md
        "bp": (counts(RORB_ENS) > 0).astype(int),
        "donor": exc.obs[donor_col].astype(str).values,
        "bkn": pd.to_numeric(exc.obs[braak_col].astype(str)
                             .str.extract(r"(\d+)")[0], errors="coerce").values,
        "logd": np.log10(np.clip(np.asarray(X.sum(1)).ravel(), 1, None)),
    })
    return df[df["bkn"].notna()].copy()


def load_seaad():
    """SEA-AD MEC excitatory nuclei with the DEPTH-NORMALISED positivity rule."""
    df = pd.read_parquet(require("seaad_mec_extract.parquet"))
    df["rp"] = (df["reln"] >= 1).astype(int)
    df["bp"] = (df["rorb"] >= 1).astype(int)
    df["bkn"] = df["Braak"].astype(str).map(BRAAK_MAP)
    df["logd"] = np.log10(df["umis"].clip(lower=1))
    exc = df[df["Class"].astype(str).str.contains("Glutamatergic", na=False)]
    exc = exc[exc["bkn"].notna()].copy()

    reln_pct, rorb_pct = 100 * exc["rp"].mean(), 100 * exc["bp"].mean()
    assert abs(reln_pct - 13.3) < 0.2 and abs(rorb_pct - 20.2) < 0.2, (
        f"positivity rule mis-calibrated: RELN+ {reln_pct:.1f}% (expect 13.3), "
        f"RORB+ {rorb_pct:.1f}% (expect 20.2). See DATA.md before changing this."
    )
    return exc


def directionality(df, cohort, subclass_col=None, subclass=None):
    """The two conditional acquisition slopes, plus the convergence interaction."""
    if subclass is not None:
        df = df[df[subclass_col].astype(str) == subclass]
    rows = []

    # RELN acquisition among RORB+ cells
    q = df[df["bp"] == 1]
    rows.append(dict(cohort=cohort, quantity="RELN acquisition within RORB+",
                     **fit(q, "rp ~ bkn + logd", "bkn")))
    # RORB acquisition among RELN+ cells
    q = df[df["rp"] == 1]
    rows.append(dict(cohort=cohort, quantity="RORB acquisition within RELN+",
                     **fit(q, "bp ~ bkn + logd", "bkn")))
    # the convergence itself -- the quantity that DOES replicate
    rows.append(dict(cohort=cohort, quantity="convergence (RELNxBraak interaction)",
                     **fit(df, "bp ~ rp * bkn + logd", "rp:bkn")))

    # marginal trends, which is what makes the two conditionals differ
    for marker, term in [("RELN+ marginal", "rp"), ("RORB+ marginal", "bp")]:
        rows.append(dict(cohort=cohort, quantity=marker,
                         **fit(df, f"{term} ~ bkn + logd", "bkn")))
    return rows


def main():
    rows = []
    rows += directionality(load_leng(), "Leng EC (discovery)")

    seaad = load_seaad()
    rows += directionality(seaad, "SEA-AD MEC (staged)")
    rows += directionality(seaad, "SEA-AD MEC, EC-IT subclass",
                           subclass_col="Subclass", subclass="EC IT")

    out = pd.DataFrame(rows)
    print("=== Conditional acquisition slopes and the convergence, by cohort ===")
    print(out.to_string(index=False))
    print("\nRead the two 'acquisition' rows per cohort together: their ORDERING")
    print("reverses between cohorts, while the 'convergence' row is positive and")
    print("significant in both. Only the latter supports a claim.")
    out.to_csv("directionality_replication.csv", index=False)


if __name__ == "__main__":
    main()
