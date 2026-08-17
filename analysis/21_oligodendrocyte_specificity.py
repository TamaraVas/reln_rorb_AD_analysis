"""
Oligodendrocyte RELN: specificity control for the neuronal convergence
======================================================================

RELN+ oligodendrocytes rise with Braak stage. This script establishes what that
observation is and is NOT, so it can be reported as a specificity control for the
neuronal RELN-RORB convergence rather than as evidence about it.

Four tests, in the order that bounds the interpretation:

  1. CLASS SPECIFICITY -- the same depth-controlled model in every non-neuronal
     subclass. Only oligodendrocytes exceed OR = 1; astrocytes, OPCs, immune,
     endothelial and VLMC/perivascular cells are flat or declining. This is the
     argument against ambient RNA: a contamination floor rising with disease
     would lift every class, not one.

  2. REGIONAL SPECIFICITY -- the same model in all four SEA-AD regions. The
     oligodendrocyte rise is BRAIN-WIDE (present in primary visual cortex, where
     the neuronal convergence cannot even be defined for lack of RELN), whereas
     the neuronal convergence is entorhinal. This double dissociation is the
     reason the observation belongs in the specificity section.

  3. MARKER SPECIFICITY -- RORB in oligodendrocytes, and the within-oligodendrocyte
     RELN x Braak interaction. Both null: these cells acquire RELN detection
     without acquiring the paired state that defines the neuronal phenotype.

  4. PATHOLOGY AXIS -- tau vs amyloid vs dementia vs APOE, singly and jointly.
     Unlike the neuronal convergence the effect is NOT tau-specific: amyloid
     tracks equally and neither survives adjustment for the other (they are
     collinear at donor level). Dementia status is the largest single effect.

  5. DEPTH -- the headline effect within library-size quintiles. It survives in
     every quintile, and sequencing depth DECLINES with Braak in these cells, so
     depth works against the effect rather than producing it.

NOTE ON POSITIVITY (see DATA.md): the SEA-AD extract marker columns are
DEPTH-NORMALISED truncated values that resemble small counts. Positivity is
`>= 1` on those columns, which reproduces the manuscript prevalences; the script
asserts this on the excitatory pool before proceeding. Leng raw counts use `> 0`.
Do not "fix" one rule to match the other.

SCOPE NOTE: this script deliberately stops at the specificity controls. Whether a
distinct RELN+ oligodendrocyte transcriptional STATE exists is a separate
question, and in these data it is confounded by sequencing depth (RELN+
oligodendrocytes carry ~2x the library size of RELN- ones, and a differential-
expression signature retains only ~14% of its effect after depth matching). That
analysis is not part of this manuscript.

Reported in manuscript section 3.6.

Inputs : data/seaad_mec_extract.parquet
         data/seaad_lec_extract.parquet
         data/seaad_hip_extract.parquet
         data/seaad_v1c_extract.parquet
Outputs: oligo_class_specificity.csv
         oligo_regional_specificity.csv
         oligo_pathology_axes.csv   (all printed to stdout)
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
import statsmodels.formula.api as smf

BRAAK_MAP = {"Braak 0": 0, "Braak I": 1, "Braak II": 2,
             "Braak III": 3, "Braak IV": 4, "Braak V": 5, "Braak VI": 6}
THAL_MAP = {"Thal 0": 0, "Thal 1": 1, "Thal 2": 2, "Thal 3": 3, "Thal 4": 4, "Thal 5": 5}
CERAD_MAP = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}

NON_NEURONAL = ["Oligodendrocyte", "Astrocyte", "OPC", "Immune",
                "Endothelial", "VLMC & Perivascular"]
REGIONS = ["MEC", "LEC", "HIP", "V1C"]
MIN_CELLS = 500          # below this a subclass/region estimate is not reported


def load_region(region):
    """Load one SEA-AD regional extract and derive the modelling columns."""
    path = require(f"seaad_{region.lower()}_extract.parquet")
    df = pd.read_parquet(path)
    df["rp"] = (df["reln"] >= 1).astype(int)     # depth-normalised rule, see DATA.md
    df["bp"] = (df["rorb"] >= 1).astype(int)
    df["bkn"] = df["Braak"].astype(str).map(BRAAK_MAP)
    df["logd"] = np.log10(df["umis"].clip(lower=1))
    return df[df["bkn"].notna()].copy()


def assert_positivity(df, region):
    """The MEC excitatory prevalences are the published calibration point."""
    if region != "MEC":
        return
    exc = df[df["Class"].astype(str).str.contains("Glutamatergic", na=False)]
    reln_pct, rorb_pct = 100 * exc["rp"].mean(), 100 * exc["bp"].mean()
    assert abs(reln_pct - 13.3) < 0.2 and abs(rorb_pct - 20.2) < 0.2, (
        f"positivity rule mis-calibrated: RELN+ {reln_pct:.1f}% (expect 13.3), "
        f"RORB+ {rorb_pct:.1f}% (expect 20.2). See DATA.md before changing this."
    )


def fit(df, formula, term):
    """Logistic fit with donor-clustered robust errors. Returns OR, CI, p."""
    groups = df["donor"].astype("category").cat.codes.values
    m = smf.logit(formula, data=df).fit(disp=0, cov_type="cluster",
                                        cov_kwds={"groups": groups})
    ci = m.conf_int()
    return (float(np.exp(m.params[term])),
            float(np.exp(ci.loc[term, 0])), float(np.exp(ci.loc[term, 1])),
            float(m.pvalues[term]))


def class_specificity(mec):
    """Test 1 -- is the rise oligodendrocyte-specific among non-neuronal cells?"""
    rows = []
    for sc in NON_NEURONAL:
        q = mec[mec["Subclass"].astype(str) == sc]
        if len(q) < MIN_CELLS:
            rows.append(dict(subclass=sc, n=len(q), note="below MIN_CELLS"))
            continue
        orr, lo, hi, p = fit(q, "rp ~ bkn + logd", "bkn")
        rows.append(dict(subclass=sc, n=len(q), donors=q["donor"].nunique(),
                         RELN_pct=round(100 * q["rp"].mean(), 2),
                         OR=round(orr, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                         p=p))
    return pd.DataFrame(rows).sort_values("OR", ascending=False)


def regional_specificity():
    """Test 2 -- oligodendrocyte rise vs neuronal convergence, all four regions."""
    rows = []
    for region in REGIONS:
        df = load_region(region)
        assert_positivity(df, region)

        oli = df[df["Subclass"].astype(str) == "Oligodendrocyte"]
        if len(oli) >= MIN_CELLS:
            orr, lo, hi, p = fit(oli, "rp ~ bkn + logd", "bkn")
            oli_stats = dict(oli_n=len(oli), oli_OR=round(orr, 3),
                             oli_ci=f"{lo:.2f}-{hi:.2f}", oli_p=p)
        else:
            oli_stats = dict(oli_n=len(oli), oli_OR=np.nan, oli_ci="", oli_p=np.nan)

        # the contrasting neuronal quantity, same region, same model family
        exc = df[df["Class"].astype(str).str.contains("Glutamatergic", na=False)]
        n_dp = int(((exc["rp"] == 1) & (exc["bp"] == 1)).sum())
        if n_dp >= 50 and exc["rp"].mean() > 0.005 and exc["bp"].mean() > 0.005:
            orr, lo, hi, p = fit(exc, "bp ~ rp * bkn + logd", "rp:bkn")
            neu_stats = dict(exc_n=len(exc), conv_OR=round(orr, 3),
                             conv_ci=f"{lo:.2f}-{hi:.2f}", conv_p=p)
        else:
            # no RELN substrate -> the convergence is undefined, not null
            neu_stats = dict(exc_n=len(exc), conv_OR=np.nan,
                             conv_ci="undefined (no substrate)", conv_p=np.nan)

        rows.append(dict(region=region, **oli_stats, **neu_stats))
    return pd.DataFrame(rows)


def marker_and_pathology(mec):
    """Tests 3-5 -- marker specificity, pathology axis, and depth robustness."""
    oli = mec[mec["Subclass"].astype(str) == "Oligodendrocyte"].copy()
    oli["thal"] = oli["Thal"].astype(str).map(THAL_MAP)
    oli["cerad"] = oli["CERAD"].astype(str).map(CERAD_MAP)
    # NOTE: exact match, not a substring test -- the negative level is literally
    # "No dementia", which CONTAINS "dementia" and would score every cell positive.
    dem_str = oli["Cognitive"].astype(str).str.strip().str.lower()
    assert set(dem_str.unique()) <= {"dementia", "no dementia"}, (
        f"unexpected Cognitive levels: {sorted(dem_str.unique())}")
    oli["dem"] = (dem_str == "dementia").astype(int)
    oli["e4"] = oli["APOE"].astype(str).str.count("4")

    rows = []
    # Test 3: is it RELN-specific, or a convergence like the neuronal one?
    orr, lo, hi, p = fit(oli, "bp ~ bkn + logd", "bkn")
    rows.append(dict(test="RORB+ vs Braak (marker specificity)", OR=round(orr, 3),
                     ci=f"{lo:.2f}-{hi:.2f}", p=p))
    orr, lo, hi, p = fit(oli, "bp ~ rp * bkn + logd", "rp:bkn")
    rows.append(dict(test="RELNxBraak interaction within oligodendrocytes",
                     OR=round(orr, 3), ci=f"{lo:.2f}-{hi:.2f}", p=p))

    # Test 4: which pathology axis, singly then jointly
    for var, label in [("bkn", "Braak tau stage"), ("thal", "Thal amyloid phase"),
                       ("cerad", "CERAD neuritic plaque"), ("dem", "Dementia (vs none)"),
                       ("e4", "APOE-e4 allele dose")]:
        q = oli[oli[var].notna()]
        orr, lo, hi, p = fit(q, f"rp ~ {var} + logd", var)
        rows.append(dict(test=label, OR=round(orr, 3), ci=f"{lo:.2f}-{hi:.2f}", p=p))

    q = oli[oli["bkn"].notna() & oli["thal"].notna()]
    for var, label in [("bkn", "Braak, amyloid-adjusted"), ("thal", "Amyloid, Braak-adjusted")]:
        orr, lo, hi, p = fit(q, "rp ~ bkn + thal + logd", var)
        rows.append(dict(test=label, OR=round(orr, 3), ci=f"{lo:.2f}-{hi:.2f}", p=p))

    # Test 5: depth. Does the effect hold within library-size quintiles, and
    # which way does depth itself move with stage?
    bins = np.quantile(oli["umis"], [0, .2, .4, .6, .8, 1.0])
    oli["depth_bin"] = pd.cut(oli["umis"], bins, include_lowest=True)
    for b, q in oli.groupby("depth_bin", observed=True):
        if q["bkn"].nunique() < 4:
            continue
        orr, lo, hi, p = fit(q, "rp ~ bkn", "bkn")
        rows.append(dict(test=f"Braak within depth quintile {int(b.left)}-{int(b.right)}",
                         OR=round(orr, 3), ci=f"{lo:.2f}-{hi:.2f}", p=p))
    rho = oli[["bkn", "umis"]].corr(method="spearman").loc["bkn", "umis"]
    rows.append(dict(test="depth vs Braak (Spearman rho, NOT an OR)",
                     OR=round(float(rho), 3), ci="", p=np.nan))
    return pd.DataFrame(rows)


def main():
    mec = load_region("MEC")
    assert_positivity(mec, "MEC")

    cls = class_specificity(mec)
    print("=== Test 1: RELN+ vs Braak by non-neuronal subclass (SEA-AD MEC) ===")
    print(cls.to_string(index=False), "\n")

    reg = regional_specificity()
    print("=== Test 2: oligodendrocyte rise (brain-wide) vs neuronal convergence "
          "(entorhinal) ===")
    print(reg.to_string(index=False), "\n")

    axes = marker_and_pathology(mec)
    print("=== Tests 3-5: marker specificity, pathology axis, depth robustness ===")
    print(axes.to_string(index=False))

    cls.to_csv("oligo_class_specificity.csv", index=False)
    reg.to_csv("oligo_regional_specificity.csv", index=False)
    axes.to_csv("oligo_pathology_axes.csv", index=False)


if __name__ == "__main__":
    main()
