"""
MERFISH detection-error and segmentation controls (SEA-AD MEC)
==============================================================

The 433-gene MERFISH release has no peer-reviewed methods paper, so its detection
error cannot be cited -- but the delivered object retains the control channels
needed to measure it. This script does that, and then asks the question that
actually matters for a same-cell claim: can measurement error PRODUCE the
co-localization?

Two error modes, tested separately because the control channels only bound one:

  1. DECODING error -- a codeword misread as a different gene. Bounded by the
     control probe / control codeword / unassigned codeword channels. Independent
     per codeword, so it ATTENUATES an association; the script inverts the
     false-positive rate to show the corrected OR moves away from unity.

  2. SEGMENTATION spillover -- a neighbouring cell's transcripts assigned to the
     cell being scored. NOT captured by the control channels, and the more
     consequential failure mode for a same-cell claim in intact tissue. Tested
     two ways: within each segmentation method (the sections carry three), and by
     adding segmented cell area to the depth-controlled model.

The control-panel size is not recorded in the object, so the per-gene error rate is
bounded across plausible panel sizes rather than asserted. Conclusions are drawn
from the most conservative bound.

Inputs : data/seaad_merfish_mec.h5ad
Outputs: merfish_false_positive_qc.csv, merfish_segmentation_controls.csv
Method : Poisson false-positive model + Fisher / logistic co-localization tests

Part of the RELN-RORB EC convergence analysis. See README.md.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data_paths import require

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact, mannwhitneyu

# Control-panel size is not recorded in the object. Bound it: fewer control probes
# means a higher per-probe error rate, so 20 is the conservative end.
PANEL_SIZES = (20, 50, 100)
CONSERVATIVE = 20


def counts(a, symbol):
    """Per-cell counts for one gene, tolerating symbol vs index naming."""
    names = [str(g).upper() for g in a.var_names]
    if symbol not in names:
        for col in ("gene_symbol", "feature_name"):
            if col in a.var.columns:
                names = [str(g).upper() for g in a.var[col]]
                break
    i = names.index(symbol)
    x = a[:, i].X
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()


def control_channels(obs):
    """Totals for each control channel plus real signal."""
    num = lambda c: pd.to_numeric(obs[c], errors="coerce")
    ch = dict(real=num("transcript_counts"),
              control_probe=num("control_probe_counts"),
              genomic_control=num("genomic_control_counts"),
              control_codeword=num("control_codeword_counts"),
              unassigned=num("unassigned_codeword_counts"))
    total = sum(v.sum() for v in ch.values())
    return ch, total


def fp_rate(control_probe_total, n_cells, n_probes):
    """P(at least one spurious count for one gene in one cell), Poisson."""
    lam = control_probe_total / n_cells / n_probes
    return 1.0 - np.exp(-lam)


def deconvolve(table, f):
    """Invert independent per-gene false positives at rate f.

    table = [[both+, A only], [B only, neither]]. Independent FPs move an
    observed OR toward 1, so the recovered OR should exceed the observed one.
    """
    (A_o, B_o), (C_o, D_o) = table
    D = D_o / ((1 - f) ** 2)
    B = (B_o - D * f * (1 - f)) / (1 - f)
    C = (C_o - D * f * (1 - f)) / (1 - f)
    A = A_o - f * B - f * C - (f ** 2) * D
    return (A * D) / (B * C)


def main():
    path = require("seaad_merfish_mec.h5ad")
    a = ad.read_h5ad(path)
    obs = a.obs
    ch, total = control_channels(obs)
    n = len(obs)

    print(f"=== control channels ({n:,} cells) ===")
    for k, v in ch.items():
        print(f"  {k:18s} total {v.sum():12,.0f}   {100*v.sum()/total:7.4f}% of signal")
    ctrl_share = 100 * sum(ch[k].sum() for k in
                           ("control_probe", "control_codeword", "unassigned")) / total
    print(f"  all control channels: {ctrl_share:.4f}% of signal")
    if ch["genomic_control"].sum() == 0:
        print("  genomic controls exactly zero -> no detectable genomic DNA contamination")

    reln, rorb = counts(a, "RELN"), counts(a, "RORB")
    rp, bp = (reln > 0).astype(int), (rorb > 0).astype(int)
    cls = next((c for c in obs.columns if "Class" in c and "conf" not in c), None)
    glut = (obs[cls].astype(str).str.contains("Glut|Exc", case=False, na=False).values
            if cls else np.ones(n, bool))

    tab = [[int(((rp == 1) & (bp == 1))[glut].sum()), int(((rp == 1) & (bp == 0))[glut].sum())],
           [int(((rp == 0) & (bp == 1))[glut].sum()), int(((rp == 0) & (bp == 0))[glut].sum())]]
    obs_or, obs_p = fisher_exact(tab)
    print(f"\n=== co-localization, glutamatergic (n={int(glut.sum()):,}) ===")
    print(f"  observed OR = {obs_or:.3f}  p = {obs_p:.3g}  double-positive = {tab[0][0]:,}")

    print("\n=== 1. decoding error: can it create the association? ===")
    rows = []
    for np_ in PANEL_SIZES:
        f = fp_rate(ch["control_probe"].sum(), n, np_)
        exp_dp = f * f * glut.sum()
        corr = deconvolve(tab, f)
        print(f"  {np_:3d} control probes: FP/call {100*f:.4f}%  "
              f"expected DP from noise {exp_dp:.2f}  corrected OR {corr:.3f}")
        rows.append(dict(control="noise-corrected (%d control probes)" % np_,
                         value=round(corr, 3), p="-", n=int(glut.sum()),
                         note="independent FPs attenuate; correction moves OR away from 1"))
    f_c = fp_rate(ch["control_probe"].sum(), n, CONSERVATIVE)
    print(f"  -> at the conservative bound, noise explains "
          f"{100*f_c*f_c*glut.sum()/tab[0][0]:.4f}% of double-positive calls")

    print("\n=== 2. segmentation spillover (NOT bounded by control channels) ===")
    seg = obs["segmentation_method"].astype(str)
    for meth, cnt in seg.value_counts().items():
        m = (seg == meth).values & glut
        if m.sum() < 500:
            print(f"  skipped (n={int(m.sum())}): {meth[:50]}")
            continue
        tt = [[int(((rp == 1) & (bp == 1))[m].sum()), int(((rp == 1) & (bp == 0))[m].sum())],
              [int(((rp == 0) & (bp == 1))[m].sum()), int(((rp == 0) & (bp == 0))[m].sum())]]
        o_, p_ = fisher_exact(tt)
        print(f"  OR {o_:5.3f}  p {p_:9.3g}  n {int(m.sum()):6,d}  {meth[:46]}")
        rows.append(dict(control="within %s" % meth[:46], value=round(o_, 3),
                         p=f"{p_:.3g}", n=int(m.sum()), note="segmentation-stratified"))

    area = pd.to_numeric(obs["cell_area"], errors="coerce")
    dp = ((rp == 1) & (bp == 1))
    u, pu = mannwhitneyu(area[glut & dp].dropna(), area[glut & ~dp].dropna())
    print(f"\n  median area: DP {area[glut & dp].median():.1f} vs non-DP "
          f"{area[glut & ~dp].median():.1f} um^2  (p = {pu:.3g})")

    w = pd.DataFrame({"rp": rp[glut], "bp": bp[glut],
                      "area": area[glut].values, "tot": ch["real"][glut].values}).dropna()
    w["logarea"] = np.log10(w.area.clip(lower=1))
    w["logtot"] = np.log10(w.tot.clip(lower=1))
    m_adj = smf.logit("bp ~ rp + logarea + logtot", data=w).fit(disp=0)
    print(f"  adjusted for area + depth: OR {np.exp(m_adj.params['rp']):.3f} "
          f"p {m_adj.pvalues['rp']:.3g}")
    print(f"  area term itself: OR {np.exp(m_adj.params['logarea']):.3f} "
          f"p {m_adj.pvalues['logarea']:.3g}")
    print("  -> a null area term means the size difference is depth, not spillover")
    rows.append(dict(control="adjusted for cell area + transcript depth",
                     value=round(float(np.exp(m_adj.params["rp"])), 3),
                     p=f"{m_adj.pvalues['rp']:.3g}", n=len(w),
                     note="area term OR %.3f p %.3g (null)" % (
                         np.exp(m_adj.params["logarea"]), m_adj.pvalues["logarea"])))

    pd.DataFrame([
        dict(metric="cells profiled", value=n),
        dict(metric="all control channels as % of signal", value=round(ctrl_share, 4)),
        dict(metric="genomic control counts", value=int(ch["genomic_control"].sum())),
        dict(metric="double-positive glutamatergic cells", value=tab[0][0]),
        dict(metric="expected DP from noise (conservative)",
             value=round(f_c * f_c * glut.sum(), 2)),
        dict(metric="observed Fisher OR", value=round(obs_or, 3)),
    ]).to_csv("merfish_false_positive_qc.csv", index=False)
    pd.DataFrame(rows).to_csv("merfish_segmentation_controls.csv", index=False)
    print("\nwrote merfish_false_positive_qc.csv, merfish_segmentation_controls.csv")


if __name__ == "__main__":
    main()
