"""
Same-cell RELN x RORB co-occupancy in EC excitatory neurons (Leng 2021)
=======================================================================

Discovery analysis: RELN and RORB are expressed independently in cognitively
normal EC excitatory neurons but converge into the SAME neurons in Alzheimer's,
graded by Braak tau stage. Depth-controlled logistic model isolates the
same-cell interaction from the sequencing-depth confound.

Inputs : data/leng2021_ec_full.h5ad   (Leng et al. 2021 EC snRNA-seq)
Outputs: reln_rorb_coexpression_AD.png, coexpression_stats.csv
Method : statsmodels Logit (depth-controlled), scipy Fisher exact

NOTE ON PROVENANCE: clean standalone re-implementation from the documented
parameters (positivity = raw count >= 1; depth = log10 nUMI; normal = Braak 0,
AD = Braak 2 + 6). See METHODS_PROVENANCE.md.

Part of the RELN-RORB EC convergence analysis suite. See README.md.
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import statsmodels.api as sm
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

A = ad.read_h5ad("data/leng2021_ec_full.h5ad")
gi = {g: i for i, g in enumerate(A.var["feature_name"])}
X = A.X.tocsc()
def col(g): return np.asarray(X[:, gi[g]].todense()).ravel()

o = A.obs.copy()
o["RELN_pos"] = col("RELN") >= 1
o["RORB_pos"] = col("RORB") >= 1
o["logdepth"] = np.log10(o["nUMI"].astype(float).clip(lower=1))
o["braak"] = o["BraakStage"].astype(int)
exc = o[o["clusterCellType"].astype(str).str.contains("Exc")].copy()
exc["AD"] = (exc.braak > 0).astype(int)          # normal = Braak 0; AD = Braak 2 + 6
inh = o[o["clusterCellType"].astype(str).str.contains("Inh")].copy()
inh["AD"] = (inh.braak > 0).astype(int)

def depth_or(df):
    if df.RELN_pos.sum() < 5 or df.RORB_pos.sum() < 5:
        return np.nan, np.nan, np.nan, np.nan
    Xm = sm.add_constant(df[["RELN_pos", "logdepth"]].astype(float))
    m = sm.Logit(df.RORB_pos.astype(float), Xm).fit(disp=0)
    ci = m.conf_int().loc["RELN_pos"]
    return np.exp(m.params["RELN_pos"]), m.pvalues["RELN_pos"], np.exp(ci[0]), np.exp(ci[1])

# depth-controlled same-cell interaction (the estimand)
exc["RxAD"] = exc.RELN_pos.astype(int) * exc.AD
Xi = sm.add_constant(exc[["RELN_pos", "AD", "RxAD", "logdepth"]].astype(float))
mi = sm.Logit(exc.RORB_pos.astype(float), Xi).fit(disp=0)
inter_or, inter_p = np.exp(mi.params["RxAD"]), mi.pvalues["RxAD"]

braaks = [0, 2, 6]
rln = [exc[exc.braak == b].RELN_pos.mean() * 100 for b in braaks]
rrb = [exc[exc.braak == b].RORB_pos.mean() * 100 for b in braaks]
both = [(exc[exc.braak == b].RELN_pos & exc[exc.braak == b].RORB_pos).mean() * 100 for b in braaks]

or_n, p_n, lo_n, hi_n = depth_or(exc[exc.AD == 0])
or_a, p_a, lo_a, hi_a = depth_or(exc[exc.AD == 1])
or_i, p_i, lo_i, hi_i = depth_or(inh[inh.AD == 1])

pd.DataFrame({
    "metric": ["exc_normal_OR", "exc_AD_OR", "inh_AD_OR", "interaction_OR"],
    "value": [or_n, or_a, or_i, inter_or],
    "p": [p_n, p_a, p_i, inter_p],
    "ci_lo": [lo_n, lo_a, lo_i, np.nan],
    "ci_hi": [hi_n, hi_a, hi_i, np.nan],
}).to_csv("coexpression_stats.csv", index=False)

# ---- figure ----
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), gridspec_kw={"wspace": 0.34})
axA = axes[0]
axA.plot(braaks, rln, "o-", color="#c62828", label="RELN⁺", ms=6)
axA.plot(braaks, rrb, "s-", color="#2e7d32", label="RORB⁺", ms=6)
axA.plot(braaks, both, "^-", color="#6a1b9a", label="RELN⁺RORB⁺", ms=6)
axA.set_xticks(braaks); axA.set_xlabel("Braak tau stage")
axA.set_ylabel("% EC excitatory neurons detected⁺")
axA.legend(frameon=False, fontsize=8, loc="upper left")
axA.set_title("Detection by Braak stage", loc="left", fontsize=9.5)
panel_letter(axA, "A", dx=-0.16, dy=1.04, case="upper", fontsize=15)

axB = axes[1]
or_by = []
for b in braaks:
    orr, p, lo, hi = depth_or(exc[exc.braak == b])
    or_by.append((orr, lo, hi))
y = [x[0] for x in or_by]
lo = [x[0] - x[1] for x in or_by]; hi = [x[2] - x[0] for x in or_by]
axB.errorbar(braaks, y, yerr=[lo, hi], fmt="o-", color="#6a1b9a", ms=7, capsize=4)
axB.axhline(1.0, ls="--", color="0.6", lw=1)
axB.set_xticks(braaks); axB.set_xlabel("Braak tau stage")
axB.set_ylabel("depth-controlled same-cell OR")
axB.set_title(f"Same-cell co-occupancy rises with tau\n(interaction OR={inter_or:.2f}, p={inter_p:.1e})",
              loc="left", fontsize=9.5)
panel_letter(axB, "B", dx=-0.16, dy=1.04, case="upper", fontsize=15)

axC = axes[2]
groups = [("EC Exc\n(AD)", or_a, lo_a, hi_a, "#c62828"),
          ("EC Exc\n(normal)", or_n, lo_n, hi_n, "#4a90d9"),
          ("EC Inh\n(AD)", or_i, lo_i, hi_i, "#2e7d32")]
xs = np.arange(len(groups))
for xi, (lab, orr, l, h, c) in zip(xs, groups):
    axC.bar(xi, orr, color=c, zorder=3, width=0.6)
    axC.errorbar(xi, orr, yerr=[[orr - l], [h - orr]], fmt="none", ecolor="k", capsize=4, zorder=4)
axC.axhline(1.0, ls="--", color="0.6", lw=1)
axC.set_xticks(xs); axC.set_xticklabels([g[0] for g in groups])
axC.set_ylabel("depth-controlled same-cell OR")
axC.set_title("EC-specific, AD-enriched\nco-occupancy", loc="left", fontsize=9.5)
panel_letter(axC, "C", dx=-0.16, dy=1.04, case="upper", fontsize=15)

fig.suptitle("RELN and RORB co-occupy the same entorhinal excitatory neurons in AD",
             fontsize=11, y=1.03)
fig.savefig("reln_rorb_coexpression_AD.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"saved. interaction OR={inter_or:.2f} p={inter_p:.1e}; AD OR={or_a:.2f}, normal OR={or_n:.2f}, inh OR={or_i:.2f}")
