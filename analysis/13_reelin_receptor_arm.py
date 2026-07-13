"""
Reelin receiver machinery in RORB+/double-positive neurons
==========================================================

Receptor/transducer arm (LRP8/VLDLR/DAB1/FYN) detection by cell group and its AD shift within RORB+ receivers.

Inputs : data/leng2021_ec_full.h5ad
Outputs: reelin_receptor_arm.png, reelin_receptor_bygroup.csv, reelin_receptor_AD_shift.csv
Method : statsmodels Logit

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib as mpl
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

A = ad.read_h5ad("data/leng2021_ec_full.h5ad")

name2i = {n: i for i, n in enumerate(A.var["feature_name"])}


def vec(gene):
    if gene not in name2i:
        return None
    x = A.X[:, name2i[gene]]
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()


obs = A.obs
exc_mask = obs["clusterCellType"].astype(str).str.contains("Exc", case=False, na=False)
Eidx = np.where(exc_mask.values)[0]
logdepth = np.log1p(obs.loc[exc_mask, "nUMI"].astype(float).values)
isAD = (obs.loc[exc_mask, "disease"] == "Alzheimer disease").values


def pos(gene):
    v = vec(gene)
    return None if v is None else (v[Eidx] >= 1)


rorb_pos = pos("RORB")
reln_pos = pos("RELN")
dp = reln_pos & rorb_pos

receptors = ["LRP8", "VLDLR", "DAB1", "FYN"]

recrows = []
for g in receptors:
    gp = pos(g)
    if gp is None:
        continue
    r_neg = 100 * gp[~rorb_pos].mean()
    r_pos = 100 * gp[rorb_pos].mean()
    r_dp = 100 * gp[dp].mean()
    enr = r_pos / max(r_neg, 1e-9)
    recrows.append((g, round(r_neg, 1), round(r_pos, 1), round(r_dp, 1), round(enr, 2)))

adrows = []
for g in receptors:
    gp = pos(g)
    sub = rorb_pos
    y = gp[sub].astype(int)
    Xd = sm.add_constant(np.column_stack([isAD[sub].astype(int), logdepth[sub]]))
    m = sm.Logit(y, Xd).fit(disp=0)
    orv, pv = np.exp(m.params[1]), m.pvalues[1]
    n_pct = 100 * gp[sub & ~isAD].mean()
    a_pct = 100 * gp[sub & isAD].mean()
    adrows.append((g, round(n_pct, 1), round(a_pct, 1), round(orv, 2), pv))

apply_figure_style()

fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw={"wspace": 0.30, "width_ratios": [1, 1]})
axA = axes[0]
x = np.arange(len(receptors))
w = 0.26
neg = [r[1] for r in recrows]
posr = [r[2] for r in recrows]
dpr = [r[3] for r in recrows]
axA.bar(x - w, neg, w, label="RORB⁻", color="#bdbdbd", zorder=3)
axA.bar(x, posr, w, label="RORB⁺", color="#1565c0", zorder=3)
axA.bar(x + w, dpr, w, label="RELN⁺RORB⁺ (DP)", color="#6a1b9a", zorder=3)
axA.set_xticks(x)
axA.set_xticklabels(receptors)
axA.set_ylabel("% EC excitatory neurons detected⁺")
axA.set_ylim(0, 108)
axA.legend(frameon=False, fontsize=7.8, loc="upper left", bbox_to_anchor=(0.0, 1.0))
axA.set_title("Reelin receptor/transducer arm is\nenriched in RORB⁺ and DP neurons", loc="left", fontsize=9.5)
panel_letter(axA, "A", dx=-0.15, dy=1.04, case="upper", fontsize=15)

axB = axes[1]
nn = [r[1] for r in adrows]
aa = [r[2] for r in adrows]
axB.bar(x - 0.2, nn, 0.4, label="Normal", color="#4a90d9", zorder=3)
axB.bar(x + 0.2, aa, 0.4, label="AD", color="#c62828", zorder=3)
axB.set_ylim(0, 98)
for xi, (g, n_, a_, orv, pv) in zip(x, adrows):
    star = "***" if pv < 1e-3 else ("**" if pv < 1e-2 else ("*" if pv < 0.05 else "n.s."))
    axB.text(xi, max(n_, a_) + 2.5, star, ha="center", va="bottom",
             fontsize=9 if pv < 0.05 else 7,
             fontweight="bold" if pv < 0.05 else "normal",
             color="k" if pv < 0.05 else "#555")
axB.set_xticks(x)
axB.set_xticklabels(receptors)
axB.set_ylabel("% RORB⁺ neurons detected⁺")
axB.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 1.0))
axB.set_title("Receptor competence rises in AD\n(within RORB⁺ receivers, depth-controlled)", loc="left", fontsize=9.5)
panel_letter(axB, "B", dx=-0.15, dy=1.04, case="upper", fontsize=15)
fig.suptitle("Reelin receiver machinery in the RORB⁺ / double-positive population (Leng 2021 EC, n=10,780)", fontsize=10.5, y=1.02)
fig.savefig("reelin_receptor_arm.png", dpi=200, bbox_inches="tight")
plt.close(fig)
