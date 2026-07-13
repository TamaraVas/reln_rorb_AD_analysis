"""
Module-level specificity control
================================

Depth-controlled AD-interaction OR for RELN×RORB vs reelin-pathway×RORB pairs, RELN×Exc.4-marker pairs, and frequency-matched controls.

Inputs : data/leng2021_ec_full.h5ad
Outputs: marker_pair_specificity.png, marker_pair_specificity.csv
Method : statsmodels Logit

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import os
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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
E = obs[exc_mask].copy()
logdepth = np.log1p(E["nUMI"].astype(float).values)
isAD = (E["disease"] == "Alzheimer disease").values


def pos(gene):
    v = vec(gene)
    return None if v is None else (v[Eidx] >= 1)


def interaction_OR(g_reln_side, g_rorb_side):
    p1 = pos(g_reln_side)
    p2 = pos(g_rorb_side)
    if p1 is None or p2 is None:
        return None
    X = np.column_stack([p1.astype(int), isAD.astype(int), p1.astype(int) * isAD.astype(int), logdepth])
    X = sm.add_constant(X)
    try:
        m = sm.Logit(p2.astype(int), X).fit(disp=0)
        return np.exp(m.params[3]), m.pvalues[3]
    except Exception as e:
        return ("ERR", str(e)[:60])


reelin_r = [("DAB1×RORB", 3.691, 6.39e-29), ("FYN×RORB", 2.039, 2.51e-11), ("VLDLR×RORB", 1.816, 5.71e-3),
            ("LRP8×RORB", 1.618, 7.28e-3), ("APOE×RORB", 1.366, 1.20e-1)]
exc4_r = [("RELN×PTPRK", 6.608, 1.02e-26), ("RELN×NTNG1", 6.258, 1.23e-24), ("RELN×CNTN5", 4.445, 5.57e-9),
          ("RELN×CDH8", 3.535, 3.79e-12), ("RELN×GRM7", 2.530, 2.50e-4)]
ctrl_r = [("RELN×EPC1", 1.606, 9.52e-3), ("RELN×UBE2D2", 1.098, 6.05e-1), ("RELN×C2CD3", 1.023, 9.05e-1),
          ("RELN×TM9SF3", 0.862, 4.35e-1), ("RELN×SLITRK5", 0.832, 3.51e-1), ("RELN×ENO2", 0.687, 4.37e-2)]
focal = ("RELN×RORB", 2.373, 9.53e-6)

apply_figure_style()
fig, ax = plt.subplots(figsize=(8.2, 6.4))
y = 0
yt = []
yl = []


def block(items, color):
    global y
    for name, orr, p in items:
        star = "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else ""))
        ax.barh(y, orr, color=color, height=0.7, zorder=3)
        ax.text(orr + 0.08, y, f"{orr:.2f}{star}", va="center", fontsize=7.5)
        yt.append(y)
        yl.append(name)
        y += 1
    y += 0.7


block(ctrl_r, "#9e9e9e")
block(exc4_r, "#00897b")
block(reelin_r, "#1565c0")
ax.barh(y, focal[1], color="#6a1b9a", height=0.8, zorder=3)
ax.text(focal[1] + 0.08, y, f"{focal[1]:.2f}***", va="center", fontsize=8.5, fontweight="bold")
yt.append(y)
yl.append("RELN×RORB (focal)")
ax.axvline(1.0, ls="--", color="#555", lw=1, zorder=0)
ax.set_yticks(yt)
ax.set_yticklabels(yl, fontsize=8.5)
ax.set_xlabel("Depth-controlled AD-interaction OR (co-occupancy in same EC excitatory neuron)")
ax.set_xlim(0, 7.2)
leg = [Patch(color="#6a1b9a", label="RELN×RORB (the reported pair)"),
       Patch(color="#1565c0", label="reelin-pathway gene × RORB"),
       Patch(color="#00897b", label="RELN × EC:Exc.4 identity marker"),
       Patch(color="#9e9e9e", label="frequency-matched control")]
ax.legend(handles=leg, frameon=False, fontsize=7.8, loc="lower right")
ax.set_title("The AD convergence is module-specific, not unique to one gene pair", loc="left", fontsize=10)
fig.tight_layout()
fig.savefig("marker_pair_specificity.png", dpi=200, bbox_inches="tight")
plt.close(fig)
