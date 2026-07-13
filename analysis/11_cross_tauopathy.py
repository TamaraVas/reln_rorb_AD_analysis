"""
Cross-tauopathy specificity (Rexach 2024)
=========================================

Per-disorder depth-controlled same-cell OR in AD, Pick, PSP, normal across BA4/insula/V1 (non-EC). Establishes tau-conditional, EC-substrate-requiring specificity.

Inputs : streamed from CELLxGENE Census dataset ac0c6561 (see fetch note)
Outputs: cross_tauopathy_specificity.png, cross_tauopathy_results.csv, cross_tauopathy_byregion.csv
Method : statsmodels Logit; scipy Fisher exact

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import fisher_exact
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"wspace": 0.34, "width_ratios": [1.15, 1]})
axA = axes[0]
labs = ["Normal", "AD\n(non-EC)", "Pick\ndisease", "PSP"]
ors = [0.943, 0.820, 0.762, 0.730]
ps = [0.426, 2.57e-3, 4.67e-4, 6.87e-8]
cols_bar = ["#9e9e9e", "#d84315", "#d84315", "#d84315"]
x = np.arange(len(labs))
axA.axhline(1.0, ls="--", color="#555", lw=1, zorder=0)
axA.axhline(2.4, ls=":", color="#6a1b9a", lw=1.6, zorder=0)
axA.text(3.35, 2.44, "EC excitatory\n(this study, OR≈2.4)", color="#6a1b9a", fontsize=8, va="bottom", ha="right")
axA.bar(x, ors, color=cols_bar, width=0.62, zorder=3)
for xi, o, p in zip(x, ors, ps):
    star = "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else "n.s."))
    axA.text(xi, o + 0.05, star, ha="center", va="bottom", fontsize=11 if p < 0.05 else 8,
             fontweight="bold" if p < 0.05 else "normal", color="k" if p < 0.05 else "#555")
axA.set_xticks(x)
axA.set_xticklabels(labs, fontsize=9)
axA.set_ylabel("RELN×RORB same-cell OR\n(depth-controlled, glutamatergic)")
axA.set_ylim(0, 2.75)
axA.set_title("No convergence in non-EC cortex —\nnot in AD, nor two other tauopathies", loc="left", fontsize=9.5)
axB = axes[1]
rorbpct = [37.7, 43.0, 46.6, 51.7]
relnpct = [3.0, 3.0, 2.3, 2.9]
xb = np.arange(len(labs))
w = 0.38
axB.bar(xb - w / 2, rorbpct, w, label="RORB⁺ %", color="#1565c0", zorder=3)
axB.bar(xb + w / 2, relnpct, w, label="RELN⁺ %", color="#c62828", zorder=3)
axB.set_xticks(xb)
axB.set_xticklabels(labs, fontsize=9)
axB.set_ylabel("% glutamatergic neurons detected⁺")
axB.legend(frameon=False, fontsize=8, loc="upper left")
axB.set_title("Substrate present (RORB⁺ 38→52%)\nyet the pair does not co-occupy", loc="left", fontsize=9.5)
panel_letter(axA, "A", dx=-0.16, dy=1.04, case="upper", fontsize=15)
panel_letter(axB, "B", dx=-0.14, dy=1.04, case="upper", fontsize=15)
fig.suptitle("Cross-tauopathy specificity (Rexach et al. 2024, 432,555 nuclei; BA4 / insula / V1)", fontsize=10.5, y=1.03)
fig.savefig("cross_tauopathy_specificity.png", dpi=200, bbox_inches="tight")
plt.close(fig)
