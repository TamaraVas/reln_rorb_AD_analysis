"""
Region × cell-class specificity of the convergence
==================================================

Fits the interaction model within each cell class across regions to show the convergence is exclusive to excitatory neurons.

Inputs : data/cellclass_all_regions.csv
Outputs: reln_rorb_region_class_specificity.png
Method : precomputed ORs (from script 04); matplotlib

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

cc = pd.read_csv("data/cellclass_all_regions.csv")

classes = ["Exc", "Inh", "Astrocyte", "Oligodendrocyte", "Immune", "OPC"]
class_lbl = {"Exc": "Excitatory", "Inh": "Inhibitory", "Astrocyte": "Astrocyte",
             "Oligodendrocyte": "Oligodendrocyte", "Immune": "Microglia", "OPC": "OPC"}
regions = ["MEC", "LEC"]
region_col = {"MEC": "#6a1b9a", "LEC": "#9c6ade"}

fig, ax = plt.subplots(figsize=(9.5, 5.0))
nreg = len(regions); w = 0.36
x = np.arange(len(classes))
for ri, reg in enumerate(regions):
    xpos = x + (ri - (nreg - 1) / 2) * w
    for xi, cl in enumerate(classes):
        row = cc[(cc.region == reg) & (cc.cellclass == cl)]
        if not len(row) or pd.isna(row.interaction_OR.values[0]):
            continue
        v = row.interaction_OR.values[0]; s = row.sig.values[0]
        ax.bar(xpos[xi], v, w, color=region_col[reg], edgecolor="#333", linewidth=0.5)
        if isinstance(s, str) and "*" in s:
            ax.text(xpos[xi], v + 0.03, s, ha="center", va="bottom", fontsize=13, fontweight="bold", color="black")

ax.axhline(1.0, ls="--", color="#888", lw=1, zorder=0)
ax.set_xticks(x); ax.set_xticklabels([class_lbl[c] for c in classes], fontsize=9.5)
ax.set_ylabel("RELN×RORB×Braak same-cell interaction OR / step", fontsize=9.5)
ax.set_title("Convergence is exclusive to excitatory neurons (entorhinal cortex)", fontsize=10.5, loc="left", pad=10)
ax.set_ylim(0, 2.0)
ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
handles = [Patch(facecolor=region_col[r], edgecolor="#333", label=f"{r} (entorhinal)") for r in regions]
ax.legend(handles=handles, loc="upper right", fontsize=9, title="Region", title_fontsize=9)
fig.savefig("reln_rorb_region_class_specificity.png", dpi=300, bbox_inches="tight")
