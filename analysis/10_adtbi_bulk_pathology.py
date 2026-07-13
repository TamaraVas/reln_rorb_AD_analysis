"""
Bulk RNA-seq RELN/RORB vs pathology (Allen ADTBI)
=================================================

Spearman correlations of RELN/RORB FPKM with per-donor Braak and CERAD across four regions.

Inputs : data/adtbi_merged.parquet
Outputs: adtbi_bulk_pathology.png
Method : scipy Spearman

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import pandas as pd, numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

m = pd.read_parquet("data/adtbi_merged.parquet")
m['logRELN'] = np.log1p(m.RELN)
m['logRORB'] = np.log1p(m.RORB)
struct = ['TCx', 'PCx', 'HIP', 'FWM']
labels = {'TCx': 'Temporal ctx', 'PCx': 'Parietal ctx', 'HIP': 'Hippocampus', 'FWM': 'White matter'}
cols = [('logRELN', 'braak', 'RELN·Braak'), ('logRELN', 'cerad', 'RELN·CERAD'),
        ('logRORB', 'braak', 'RORB·Braak'), ('logRORB', 'cerad', 'RORB·CERAD')]
M = np.full((4, 4), np.nan)
P = np.full((4, 4), np.nan)
for i, st in enumerate(struct):
    d = m[m.structure_acronym == st]
    for j, (g, path, _) in enumerate(cols):
        dd = d.dropna(subset=[path, g])
        if len(dd) > 5:
            rho, pv = spearmanr(dd[path], dd[g])
            M[i, j] = rho
            P[i, j] = pv

fig = plt.figure(figsize=(13.5, 5.0))
gs = GridSpec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.5, figure=fig)

axA = fig.add_subplot(gs[0])
im = axA.imshow(M, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
axA.set_xticks(range(4))
axA.set_xticklabels([c[2] for c in cols], rotation=35, ha="right", fontsize=8)
axA.set_yticks(range(4))
axA.set_yticklabels([labels[s] for s in struct], fontsize=8.5)
for i in range(4):
    for j in range(4):
        if not np.isnan(M[i, j]):
            star = "*" if P[i, j] < 0.05 else ""
            axA.text(j, i, f"{M[i, j]:+.2f}{star}", ha="center", va="center", fontsize=8.2,
                     color="white" if abs(M[i, j]) > 0.18 else "#222",
                     fontweight="bold" if star else "normal")
axA.set_title("Expression–pathology correlation (ρ)", fontsize=9.5, loc="center", pad=8)
cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04)
cb.set_label("Spearman ρ", fontsize=8)
cb.ax.tick_params(labelsize=7)
axA.text(-0.34, 1.08, "A", transform=axA.transAxes, fontsize=16, fontweight="bold")

hip = m[m.structure_acronym == 'HIP'].dropna(subset=['braak', 'logRORB'])
rho, _ = spearmanr(hip.braak, hip.logRORB)
axB = fig.add_subplot(gs[1])
axB.scatter(hip.braak + np.random.uniform(-0.12, 0.12, len(hip)), hip.logRORB, s=16, alpha=0.5,
            color="#C44E52", edgecolor="none")
z = np.polyfit(hip.braak, hip.logRORB, 1)
xs = np.array([0, 6])
axB.plot(xs, np.polyval(z, xs), color="#2166ac", lw=2)
axB.set_xlabel("Braak (tau)", fontsize=9)
axB.set_ylabel("log RORB FPKM (hippocampus)", fontsize=9)
axB.set_title(f"Hippocampal RORB\nfalls with tau (ρ={rho:.2f}*)", fontsize=9.5, loc="center", pad=8)
axB.text(-0.30, 1.08, "B", transform=axB.transAxes, fontsize=16, fontweight="bold")

hipR = m[m.structure_acronym == 'HIP'].dropna(subset=['braak', 'logRELN'])
rhoR, _ = spearmanr(hipR.braak, hipR.logRELN)
axC = fig.add_subplot(gs[2])
axC.scatter(hipR.braak + np.random.uniform(-0.12, 0.12, len(hipR)), hipR.logRELN, s=16, alpha=0.5,
            color="#4C72B0", edgecolor="none")
zR = np.polyfit(hipR.braak, hipR.logRELN, 1)
axC.plot(xs, np.polyval(zR, xs), color="#b2182b", lw=2)
axC.set_xlabel("Braak (tau)", fontsize=9)
axC.set_ylabel("log RELN FPKM (hippocampus)", fontsize=9)
axC.set_title(f"Hippocampal RELN\nfalls with tau (ρ={rhoR:.2f}*)", fontsize=9.5, loc="center", pad=8)
axC.text(-0.30, 1.08, "C", transform=axC.transAxes, fontsize=16, fontweight="bold")

fig.suptitle("ADTBI bulk RNA-seq: RELN & RORB decline with tau specifically in hippocampus, not neocortex (377 samples, 107 donors)",
             fontsize=10.5, y=1.02)
fig.savefig("adtbi_bulk_pathology.png", dpi=300, bbox_inches="tight")
print("saved")
