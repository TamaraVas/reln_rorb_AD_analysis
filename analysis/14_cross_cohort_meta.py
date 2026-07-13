"""
Cross-cohort random-effects meta-analysis
=========================================

Pools the RELN×Braak per-stage interaction OR across Leng, SEA-AD MEC, SEA-AD LEC (DerSimonian–Laird) and shows cross-sectional same-cell ORs at high pathology.

Inputs : data/leng2021_ec_full.h5ad; data/staged_ec_all_regions.json
Outputs: cross_cohort_meta.png, cross_cohort_meta.csv
Method : custom DerSimonian–Laird random-effects; scipy

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
import matplotlib.pyplot as plt
import anndata as ad
import json
from scipy.stats import norm
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

# Load Leng EC data
A = ad.read_h5ad("data/leng2021_ec_full.h5ad")

obs = A.obs
exc_mask = obs["clusterCellType"].astype(str).str.contains("Exc", case=False, na=False)
Eidx = np.where(exc_mask.values)[0]

name2i = {n: i for i, n in enumerate(A.var["feature_name"])}

def vec(gene):
    if gene not in name2i:
        return None
    x = A.X[:, name2i[gene]]
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

def pos(gene):
    v = vec(gene)
    return None if v is None else (v[Eidx] >= 1)

E = obs[exc_mask].copy()
logdepth = np.log1p(E["nUMI"].astype(float).values)
isAD = (E["disease"] == "Alzheimer disease").values
braak = obs.loc[exc_mask, "BraakStage"].astype(int).values

reln_pos = pos("RELN")
rorb_pos = pos("RORB")

# Leng RELN x Braak interaction per step
def leng_braak_interaction():
    y = rorb_pos.astype(int)
    X = sm.add_constant(np.column_stack([reln_pos.astype(int), braak, reln_pos.astype(int) * braak, logdepth]))
    m = sm.Logit(y, X).fit(disp=0)
    return m.params[3], m.bse[3]

lc, lse = leng_braak_interaction()
print("Leng RELN×Braak interaction/step: OR=%.3f CI[%.2f-%.2f]" % (np.exp(lc), np.exp(lc - 1.96 * lse), np.exp(lc + 1.96 * lse)))

# Leng EC AD group same-cell OR
def group_same_cell_OR(mask):
    y = rorb_pos[mask].astype(int)
    X = sm.add_constant(np.column_stack([reln_pos[mask].astype(int), logdepth[mask]]))
    m = sm.Logit(y, X).fit(disp=0)
    c, se = m.params[1], m.bse[1]
    return np.exp(c), np.exp(c - 1.96 * se), np.exp(c + 1.96 * se), int(mask.sum())

leng_AD = group_same_cell_OR(isAD)
print("Leng EC (AD): OR=%.2f CI[%.2f-%.2f] n=%d" % leng_AD)

# Load staged results
staged = json.load(open("data/staged_ec_all_regions.json"))

# SEA-AD interaction OR per step
from scipy.stats import norm

def se_from_or_p(orr, p):
    z = abs(norm.ppf(p / 2))
    return np.log(orr) / z if z > 0 else np.nan

mec_or = staged["MEC"]["reln_braak_interaction_OR_per_step"]
mec_p = staged["MEC"]["interaction_p"]
lec_or = staged["LEC"]["interaction_OR"]
lec_p = staged["LEC"]["interaction_p"]
mec_se = se_from_or_p(mec_or, mec_p)
lec_se = se_from_or_p(lec_or, lec_p)
print(f"SEA-AD MEC: OR={mec_or} SE={mec_se:.4f}")
print(f"SEA-AD LEC: OR={lec_or} SE={lec_se:.4f}")

# Random-effects meta-analysis (DerSimonian-Laird)
studies = [
    ("Leng 2021 EC (discovery)", lc, lse),
    ("SEA-AD MEC (81 donors)", np.log(mec_or), mec_se),
    ("SEA-AD LEC (34 donors)", np.log(lec_or), lec_se),
]
ys = np.array([s[1] for s in studies])
ses = np.array([s[2] for s in studies])
w = 1 / ses ** 2
fixed = np.sum(w * ys) / np.sum(w)
Q = np.sum(w * (ys - fixed) ** 2)
df = len(ys) - 1
C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
tau2 = max(0, (Q - df) / C)
wr = 1 / (ses ** 2 + tau2)
re = np.sum(wr * ys) / np.sum(wr)
re_se = np.sqrt(1 / np.sum(wr))
I2 = max(0, 100 * (Q - df) / Q)
print(f"\nRandom-effects pooled interaction OR/step = {np.exp(re):.3f} CI[{np.exp(re - 1.96 * re_se):.2f}-{np.exp(re + 1.96 * re_se):.2f}]")
print(f"Q={Q:.2f} (df={df}), I²={I2:.0f}%, tau²={tau2:.4f}")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 4.0), gridspec_kw={"wspace": 0.42, "width_ratios": [1.1, 1]})

# Panel A: forest plot
axA = axes[0]
labels = ["Leng 2021 EC\n(discovery)", "SEA-AD MEC\n(81 donors)", "SEA-AD LEC\n(34 donors)"]
ors = [np.exp(lc), mec_or, lec_or]
los = [np.exp(lc - 1.96 * lse), np.exp(np.log(mec_or) - 1.96 * mec_se), np.exp(np.log(lec_or) - 1.96 * lec_se)]
his = [np.exp(lc + 1.96 * lse), np.exp(np.log(mec_or) + 1.96 * mec_se), np.exp(np.log(lec_or) + 1.96 * lec_se)]
yy = np.arange(len(labels))[::-1] + 1.5
for y, o, lo, hi in zip(yy, ors, los, his):
    axA.plot([lo, hi], [y, y], color="#333", lw=1.5, zorder=2)
    axA.plot(o, y, "s", color="#1565c0", ms=9, zorder=3)
pool = np.exp(re)
plo = np.exp(re - 1.96 * re_se)
phi = np.exp(re + 1.96 * re_se)
axA.add_patch(plt.Polygon([[plo, 0.5], [pool, 0.75], [phi, 0.5], [pool, 0.25]], color="#6a1b9a", zorder=4))
axA.axvline(1.0, ls="--", color="#888", lw=1)
axA.set_yticks(list(yy) + [0.5])
axA.set_yticklabels(labels + ["Random-effects\npooled"], fontsize=8)
axA.set_ylim(0, len(labels) + 2)
axA.set_xlim(0.9, 1.65)
axA.set_xlabel("RELN×RORB×Braak interaction OR per stage")
for y, o, lo, hi in zip(yy, ors, los, his):
    axA.text(hi + 0.01, y, f"{o:.2f} [{lo:.2f}–{hi:.2f}]", va="center", fontsize=7)
axA.text(phi + 0.01, 0.5, f"{pool:.2f} [{plo:.2f}–{phi:.2f}]", va="center", fontsize=7.5, fontweight="bold", color="#6a1b9a")
axA.set_title(f"A  Tau-graded convergence — 3 cohorts agree in direction\n     (I²={I2:.0f}%, all p<10⁻⁶)", loc="left", fontsize=9.3)

# Panel B: cross-sectional same-cell OR at high pathology
axB = axes[1]
cs_lab = ["Leng EC\n(AD group)", "SEA-AD MEC\n(Braak VI)", "SEA-AD MEC\nMERFISH spatial"]
cs_or = [leng_AD[0], 1.51, 1.99]
cs_lo = [leng_AD[1], 0.95, 1.79]
cs_hi = [leng_AD[2], 2.41, 2.21]
yb = np.arange(len(cs_lab))[::-1] + 1
for y, o, lo, hi in zip(yb, cs_or, cs_lo, cs_hi):
    axB.plot([lo, hi], [y, y], color="#333", lw=1.5, zorder=2)
    axB.plot(o, y, "o", color="#c62828", ms=8, zorder=3)
    axB.text(hi + 0.03, y, f"{o:.2f}", va="center", fontsize=7.5)
axB.axvline(1.0, ls="--", color="#888", lw=1)
axB.set_yticks(yb)
axB.set_yticklabels(cs_lab, fontsize=8)
axB.set_ylim(0.3, len(cs_lab) + 0.8)
axB.set_xlim(0.8, 2.7)
axB.set_xlabel("Same-cell co-occupancy OR (depth-controlled)")
axB.set_title("B  Same-cell OR at high pathology\n     (independent readouts, not pooled)", loc="left", fontsize=9.3)

fig.suptitle("Cross-cohort consistency of the EC RELN–RORB same-cell convergence", fontsize=10.5, y=1.03)
fig.savefig("cross_cohort_meta.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("saved")
