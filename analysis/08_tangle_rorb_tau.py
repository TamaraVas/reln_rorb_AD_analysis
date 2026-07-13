"""
RORB in tangle-bearing neurons (Otero-Garcia)
=============================================

Within-donor comparison of RORB/RELN detection in AT8+ tangle-bearing vs tangle-free prefrontal neurons.

Inputs : data/oterogarcia_tangle_cells.parquet
Outputs: tangle_rorb_tau_amyloid.png
Method : statsmodels Logit with donor fixed effects

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
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import fisher_exact
import anndata as ad
import matplotlib.pyplot as plt
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

tang = pd.read_parquet("data/oterogarcia_tangle_cells.parquet").reset_index(drop=True)


tang["tau"] = (tang.SORT == "AT8").astype(int)
amap = {"No": 0, "DP (C0)": 1, "C3": 2}
tang["amyloid"] = tang.Amyloid.map(amap)
tang["rorb_pos"] = tang.RORB_pos.astype(int)
tang["reln_pos"] = tang.RELN_pos.astype(int)
tang["logn"] = np.log(tang.nCount_RNA.clip(lower=1))

# Within-donor tau effect (donor fixed effects)
donor_braak = tang.groupby("donor", observed=True)["tau"].nunique() if hasattr(tang, "donor") else pd.Series(dtype=int)

# Identify AD donors with both tangle-bearing and tangle-free cells
try:
    ad_donors = tang.groupby("donor", observed=True).tau.nunique()
    ad_donors = ad_donors[ad_donors > 1].index.tolist()
except Exception:
    ad_donors = []

w = tang[tang.donor.isin(ad_donors)].copy() if ad_donors else tang.copy()

# Build figure: tangle RORB + RELN associations
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: RORB+ % by tangle status
axA = axes[0]
sort_order = tang.SORT.unique()
rorb_by_sort = tang.groupby("SORT")["rorb_pos"].mean() * 100
reln_by_sort = tang.groupby("SORT")["reln_pos"].mean() * 100

for yi, (gene, vals, col) in enumerate([("RORB⁺", rorb_by_sort, "#C44E52"), ("RELN⁺", reln_by_sort, "#4C72B0")]):
    groups = sorted(vals.index, key=lambda x: 0 if x == "MAP2" else 1)
    xvals = [vals.get(g, 0) for g in groups]
    axA.barh([yi * 3 + gi for gi in range(len(groups))], xvals, color=col, alpha=0.7 if gi == 0 else 0.4,
             label=gene if gi == 0 else "")

axA.set_xlabel("% cells positive")
axA.set_title("A  RORB and RELN by tangle status\n(Otero-Garcia prefrontal cortex)", loc="left")

# Simpler version: bar chart of RORB+ % by SORT group
axA.cla()
groups_sorted = ["MAP2", "AT8"]
rorb_vals = [tang[tang.SORT == g]["rorb_pos"].mean() * 100 for g in groups_sorted if g in tang.SORT.values]
reln_vals = [tang[tang.SORT == g]["reln_pos"].mean() * 100 for g in groups_sorted if g in tang.SORT.values]

x = np.arange(len(groups_sorted))
w_bar = 0.35
axA.bar(x - w_bar/2, rorb_vals, w_bar, label="RORB⁺", color="#C44E52")
axA.bar(x + w_bar/2, reln_vals, w_bar, label="RELN⁺", color="#4C72B0")
axA.set_xticks(x)
axA.set_xticklabels(["Tangle-free\n(MAP2)", "Tangle-bearing\n(AT8)"])
axA.set_ylabel("% positive cells")
axA.legend(fontsize=8)
axA.set_title("A  RORB and RELN fall in\ntangle-bearing neurons (PFC)", loc="left", fontsize=9)

# Panel B: within-donor ORs
axB = axes[1]
results = []
for gene, label, col in [("rorb_pos", "RORB⁺", "#C44E52"), ("reln_pos", "RELN⁺", "#4C72B0")]:
    if len(w) > 100 and w[gene].sum() > 10:
        try:
            D = pd.get_dummies(w.donor, prefix="d", drop_first=True).astype(float)
            import statsmodels.api as sm2
            Xfe = pd.concat([w[["tau", "logn"]].astype(float).reset_index(drop=True),
                             D.reset_index(drop=True)], axis=1)
            Xfe = sm2.add_constant(Xfe)
            mfe = sm2.Logit(w[gene].astype(float).values, Xfe.values).fit(disp=0)
            ti = list(Xfe.columns).index("tau")
            orr = np.exp(mfe.params[ti])
            lo, hi = np.exp(mfe.conf_int()[ti])
            p = mfe.pvalues[ti]
            results.append((label, orr, lo, hi, p, col))
        except Exception:
            results.append((label, np.nan, np.nan, np.nan, np.nan, col))
    else:
        results.append((label, np.nan, np.nan, np.nan, np.nan, col))

axB.axvline(1, color="#888", ls="--", lw=1)
for yi, (label, orr, lo, hi, p, col) in enumerate(results):
    if not np.isnan(orr):
        axB.errorbar(orr, yi, xerr=[[orr - lo], [hi - orr]], fmt="o", color=col, ms=8, capsize=4, lw=2)
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        axB.annotate(f"OR {orr:.2f} ({star})", (hi, yi), xytext=(5, 0),
                     textcoords="offset points", va="center", fontsize=8, color=col)

axB.set_yticks(range(len(results)))
axB.set_yticklabels([r[0] for r in results])
axB.set_xlabel("Within-donor OR (tangle-bearing vs tangle-free)")
axB.set_title("B  Within-donor tau effect\n(donor fixed effects, 8 AD donors)", loc="left", fontsize=9)
axB.set_xlim(0.5, 1.5)

fig.suptitle("RORB and RELN are reduced in tangle-bearing neurons (Otero-Garcia PFC, within-donor)",
             fontsize=10, y=1.02)
fig.tight_layout()
fig.savefig("tangle_rorb_tau_amyloid.png", dpi=300, bbox_inches="tight")
print("saved tangle_rorb_tau_amyloid.png")
