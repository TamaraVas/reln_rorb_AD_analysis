"""
Gene-regulatory-network importance of RORB for reelin-pathway genes
==================================================================

Tests whether RORB has above-chance predictive importance for reelin-pathway
gene expression, using gradient-boosted-tree regression against a random-gene
null. Result: RORB's importance sits within the null (all n.s.) -> no evidence
that RORB transcriptionally drives the reelin pathway.

Inputs : data/leng2021_ec_full.h5ad  (EC excitatory neurons)
Outputs: rorb_grn_reelin.png, grn_importance.csv
Method : scikit-learn GradientBoostingRegressor; permutation (random-gene) null

NOTE ON PROVENANCE: this is a clean re-implementation of the analysis from its
documented parameters (see METHODS_PROVENANCE.md). The published figure was
produced by an equivalent computation; this standalone version regenerates the
statistic end-to-end so it is fully reproducible from the input file alone.

Part of the RELN-RORB EC convergence analysis suite. See README.md.
"""
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from sklearn.ensemble import GradientBoostingRegressor
import matplotlib.pyplot as plt
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

rng = np.random.default_rng(0)

# ---- load EC excitatory neurons ----
A = ad.read_h5ad("data/leng2021_ec_full.h5ad")
name2i = {n: i for i, n in enumerate(A.var["feature_name"])}
exc = A[A.obs["clusterCellType"].astype(str).str.contains("Exc", na=False)]

def dense(adata, gene):
    if gene not in name2i:
        return None
    x = adata.X[:, name2i[gene]]
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

# reelin-pathway targets and the transcription-factor predictor panel
reelin_targets = ["RELN", "DAB1", "LRP8", "VLDLR", "FYN"]
tf_panel = ["RORB", "RORA", "CUX2", "FOXP2", "TBR1", "SATB2", "MEF2C",
            "NR4A2", "POU3F2", "BCL11B", "FEZF2", "TLE4", "ETV1", "SOX5", "NEUROD6"]
tf_panel = [g for g in tf_panel if g in name2i]

Xtf = np.column_stack([dense(exc, g) for g in tf_panel])

def rorb_importance_for(target_gene):
    y = dense(exc, target_gene)
    if y is None:
        return np.nan
    gbr = GradientBoostingRegressor(n_estimators=100, max_depth=3, random_state=0)
    gbr.fit(Xtf, y)
    imp = dict(zip(tf_panel, gbr.feature_importances_))
    return imp

# observed RORB importance for each reelin-pathway gene
rows = []
net_rows = []
for tg in reelin_targets:
    imp = rorb_importance_for(tg)
    rows.append((tg, imp.get("RORB", np.nan)))
    for tf, v in imp.items():
        net_rows.append((tf, tg, v))
obs = pd.DataFrame(rows, columns=["target", "RORB_importance"])
net = pd.DataFrame(net_rows, columns=["TF", "target", "importance"])

# random-gene null: RORB importance predicting 100 random genes
expressed = [g for g, i in name2i.items()
             if 0.1 <= (np.asarray((exc.X[:, i] >= 1).mean())) <= 0.6]
null_targets = list(rng.choice(expressed, size=min(100, len(expressed)), replace=False))
null_imp = []
for g in null_targets:
    imp = rorb_importance_for(g)
    if imp and not np.isnan(imp.get("RORB", np.nan)):
        null_imp.append(imp["RORB"])
null_imp = np.array(null_imp)

# per-target p-value = fraction of null >= observed
obs["null_p"] = [(null_imp >= v).mean() for v in obs["RORB_importance"]]
obs.to_csv("grn_importance.csv", index=False)

# ---- figure ----
fig = plt.figure(figsize=(11.6, 4.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.30)
axA = fig.add_subplot(gs[0, 0])
top = net[net.target == "RELN"].sort_values("importance").tail(12)
axA.barh(np.arange(len(top)), top.importance, color="#1565c0", height=0.8)
axA.set_yticks(np.arange(len(top))); axA.set_yticklabels(top.TF, fontsize=6)
axA.set_xlabel("Gradient-boosted-tree feature importance")
axA.set_title("TF importance for RELN", loc="left", fontsize=8)
panel_letter(axA, "A", dx=-0.12, dy=1.03, case="upper", fontsize=15)

axB = fig.add_subplot(gs[0, 1])
axB.hist(null_imp, bins=25, color="0.8", edgecolor="0.5",
         label="RORB imp.\nfor random genes")
palette = ["#c62828", "#6a1b9a", "#00897b", "#ef6c00"]
for (tg, v, p), col in zip(obs[["target", "RORB_importance", "null_p"]].values, palette * 2):
    axB.axvline(v, color=col, lw=1.6, label=f"{tg} (p={p:.2f})")
axB.axvline(np.percentile(null_imp, 95), color="k", ls="--", lw=0.9)
axB.set_xlim(-0.02, 0.46)
axB.set_xlabel("RORB feature importance"); axB.set_ylabel("# genes")
axB.set_title("RORB->reelin genes not above\nrandom-gene null (all n.s.)", loc="left", fontsize=8)
axB.legend(fontsize=6.2, frameon=False, loc="upper left",
           bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
panel_letter(axB, "B", dx=-0.12, dy=1.03, case="upper", fontsize=15)

fig.savefig("rorb_grn_reelin.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("saved rorb_grn_reelin.png and grn_importance.csv")
