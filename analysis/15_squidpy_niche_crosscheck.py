"""
Spatial niche cross-check with squidpy (standard-tool validation)
=================================================================

Independent corroboration of the custom kNN neighborhood test in
07_mec_spatial_colocalization.py, using squidpy's peer-reviewed
neighborhood-enrichment permutation test (Palla et al. 2022, Nat Methods).

Confirms: RELN+RORB+ double-positive (DP) neurons spatially self-cluster
(z ~ 42) and are positively associated with RORB+ neighbors, while excluded
from marker-negative regions.

Inputs : data/seaad_merfish_mec.h5ad  (SEA-AD MERFISH MEC, 433-gene panel)
Outputs: squidpy_niche_crosscheck.png, squidpy_nhood_zscores.csv
Method : squidpy spatial_neighbors (Delaunay) + nhood_enrichment (1000 perms)

Part of the RELN-RORB EC convergence analysis suite. See README.md.
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/nb")
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import warnings; warnings.filterwarnings("ignore")

A = ad.read_h5ad("data/seaad_merfish_mec.h5ad")
gi = {g: i for i, g in enumerate(A.var_names)}

def gv(g):
    x = A.X[:, gi[g]]
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

reln = gv("RELN") >= 1
rorb = gv("RORB") >= 1
glut = (A.obs["Class_scANVI"].astype(str) == "Neuronal: Glutamatergic").values
G = A[glut].copy()
r, b = reln[glut], rorb[glut]
lab = np.where(r & b, "RELN+RORB+ (DP)",
       np.where(r & ~b, "RELN+ only",
        np.where(~r & b, "RORB+ only", "neither")))
G.obs["copair"] = pd.Categorical(lab, categories=["RELN+RORB+ (DP)", "RELN+ only", "RORB+ only", "neither"])

# spatial graph + neighborhood enrichment (standard squidpy pipeline)
sq.gr.spatial_neighbors(G, coord_type="generic", delaunay=True)
sq.gr.nhood_enrichment(G, cluster_key="copair", seed=0, n_perms=1000, show_progress_bar=False)

cats = list(G.obs["copair"].cat.categories)
zdf = pd.DataFrame(G.uns["copair_nhood_enrichment"]["zscore"], index=cats, columns=cats)
zdf.to_csv("squidpy_nhood_zscores.csv")

fig, ax = plt.subplots(figsize=(6.2, 5.2))
norm = TwoSlopeNorm(vmin=-55, vcenter=0, vmax=55)
im = ax.imshow(zdf.values, cmap="RdBu_r", norm=norm)
short = ["DP", "RELN+", "RORB+", "neither"]
ax.set_xticks(range(4)); ax.set_yticks(range(4))
ax.set_xticklabels(short); ax.set_yticklabels(short)
for i in range(4):
    for j in range(4):
        v = zdf.values[i, j]
        ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                color="white" if abs(v) > 28 else "black", fontweight="bold")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("neighborhood-enrichment z-score")
ax.set_title("Squidpy neighborhood enrichment (MERFISH MEC, 1000 perms)", loc="left", fontsize=9.5)
fig.tight_layout()
fig.savefig("squidpy_niche_crosscheck.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("DP self-enrichment z =", round(zdf.loc["RELN+RORB+ (DP)", "RELN+RORB+ (DP)"], 1))
