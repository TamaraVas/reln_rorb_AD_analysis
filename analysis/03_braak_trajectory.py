"""
Braak-graded trajectory of the same-cell convergence (Leng)
===========================================================

Per-Braak-stage same-cell OR and the RELN×Braak interaction, with the EC:Exc.4 reelin+ subtype trajectory.

Inputs : data/leng2021_ec_full.h5ad
Outputs: reln_rorb_braak_trajectory.png
Method : statsmodels Logit with donor-clustered SE

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact
import re
import anndata as ad
import scipy.sparse as sp
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

# Load Leng EC data
import anndata as ad
import scipy.sparse as sp

ec = ad.read_h5ad("data/leng2021_ec_full.h5ad")

def counts(adata, ens):
    x = adata[:, ens].X
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()

reln = counts(ec, 'ENSG00000189056')
rorb = counts(ec, 'ENSG00000198963')

obs = ec.obs.copy()
obs['RELN_c'] = reln
obs['RORB_c'] = rorb
obs['RELN_pos'] = (reln >= 1).astype(int)   # int, not bool — smf.logit needs a 1-D numeric endog
obs['RORB_pos'] = (rorb >= 1).astype(int)
obs['logdepth'] = np.log10(obs['nUMI'].astype(float).clip(lower=1))

# Braak staging (BraakStage is a categorical with values '0','2','6')
obs['braak_n'] = obs['BraakStage'].astype(int)

exc = obs[obs['clusterCellType'] == 'Exc'].copy()
exc = exc.dropna(subset=['braak_n'])

def depth_ctrl_OR(sub):
    if sub.RELN_pos.sum() < 5 or sub.RORB_pos.sum() < 5:
        return np.nan, np.nan, np.nan, np.nan
    try:
        m = smf.logit("RORB_pos ~ RELN_pos + logdepth", data=sub).fit(disp=0)
        OR = np.exp(m.params['RELN_pos'])
        lo, hi = np.exp(m.conf_int().loc['RELN_pos'])
        p = m.pvalues['RELN_pos']
        return OR, lo, hi, p
    except Exception:
        return np.nan, np.nan, np.nan, np.nan

def stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

braak_stages = [0, 2, 6]
braak_labels = ["0", "II", "VI"]

per_stage = []
for b in braak_stages:
    sub = exc[exc.braak_n == b]
    OR, lo, hi, p = depth_ctrl_OR(sub)
    # raw Fisher
    tab = pd.crosstab(sub.RELN_pos.astype(int), sub.RORB_pos.astype(int))
    if tab.shape == (2, 2):
        raw_OR, raw_p = fisher_exact(tab)
    else:
        raw_OR, raw_p = np.nan, np.nan
    per_stage.append({
        'braak': b, 'label': f"Braak {'0' if b==0 else 'II' if b==2 else 'VI'}",
        'n': len(sub),
        'OR': OR, 'lo': lo, 'hi': hi, 'p': p,
        'raw_OR': raw_OR, 'raw_p': raw_p,
        'reln_pct': 100 * sub.RELN_pos.mean(),
        'rorb_pct': 100 * sub.RORB_pos.mean(),
    })

fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

# Panel A: depth-controlled OR by Braak
axA = axes[0]
x = np.arange(len(braak_stages))
ors = [r['OR'] for r in per_stage]
los = [r['lo'] for r in per_stage]
his = [r['hi'] for r in per_stage]
ps = [r['p'] for r in per_stage]
cols = ['#4C72B0', '#DD8452', '#C44E52']

axA.axhline(1.0, color='#999', ls='--', lw=1, zorder=0)
for xi, (o, lo, hi, p, c) in enumerate(zip(ors, los, his, ps, cols)):
    if not np.isnan(o):
        axA.errorbar(xi, o, yerr=[[o-lo], [hi-o]], fmt='o', color=c, ms=8, capsize=4, lw=2)
        s = stars(p)
        if s and s != 'ns':
            axA.annotate(s, (xi, hi), xytext=(0, 8), textcoords='offset points',
                        ha='center', fontsize=15, fontweight='bold', color='black')

axA.set_xticks(x)
axA.set_xticklabels(braak_labels)
axA.set_xlabel('Braak stage')
axA.set_ylabel('Depth-controlled same-cell OR\n(RELN∩RORB)')
axA.set_title('A  Co-occupancy rises with Braak stage\n     (depth-controlled, Leng EC cohort)', loc='left', fontsize=8)
axA.set_yscale('log')
axA.set_ylim(0.3, 3.5)

# Panel B: raw vs depth-controlled ORs
axB = axes[1]
raw_ors = [r['raw_OR'] for r in per_stage]
axB.plot(braak_labels, raw_ors, 'o--', color='#888', ms=7, label='Raw OR', lw=1.5)
axB.plot(braak_labels, ors, 's-', color='#6a1b9a', ms=7, label='Depth-controlled OR', lw=1.8)
axB.axhline(1.0, color='#999', ls=':', lw=1)
axB.set_xlabel('Braak stage')
axB.set_ylabel('Same-cell odds ratio')
axB.set_title('B  Raw vs depth-controlled OR\n     (depth control attenuates but preserves direction)', loc='left', fontsize=8)
axB.legend(fontsize=7)
axB.set_ylim(0.5, 2.8)

# Panel C: detection rates
axC = axes[2]
reln_pcts = [r['reln_pct'] for r in per_stage]
rorb_pcts = [r['rorb_pct'] for r in per_stage]
w = 0.35
xi = np.arange(len(braak_stages))
axC.bar(xi - w/2, reln_pcts, w, label='RELN⁺', color='#C44E52', alpha=0.85)
axC.bar(xi + w/2, rorb_pcts, w, label='RORB⁺', color='#55A868', alpha=0.85)
axC.set_xticks(xi)
axC.set_xticklabels(braak_labels)
axC.set_xlabel('Braak stage')
axC.set_ylabel('% excitatory neurons detected')
axC.set_title('C  Detection rises with Braak stage\n     (motivates depth control)', loc='left', fontsize=8)
axC.legend(fontsize=7)

for ax, L in zip(axes, ['A', 'B', 'C']):
    # normalize_panels style: strip letter from title, add bold outside
    t = ax.get_title(loc='left')
    m = re.match(r'^\s*([A-Za-z])[\).:]?\s{1,4}(\S.*)', t, re.S)
    if m:
        ax.set_title(m.group(2), loc='left')
    panel_letter(ax, L, dx=-0.12, dy=1.03, case='upper', fontsize=15)

fig.suptitle('RELN/RORB same-cell co-occupancy by Braak stage (Leng 2021 entorhinal cortex; depth-controlled)',
             fontsize=10, y=1.02)
fig.tight_layout()
fig.savefig('reln_rorb_braak_trajectory.png', dpi=300, bbox_inches='tight')
