"""
APOE-ε4 dose-response at matched pathology
==========================================

Same-cell OR by APOE isoform group at matched Braak V–VI (SEA-AD MEC).

Inputs : values inlined from data/apoe_isoform_results.json
Outputs: apoe_isoform_convergence.png
Method : statsmodels Logit; matplotlib

Part of the RELN-RORB EC convergence analysis suite. See README.md for the
data-acquisition steps that populate data/, and METHODS_PROVENANCE.md for the
full library-vs-custom-code breakdown. Figure styling is applied via the shared
figure_style module (a thin matplotlib rcParams helper; not required to
reproduce the statistics).
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from figure_style import apply_figure_style, panel_letter
apply_figure_style()

def stars(p):
    return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"

groups=["ε2\n(protective)","ε3/ε3\n(reference)","ε4\n(risk)"]
ors=[0.174,1.004,1.749]; los=[0.08,0.64,1.11]; his=[0.37,1.58,2.76]; ps=[0.0004,0.985,0.016]
x=[0,1,2]; cols=["#4C72B0","#999","#C44E52"]

fig,ax=plt.subplots(figsize=(7,4.6))
ax.axhline(1,color="#999",ls="--",lw=1,zorder=0)
for xi,o,lo,hi,c,p in zip(x,ors,los,his,cols,ps):
    ax.errorbar(xi,o,yerr=[[o-lo],[hi-o]],fmt="D",color=c,ms=9,capsize=5,lw=2)
    s=stars(p)
    if s=="ns":
        ax.annotate(s,(xi,hi*1.06),ha="center",fontsize=9,fontweight="normal",color="#555")
    else:
        ax.annotate(s,(xi,hi*1.06),ha="center",fontsize=14,fontweight="bold",color="black")
ax.plot(x,ors,color="#666",lw=1,ls=":",zorder=0)
ax.set_xticks(x); ax.set_xticklabels(groups)
ax.set_ylabel("Same-cell RELN∩RORB OR at Braak V–VI\n(pathology-matched, depth-controlled)")
ax.set_yscale("log")
ax.set_ylim(0.06,3.8)
ax.set_xlim(-0.4,2.4)
ax.set_title("APOE isoform tracks RELN/RORB convergence along the AD-risk axis\n(ε2→ε3→ε4 trend OR 2.15/step, p=0.006)",fontsize=8.8)
fig.tight_layout()
fig.savefig("apoe_isoform_convergence_v2.png",dpi=300,bbox_inches="tight")
