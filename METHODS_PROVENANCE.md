# Methods Provenance — libraries vs. custom code

Every analysis in this repository was implemented as **custom code written for
this project**, built on top of **established, peer-reviewed statistical
libraries**. No off-the-shelf end-to-end pipeline was used for the core
same-cell co-occupancy question (there is no standard named pipeline for "do
two markers occupy the same cell more with disease stage"). The statistics
themselves come entirely from standard libraries; the code that wires them to
this hypothesis is bespoke.

This table is the transparency record. For each analysis: the statistical
engine (library, citable), the custom logic layered on top, and whether a
named alternative tool exists.

| # | Analysis | Statistical engine (library) | Custom logic (written here) | Named alternative tool |
|---|----------|------------------------------|------------------------------|------------------------|
| 01 | Same-cell co-occupancy (Leng) | `statsmodels` Logit; `scipy.stats` fisher_exact | positivity call (count≥1); depth-controlled interaction model | none standard |
| 02 | Region×class specificity | (consumes 04 output); `matplotlib` | grouped-OR presentation | none |
| 03 | Braak trajectory | `statsmodels` Logit, donor-clustered SE | per-stage stratification; subtype trajectory | none standard |
| 04 | Staged SEA-AD replication | `statsmodels` Logit, cluster-robust SE | per-region×stage loop; EC-IT subclass; APOE panel | none standard |
| 05 | Cell-type specificity (EC) | (consumes 04 output); `matplotlib` | class-restricted OR presentation | none |
| 06 | APOE-ε4 dose | `statsmodels` Logit | isoform-group model at matched pathology | PLINK-family (genetics) — not applicable here |
| 07 | Spatial co-localization | `scipy.stats` fisher_exact; `statsmodels` Logit | **custom kNN neighborhood + laminar-axis test** | **squidpy / Giotto / CellCharter** |
| 08 | Tangle-bearing RORB | `statsmodels` Logit, donor FE | within-donor AT8± contrast | none standard |
| 09 | GRN importance | `scikit-learn` GradientBoostingRegressor | random-gene permutation null | SCENIC / GRNBoost2 / arboreto |
| 10 | Bulk pathology corr | `scipy.stats` spearmanr | region-wise partial correlations | none |
| 11 | Cross-tauopathy (Rexach) | `statsmodels` Logit; `scipy.stats` fisher | per-disorder×region loop | none standard |
| 12 | Marker-pair specificity | `statsmodels` Logit | module vs. frequency-matched-control design | none |
| 13 | Receptor arm | `statsmodels` Logit | receptor-group comparison + AD shift | CellPhoneDB (ligand-receptor) — different question |
| 14 | Cross-cohort meta | `scipy` (custom DerSimonian–Laird) | random-effects pooling of interaction ORs | `metafor` (R), `PythonMeta` |
| 15 | **Spatial niche cross-check** | **`squidpy` nhood_enrichment** | none — standard tool, as validation | this IS the named tool |

## The one place a named tool now backs a custom test

The spatial niche claim (script 07) originally used a **custom** kNN
neighborhood-enrichment test. Script **15** re-runs the same question through
**squidpy's** peer-reviewed `nhood_enrichment` permutation test
(Palla et al. 2022, *Nat Methods*). The two agree: double-positive neurons
self-cluster (squidpy z ≈ 42; custom z ≈ 25) and associate positively with
RORB⁺ neighbors while being excluded from marker-negative regions. The custom
result is therefore corroborated by an independent, standard implementation.

## What this means for interpretation

- The **methods are standard and citable** (statsmodels, scipy, scikit-learn,
  squidpy); the **pipelines are bespoke** for this hypothesis.
- This is normal for a computational hypothesis-generating study, but it puts
  weight on **code availability** — which is why this repository exists. Each
  script is self-contained, parameters are in its docstring, and results
  replicate across independent cohorts run through the same code.
- Scripts **01, 03, 10, 12, 13, 14** have been run end-to-end from their
  input files and reproduce the published values exactly (e.g. same-cell
  interaction OR = 2.37, p = 9.5×10⁻⁶; cross-cohort pooled OR = 1.26).
  Scripts 01 and 09 are clean standalone re-implementations from documented
  parameters (the lineage-extracted originals were patch-chains); all others
  are the exact analysis code with data paths parameterised.
- The remaining scripts are static-lint-clean (every used module imported)
  and compile, but have not each been run against their own large inputs.
  A fully independent re-implementation of every analysis remains a
  reasonable reviewer request.

## Library versions (as run)

- Python 3.11
- statsmodels, scipy, scikit-learn, anndata, h5py (analysis)
- cellxgene-census (data access)
- squidpy 1.8.2, scanpy (spatial cross-check)
- matplotlib (figures)

Exact pinned versions are in `environment.yml`.
