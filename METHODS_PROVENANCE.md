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

## Smoke-test bug fixes

The scripts were run end-to-end as a packaging smoke-test, which surfaced and
fixed four bugs in the lineage-extracted code before release:

1. **Depth column** — code referenced `obs['total_counts']`, which is absent
   from the Leng h5ad; corrected to `obs['nUMI']` (scripts 01, 03).
2. **Braak mapping** — a string-keyed map (`'Braak 0'`) applied to
   integer-valued categories (`'0'`/`'2'`/`'6'`) yielded all-NaN; replaced with
   a direct integer cast (scripts 01, 03).
3. **Wrong input file** — script 01 (the discovery figure) had been extracted
   pointing at the MERFISH file, which carries no Braak labels; rewritten as a
   clean standalone that loads the Leng EC h5ad.
4. **Boolean endog** — `smf.logit` rejects a boolean outcome column (reads it as
   a 2-column categorical), which a bare `except` silently turned into NaN and
   left a figure panel empty; fixed by casting positivity flags to int
   (script 03).

After these fixes all six run scripts reproduce the published values exactly.

## Library versions (as run)

- Python 3.11
- statsmodels, scipy, scikit-learn, anndata, h5py (analysis)
- cellxgene-census (data access)
- squidpy 1.8.2, scanpy (spatial cross-check)
- matplotlib (figures)

Exact pinned versions are in `environment.yml`.

---

## Addendum — revision-round analyses (scripts 16–19)

Added during peer review. All four re-apply methods already described above; none
introduces new methodology.

### 16 — Artifact-removal robustness controls
`statsmodels` Logit as elsewhere. What is custom is the *design* of the controls:
they **remove** an artifact rather than adjusting for it, and each therefore costs
statistical power by construction. Binomial thinning (`numpy.random.binomial`)
downsamples every nucleus to a common depth so that no depth covariate is needed;
the ≥2-transcript threshold discards any call a single ambient molecule could
produce; the donor-label permutation makes no parametric assumption about the
variance estimator at all. Reporting a *larger* p-value after these controls is the
expected outcome — the gain is the elimination of a rival explanation, not a
smaller p.

### 17 — Clustering circularity control
`scanpy` preprocessing + `sklearn.cluster.KMeans` + `scipy.stats.hypergeom`, the
same stack as script 03. The one methodological point: double-positive labels are
fixed from the **raw counts before feature selection**, so removing RELN/RORB from
the highly-variable-gene set changes only the clustering, never the labels being
counted. Absolute fold-enrichment varies with preprocessing choices (scaling,
HVG selection) across roughly 4.4–5.3×; the *comparison* between the two
configurations is like-for-like within a single run and is what the control tests.

### 18 — Reciprocal test in an independent onset cohort
`statsmodels` Logit, depth-controlled, donor-clustered — identical to the primary
model. Two data-handling decisions matter more than the statistics: the published
object is integrated with our own discovery cohort and must be subset to the
originating study's nuclei, and the original authors' own sample exclusions are
honoured. Power is low (4 donors per group after exclusions); the result is
reported as consistent-with, not proof-of, an absence.

### 19 — Within-donor regional spatial control
`statsmodels` Logit plus a label-shuffle kNN permutation (`sklearn`
`NearestNeighbors`), the same custom niche test as script 07 and cross-checked
against squidpy in script 15. The permutation holds cell positions and the
double-positive count fixed, so the null preserves both tissue geometry and marker
prevalence. Both sections derive from a single donor, confirmed by the data
generators, which converts a between-section comparison into a within-donor
regional control.

### A note on positivity thresholds
The single most important thing to understand before modifying any script: the
positivity rule is **not uniform across cohorts**, by design. Leng and MERFISH
threshold genuine raw counts; the SEA-AD extracts threshold depth-normalised,
truncated values that *look* like small counts. See the README "Statistical
approach" section and the warning block in DATA.md.

### `20_second_rorb_cluster.py`

Library code: scikit-learn (PCA, KMeans), scipy.stats (hypergeom, fisher_exact, spearmanr), statsmodels (logistic regression). Custom code: the conditional framing (restricting to RORB⁺ neurons before comparing RELN rates), the hypergeometric *lower*-tail depletion test, and the depth-quintile stratification.

Two design points are deliberate. First, the depth control is necessary rather than decorative: the two clusters differ ~2× in sequencing depth and Leng positivity is a raw-count threshold, so the conditional gap must be shown to survive depth adjustment before it can be interpreted. Second, subclasses resting on fewer than ten double-positive cells are reported as unestimable rather than quoted — maximum-likelihood separation drives such estimates to zero or infinity with no meaningful confidence interval.


### `21_oligodendrocyte_specificity.py`

statsmodels `Logit` with donor-clustered robust covariance throughout; the custom part is the control structure rather than any estimator. The script is written so that the bounding tests run *before* the headline is stated: class specificity across the six non-neuronal subclasses, then regional specificity across all four SEA-AD regions, then marker specificity, then the pathology axis, then depth stratification.

Three design points are deliberate. First, the class-specificity sweep is the substantive control for ambient RNA — a contamination floor rising with disease would lift every non-neuronal class, and it lifts one, so the comparison does work that no single-class model can. Second, the regional sweep reports the *contrasting* neuronal quantity in the same table, because the point being made is a dissociation and a reader needs both halves; where a region has no RELN substrate the convergence is labelled undefined rather than null, since a missing denominator is not a negative result. Third, depth is reported two ways — within quintiles and as the depth-versus-stage correlation — because the sign of that correlation is what rules the artifact out: depth *declines* with Braak in these cells, so it works against the effect rather than producing it.

Validated by module import against known-good data, which caught a substring test that scored the level `"No dementia"` as dementia. The script now asserts its expected level set.

### `22_directionality_replication.py`

statsmodels `Logit`, donor-clustered, on subsets conditioned by marker status. No custom estimator; the contribution is the tabulation.

The script exists to make a *negative* result auditable, so it deliberately reports quantities that are not claims: both conditional acquisition slopes, the convergence interaction, and both markers' marginal trends, in each cohort. The marginals are included because they explain the reversal — the two conditional slopes differ only by the difference in the marginals, and those differ between cohorts, so the reversal is arithmetic rather than contradictory. The docstring states the identifiability problem (acquisition versus differential survival predict the same cross-sectional asymmetry) that no test in the script can resolve, so the limitation travels with the code rather than living only in the manuscript.


### `23_census_dataset_provenance.py`

No custom estimator; the contribution is provenance recovery. The cross-disorder row of Table 1
previously read "many (public atlases)" and no dataset identifiers were recorded anywhere in the
codebase, which made the row's quality-control claim uncitable — Census standardises schema and
ontology but does not re-run QC, so provenance devolves to each constituent study.

Three design points. First, `is_primary_data == True` is load-bearing rather than hygienic:
Census carries the same nuclei under multiple collections, and without that filter every count is
inflated by duplication. Second, the script re-derives the allocortex coverage claim from the
build rather than asserting it, and warns if coverage differs — that claim is a property of a
versioned resource and can legitimately change as datasets are added. Third, the Census build is
pinned rather than floated to "latest", for the same reason.

The proxy handling in `open_census` is not incidental: TileDB's S3 client uses raw sockets and
fails DNS resolution behind a proxy rather than honouring HTTP(S)_PROXY the way most Python HTTP
clients do, so the proxy must be passed into the TileDB config explicitly.


### `24_merfish_detection_controls.py`

The 433-gene MERFISH release has no peer-reviewed methods paper, so its detection error cannot be
cited. It can, however, be measured: the delivered object retains the control channels. The
script's structure follows from a distinction the control channels alone do not make.

Decoding error — a codeword misread as another gene — is what the control probe, control codeword
and unassigned codeword channels bound. Because it is independent per codeword, it *attenuates* an
association rather than creating one, so the script inverts the false-positive rate and reports the
corrected odds ratio; that the correction moves the estimate away from unity is the substantive
result, not the raw error rate.

Segmentation spillover is not captured by any control channel and is the more consequential failure
mode for a same-cell claim in intact tissue. It is tested two ways: stratifying by segmentation
method (the sections carry three, and the membrane-based boundary stain is least prone to spillover),
and adding segmented cell area to the depth-controlled model. Double-positive cells *are* larger,
which is the signature spillover would produce, so the null area term after depth adjustment is the
test that discriminates the two explanations.

The control-panel size is not recorded in the object, so the per-gene rate is bounded across
plausible panel sizes rather than asserted, and every conclusion is drawn from the conservative end.
