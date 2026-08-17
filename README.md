# RELN–RORB convergence in the Alzheimer's entorhinal cortex

Analysis code for the study of whether **Reelin (RELN)** and **RORB** co-occupy
the *same* entorhinal-cortex (EC) excitatory neurons in Alzheimer's disease (AD),
and whether that co-occupancy tracks tau pathology.

**Headline finding.** RELN and RORB are expressed independently in EC excitatory
neurons in cognitively normal brain, but converge into the *same* neurons with
advancing Braak tau stage. The convergence is specific to EC excitatory neurons,
tau-graded, module-level (the whole reelin-signalling program co-recruits with
the RORB⁺ EC:Exc.4 identity), amplified by APOE-ε4 dose, and absent from
non-entorhinal cortex in AD and two other primary tauopathies.

All analyses are **secondary analyses of public data**. This is a computational,
hypothesis-generating study.

## What's here

Each script in `analysis/` is a self-contained analysis that reads from `data/`
and writes a figure + result table to the working directory. Run them
independently; there is no master pipeline. Scripts are numbered in narrative
order, not strict dependency order (02 and 05 consume `cellclass_all_regions.csv`
and 14 consumes `staged_ec_all_regions.json`, both produced by 04).

| Script | What it does |
|--------|--------------|
| `01_coexpression_same_cell.py` | Discovery: same-cell RELN×RORB co-occupancy, AD vs normal (Leng EC) |
| `02_region_class_specificity.py` | Region × cell-class specificity |
| `03_braak_trajectory.py` | Braak-graded trajectory + EC:Exc.4 subtype |
| `04_staged_ec_replication.py` | Independent replication across 4 SEA-AD regions → produces derived tables |
| `05_celltype_specificity.py` | Cell-type specificity in the two entorhinal regions |
| `06_apoe_isoform.py` | APOE-ε4 dose-response at matched pathology |
| `07_mec_spatial_colocalization.py` | In-situ spatial co-localization (MERFISH), custom niche test |
| `08_tangle_rorb_tau.py` | RORB in tangle-bearing neurons (Otero-Garcia) |
| `09_grn_rorb_reelin.py` | Gene-regulatory-network importance of RORB (negative result) |
| `10_adtbi_bulk_pathology.py` | Bulk RNA-seq RELN/RORB vs pathology (Allen ADTBI) |
| `11_cross_tauopathy.py` | Cross-tauopathy specificity (Rexach 2024: AD/Pick/PSP) |
| `12_marker_pair_specificity.py` | Module-level specificity control |
| `13_reelin_receptor_arm.py` | Reelin receiver machinery in RORB⁺/double-positive neurons |
| `14_cross_cohort_meta.py` | Cross-cohort random-effects meta-analysis |
| `15_squidpy_niche_crosscheck.py` | Spatial niche re-run through squidpy (standard-tool validation of 07) |
| `16_robustness_controls.py` | Artifact-removal controls: depth thinning, ambient threshold, LODO, donor permutation |
| `17_clustering_circularity.py` | Re-runs the subtype clustering with RELN/RORB removed from the feature set |
| `18_rodriguez_reciprocal.py` | Reciprocal same-cell test in the independent onset cohort (GSE287652) |
| `19_within_donor_regional.py` | Within-donor MEC-vs-hippocampus spatial control (Figure 9F) |
| `20_second_rorb_cluster.py` | RORB-alone negative control: the second RORB-high cluster (EC:Exc.5) is RELN-negative and DP-depleted; depth-controlled and replicated in SEA-AD MEC | §3.2, Supp. Fig. S10 |
| `21_oligodendrocyte_specificity.py` | Oligodendrocyte RELN as a specificity control: class-specific among glia but **brain-wide**, RELN-only, tracks general burden not tau | §3.6 |
| `22_directionality_replication.py` | Why no directional claim is made: the acquisition asymmetry reverses between cohorts while the convergence replicates | §3.12, §4.1 |
| `23_census_dataset_provenance.py` | Enumerates the constituent CELLxGENE Census datasets behind the cross-disorder comparison with their collection DOIs; re-derives the allocortex coverage claim | §3.11, Supp. Table S1 |
| `24_merfish_detection_controls.py` | Measures MERFISH detection error from the object's own control channels; shows decoding noise attenuates rather than creates the co-localization, and tests segmentation spillover | §3.8 |
| `_data_paths.py` | Shared input-path resolution and preflight checks (`require()`); honours `RELN_RORB_DATA` |
| `figure_style.py` | Minimal matplotlib styling helpers (appearance only) |

## Method provenance

Every analysis is **custom code built on established statistical libraries**
(statsmodels, scipy, scikit-learn, squidpy). There is no off-the-shelf pipeline
for the core same-cell question. See **[METHODS_PROVENANCE.md](METHODS_PROVENANCE.md)**
for a per-analysis breakdown of which library does the statistics, what custom
logic sits on top, and whether a named alternative tool exists.

The one custom spatial test (script 07) is independently corroborated by
script 15 using **squidpy**'s peer-reviewed neighborhood-enrichment test.

## Setup

```bash
conda env create -f environment.yml
conda activate reln-rorb
# or: pip install -r requirements.txt
```

Most scripts need only `statsmodels`, `scipy`, `scikit-learn`, `anndata`,
`pandas`, `matplotlib`. Script 15 additionally needs `squidpy`.

## Data

The `data/` directory is **not included** (large public files + one
controlled-derived extract). See **[DATA.md](DATA.md)** for exactly how to
obtain or regenerate each input, with accessions and DOIs.

## Reproducing a result

```bash
cd analysis
python 01_coexpression_same_cell.py     # writes reln_rorb_coexpression_AD.png + stats CSV
```

Each script prints its key statistics to stdout and writes its figure/table to
the current directory.

## Statistical approach (shared across scripts)

- **Positivity is defined per cohort**, so that marker prevalence is comparable across
  datasets whose sequencing depth differs ~14-fold. This matters for reading the code:
    - **Leng EC** (median 2,735 UMIs; scripts 01, 03, 09, 12, 13, 14) — the h5ad `X` holds
      **raw integer counts**, and a gene is positive when its count ≥ 1
      (RELN⁺ 20.3%, RORB⁺ 20.5% of EC excitatory nuclei).
    - **SEA-AD** (median 17,891 UMIs in MEC; script 04 and the extracts it reads) — the
      `reln`/`rorb` columns of `seaad_*_extract.parquet` hold **depth-normalised,
      truncated expression values, NOT raw counts**. The `> 0` test in script 04 is
      therefore a normalised-expression threshold corresponding to roughly **≥7 raw
      transcripts** at that cohort's median depth (RELN⁺ 13.3%, RORB⁺ 20.2%).
      Applying Leng's raw-count rule to SEA-AD would mark a far larger fraction positive
      and would not preserve prevalence matching between cohorts. See DATA.md.
    - **MERFISH** (script 07) — raw per-cell transcript counts from the 433-gene panel.
  Effect sizes depend on threshold stringency; the direction of the interaction does not
  (threshold sweep in script 16).
- Same-cell co-occupancy is tested by **Fisher's exact test** and by
  **depth-controlled logistic regression**: `RORB⁺ ~ RELN⁺ × disease(or Braak)
  + log₁₀(library size)`, with the interaction term as the estimand.
- Sequencing depth is a required covariate (detection of both genes scales with
  library size, and depth co-varies with disease stage).
- Per-stage and cross-region estimates use **donor-clustered** robust SEs.
- Multiple-testing correction is Benjamini–Hochberg within each test battery.

## Changes since the first public release (v1 → v3)

**Read this if you cloned v1.** Everything below is cumulative — v2.0 and v2.1 were
prepared but never published, so this is the full diff against the only release that
has been on GitHub. File-level summary: **11 files added, 4 modified, 0 removed**;
the 19 other files are byte-identical to v1.

### 1. A documentation error corrected — read this before modifying any script

v1's *Statistical approach* section stated that "a gene is scored positive at a raw
count ≥ 1" as though one rule applied to every cohort. **That was wrong, and it
mattered.** The two cohorts store marker values differently:

- **Leng EC** (`leng2021_ec_full.h5ad`) holds genuine **raw counts** → positivity is `> 0`.
- **SEA-AD extracts** (`seaad_*_extract.parquet`) hold **depth-normalised truncated**
  values that *resemble* small counts → positivity is `>= 1` on that column, which at
  SEA-AD depth corresponds to roughly ≥ 7 raw transcripts.

Applying the Leng rule to SEA-AD scores far more nuclei positive and destroys the
prevalence matching between cohorts. The code in v1 was already correct — only the
prose was wrong — so **no v1 result changes**. The rule is now documented in three
places deliberately (this README, a warning block in `DATA.md`, and an inline comment
at the thresholding lines), each with the diagnostic a skeptic can run: the rule
reproduces RELN⁺ 13.3% / RORB⁺ 20.2% of MEC excitatory nuclei. Scripts 20–22 assert
this on load and fail loudly rather than proceed on a mis-calibrated threshold.

### 2. Donor-clustered inference is now the reporting standard

v1 reported some cell-level p-values. With 10⁵–10⁶ nuclei from ~10²  donors these are
anti-conservative by orders of magnitude — the effective sample size is the donor
count, not the cell count. Every script now reports donor-clustered robust standard
errors (`cov_type="cluster"`), and `04_staged_ec_replication.py` was modified to make
the model specification explicit: **RORB-positivity is the outcome**, RELN-positivity
and Braak stage the predictors. The reverse specification gives a different
interaction OR (1.24 vs 1.48), so the direction is now stated rather than implied.

### 3. Seven new analyses (scripts 16–22)

Four were added in response to a review round, three in later work:

| Script | What it tests | Why it exists |
|---|---|---|
| `16_robustness_controls.py` | Binomial depth thinning to common depth, ambient-RNA floor with a ≥2 threshold, leave-one-donor-out, donor-label permutation | Removes rival explanations rather than adjusting for them — so its p-values are *larger* than the primary estimate's by construction |
| `17_clustering_circularity.py` | Re-runs subtype clustering with RELN and RORB removed from the feature set entirely | The double-positive cluster was defined on features that can include the two genes being counted |
| `18_rodriguez_reciprocal.py` | Reciprocal same-cell test in an independent onset cohort (GSE287652) | Note: the published object is integrated with this study's discovery cohort, so it **must be subset** or the replication is partly circular |
| `19_within_donor_regional.py` | Within-donor MEC-vs-hippocampus spatial control | Removes between-donor variation from the regional claim |
| `20_second_rorb_cluster.py` | The second RORB-high cluster (EC:Exc.5) is RELN-negative and double-positive-*depleted* | Internal negative control: rules out "RORB is just a proxy for a RELN-adjacent state" |
| `21_oligodendrocyte_specificity.py` | RELN⁺ oligodendrocytes rise with Braak — but brain-wide, RELN-only, and tracking general burden rather than tau | Specificity control. The rise is present in primary visual cortex, where the neuronal convergence cannot be defined; the *pairing* is entorhinal-only |
| `22_directionality_replication.py` | The conditional acquisition asymmetry **reverses** between cohorts while the convergence replicates | Documents why the manuscript withdraws its directional claim |

### 4. Two claims withdrawn

Both were in the manuscript, neither was ever in v1's code, but the scripts now
document why they are gone:

- **The directional refinement.** v1-era work read the acquisition asymmetry as "a
  stable RORB⁺ identity acquiring RELN." Two problems: acquisition and differential
  survival of RELN⁺RORB⁻ neurons predict the same cross-sectional pattern and cannot
  be separated in autopsy tissue; and the direction reverses in the larger cohort
  (script 22). Only the convergence itself replicates.
- **Pooled cross-cohort magnitude.** Between-cohort heterogeneity is too high to
  report a single pooled effect; `14_cross_cohort_meta.py` reports per-cohort
  estimates with heterogeneity rather than a pooled diamond.

### 5. Documentation added

- `CHANGELOG.md` — new file, per-release detail.
- `DATA.md` — adds the GEO accession for the onset cohort (GSE287652) and the
  positivity warning block; documents which marker columns carry which semantics.
- `METHODS_PROVENANCE.md` — entries for scripts 16–22, including the design decisions
  that are not obvious from the code (why robustness controls cost power, why the
  onset cohort must be subset, why subclasses with fewer than ten double-positive
  cells are reported as unestimable rather than quoted).

### Verification state of this release

Every script that has obtainable inputs has been **executed cold from a clean
`data/` directory** — as a subprocess, not an interactive import — and reproduces
the manuscript values:

| Script | Cold-run result |
|---|---|
| `16_robustness_controls.py` | primary OR 1.481 (p = 1.5×10⁻³); LODO range 1.406–1.621; thinning 1.394 |
| `17_clustering_circularity.py` | 5.28-fold with markers in the feature set, 5.06-fold with them removed |
| `19_within_donor_regional.py` | MEC OR 1.991 (p = 1.1×10⁻³⁸) vs HPF 1.235 (p = 3.3×10⁻³) |
| `20_second_rorb_cluster.py` | depth-controlled OR 68.0 [45.9–100.6]; L2/3 IT correctly flagged unestimable (7 DP cells) |
| `21_oligodendrocyte_specificity.py` | oligodendrocyte 1.207 (p = 0.004), all other glia 0.88–0.96; four-region 1.19–1.30 |
| `22_directionality_replication.py` | asymmetry reverses (Leng 1.430/1.076 vs SEA-AD 1.147/1.545) |

`18_rodriguez_reciprocal.py` (renamed from `18_roussarie_reciprocal.py` — see CHANGELOG) could not be run end-to-end: its input derives from a
~1.5 GB Seurat object that is not redistributable and must be downloaded from GEO
(GSE287652) and reduced by the R snippet in its docstring. Its preflight check and
model code are exercised; the reduction step is not.

**Missing-input behaviour is now checked, not tracebacked.** All seven scripts share
`analysis/_data_paths.py`, which resolves inputs through `require()` and exits with the
filename, its approximate size, and the `DATA.md` section documenting its source —
rather than a pandas or h5py traceback that names nothing. Verified: with an empty
`data/`, all seven exit 1 with zero traceback lines and a named provenance pointer.
Set `RELN_RORB_DATA` to read inputs from elsewhere:

```bash
RELN_RORB_DATA=/scratch/reln_rorb/data python analysis/16_robustness_controls.py
```

## Citation & license

See [LICENSE](LICENSE) (MIT for code). Underlying datasets retain their own
licenses and access terms (see DATA.md). If you use this code, please cite the
associated manuscript and the primary datasets (Leng et al. 2021; SEA-AD;
Rexach et al. 2024; Otero-Garcia et al.; Allen ADTBI).
