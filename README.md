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
| `25_multiregion_extract.py` | Streams the SEA-AD 2026 ten-region release (161 GB), extracts RELN/RORB + per-nucleus QC to small parquets; refuses the 2022 Reference-MTG object by filename | §2, DATA.md |
| `26_regional_subclass_analysis.py` | Ten-region convergence in two donor cohorts, then **within subclass**: the convergence is entorhinal-IT, MTG's region-level signal is in the L4 IT control. Includes the same-donor cross-region correlation | §3.6 |
| `27_cps_pathology_axes.py` | Continuous tau vs amyloid pseudo-progression axes, mutually adjusted; tau survives brain-wide, amyloid does not | §3.5 |
| `28_qc_filter_sensitivity.py` | Mitochondrial and doublet filter curves per region and cohort; doublet removal does not weaken the co-occupancy | §3.5 |
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

## Changes since the first public release (v1 → v4)

**Read this if you cloned v1 or v3.0.** Diff against the published v3.0 tag: **4 files added, 6 modified, none removed**; the other 28 are byte-identical.

**Read this if you cloned v1 or v3.0.** v2.0 and v2.1 were never published. v3.0 was
published and tagged. This release adds the ten-region reanalysis **and** a threshold
patch that was prepared after v3.0 but never uploaded, so v4.0 carries both.

### 1. Braak threshold made explicit (prepared after v3.0, first shipped here)

`01_coexpression_same_cell.py` dichotomised AD status as `braak > 0`. A code reviewer
correctly flagged that this would count Braak I as AD, which is wrong —
transentorhinal-only tau is common in cognitively normal ageing.

**No published result changes.** Verified across every cohort: no single-cell dataset
here contains a Braak I donor (Leng 0/2/6; SEA-AD 0, II–VI across all ten regions and
84 donors), so `braak > 0` and `braak >= 2` select identical nuclei — 0 of 10,780 in
Leng change label, and every reported value is byte-identical. The only Braak I
material is in the ADTBI bulk cohort (40 of 377 samples), which is analysed with
continuous correlations and applies no dichotomy.

The threshold is now written `braak >= 2` with the rationale in a comment, plus an
assertion that fails loudly if Braak I ever appears in the input. The guarantee comes
from the code rather than from which cohorts happened to be used.

### 2. Ten-region reanalysis (SEA-AD 2026-06-22 release)

Four new scripts, 25–28. The scientific result that changed the manuscript: a
region-level model flags two regions (MEC and MTG), but refitting **within subclass**
separates them completely. MEC's double-positives are 97.7% entorhinal-IT, where
24.19% of RORB⁺ nuclei are RELN⁺; MTG's are 84.6% L4 IT, where the figure is 0.44%
and the within-subclass interaction is not significant. L4 IT is this study's own
RORB-high/RELN-negative negative control. MTG is therefore a **contrast**, not a
replication — a detection-floor drift in a very large RORB⁺ population, significant
only because its denominator is 148,000 nuclei.

Supporting results: the same-donor cross-region correlation is null (rho = −0.16,
n = 74, CI excludes rho > 0.07), so MEC and MTG are not one donor-level process;
continuous pathology axes separate tau from amyloid brain-wide where the ordinal
variables could not; and the interaction survives every mitochondrial and doublet
threshold tested.

### 3. Guards added because these mistakes were actually made

- **Filename guard (script 25).** `MTG/RNAseq/` holds two `*_final-nuclei.h5ad`
  objects; the 2022 neurotypical **Reference** was downloaded before this was caught.
  Both `SEAAD_` and `2026-06-22` are now required in every key.
- **`MIN_DP = 10` (scripts 26, 28).** Hippocampus has 6 double-positive nuclei in
  119,096 and returns OR 8.3 with a "significant" negative interaction if left
  unguarded. Such cells report `UNESTIMABLE`, never an odds ratio.
- **Braak map assertion (script 26).** An unmapped stage — Braak I above all — raises
  rather than collapsing into a neighbouring stage.
- **Donor-level standardisation (script 27).** CPS is a donor property, so "per SD"
  uses donor moments; nucleus-weighting lets high-yield donors dominate the SD and
  shifts the per-SD OR. Documented inline with both values.

### 4. Documentation

`DATA.md` gains the 2026 release section: file layout, the 161 GB bandwidth warning,
the two file-selection traps, the finding that the `all-nuclei` objects do **not**
contain QC-failed expression, and an explicit warning not to mix the normalised and
raw-UMI positivity rules. The per-cohort Braak table now states that stage I is
absent across all ten regions, verified per region.

### Verification state

Scripts 26–28 were run cold as subprocesses from an empty `data/` directory (all
three exit with named provenance, no tracebacks) and again with data present, where
they reproduce the manuscript values exactly. Script 25's filename guard and live
region resolution were tested directly; its 161 GB download path was exercised once
during the original analysis but is not re-run on every change.

## Citation & license

See [LICENSE](LICENSE) (MIT for code). Underlying datasets retain their own
licenses and access terms (see DATA.md). If you use this code, please cite the
associated manuscript and the primary datasets (Leng et al. 2021; SEA-AD;
Rexach et al. 2024; Otero-Garcia et al.; Allen ADTBI).
