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

- A gene is scored **positive** if its raw count ≥ 1 in a cell.
- Same-cell co-occupancy is tested by **Fisher's exact test** and by
  **depth-controlled logistic regression**: `RORB⁺ ~ RELN⁺ × disease(or Braak)
  + log₁₀(library size)`, with the interaction term as the estimand.
- Sequencing depth is a required covariate (detection of both genes scales with
  library size, and depth co-varies with disease stage).
- Per-stage and cross-region estimates use **donor-clustered** robust SEs.
- Multiple-testing correction is Benjamini–Hochberg within each test battery.

## Citation & license

See [LICENSE](LICENSE) (MIT for code). Underlying datasets retain their own
licenses and access terms (see DATA.md). If you use this code, please cite the
associated manuscript and the primary datasets (Leng et al. 2021; SEA-AD;
Rexach et al. 2024; Otero-Garcia et al.; Allen ADTBI).
