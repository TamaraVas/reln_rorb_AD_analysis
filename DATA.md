# Data acquisition

The `data/` directory is not committed (file sizes + access terms). Below is how
to obtain or regenerate every input each analysis script expects. Filenames match
the paths hard-coded (as relative `data/...`) in the scripts.

## Public single-cell / spatial datasets

### `data/leng2021_ec_full.h5ad`  — Leng et al. 2021 entorhinal cortex
- **Source:** CELLxGENE Census, dataset `2727d83a-0af0-443a-bff8-58dc7028289a`
- **Paper:** Leng et al. 2021, *Nat Neurosci*, DOI 10.1038/s41593-020-00764-7 (PMID 33432193)
- **Fetch (Python):**
  ```python
  import cellxgene_census
  with cellxgene_census.open_soma(census_version="2025-01-30") as c:
      adata = cellxgene_census.get_anndata(
          c, organism="Homo sapiens",
          obs_value_filter="dataset_id=='2727d83a-0af0-443a-bff8-58dc7028289a'")
  adata.write_h5ad("data/leng2021_ec_full.h5ad")
  ```
- Used by scripts 01, 03, 09, 12, 13, 14.

### `data/seaad_merfish_mec.h5ad` — SEA-AD MERFISH, medial entorhinal
- **Source:** Allen Institute SEA-AD; Allen Brain Cell (ABC) Atlas MERFISH release.
  Specimen `1444201261` (MEC), 433-gene panel, with scANVI annotations + spatial coords.
- **Access:** https://portal.brain-map.org/explore/seattle-alzheimers-disease
  (ABC Atlas `abc_atlas_access` tooling). This file carries no donor/Braak label.
- Used by scripts 07, 15.

### SEA-AD multi-region snRNA-seq per-nucleus extracts
`data/seaad_{mec,lec,hip,v1c}_extract.parquet` — one row per nucleus with columns
`reln`, `rorb`, `umis`, Braak, CERAD, APOE genotype, cell class/subclass.

> **IMPORTANT — column semantics.** The `reln` and `rorb` columns hold
> **depth-normalised, truncated expression values, not raw UMI counts.** They are
> small integers (0–5) because the normalised value is truncated, which makes them
> *look* like counts. Thresholding at `> 0` (as script 04 does) is a
> normalised-expression threshold equivalent to roughly **≥7 raw transcripts** at this
> cohort's median depth of ~17,891 UMIs, and yields RELN⁺ 13.3% / RORB⁺ 20.2% of MEC
> excitatory nuclei. Diagnostic if you regenerate these files: under a genuine
> raw-count rule, detection rises steeply with library size; in these columns RELN
> detection is flat across depth deciles and RORB detection *falls* (30.8% → 14.2%).
> The Leng h5ad (`leng2021_ec_full.h5ad`) is different — that one *is* raw counts.
- **Source:** SEA-AD multi-region 10x snRNA-seq (2026 release), streamed per-donor
  from the public bucket `s3://sea-ad-single-cell-profiling/` and reduced to the
  RELN/RORB/depth/pathology columns. The reduction script is `04`'s data-prep
  (documented in its docstring); the raw per-donor h5ad objects are the ABC Atlas
  SEA-AD-Multiregion-10X release.
- Regions/donors: MEC 81, LEC 34, HIP 40, V1C 43.
- Used by script 04 (which then writes the two derived files below).
- **Per-region files.** Scripts 21 and 22 read more than MEC, so the four regional
  reductions are named explicitly: `data/seaad_mec_extract.parquet`,
  `data/seaad_lec_extract.parquet`, `data/seaad_hip_extract.parquet`,
  `data/seaad_v1c_extract.parquet`. All four carry the same columns
  (`reln`, `rorb`, `umis`, `donor`, `Braak`, `Thal`, `CERAD`, `ADNC`, `Cognitive`,
  `APOE`, `Subclass`, `Class`) under the same depth-normalised positivity rule
  documented in the warning above. `Cognitive` has exactly two levels, `Dementia`
  and `No dementia` — match it exactly, never with a substring test, since the
  negative level contains the positive one as a substring.

### `data/oterogarcia_tangle_cells.parquet` — Otero-Garcia et al. tangle-sorted somata
- **Source:** Otero-Garcia et al., single-soma RNA-seq of AT8⁺ tangle-bearing vs
  tangle-free prefrontal neurons. Per-cell RELN/RORB detection, AT8 status, donor,
  pathology. Obtain from the paper's GEO deposit.
- Used by script 08.

### `data/adtbi_merged.parquet` — Allen Aging, Dementia & TBI bulk RNA-seq
- **Source:** Allen Institute Aging, Dementia & TBI study (https://aging.brain-map.org).
  377 RNA-seq samples, 107 donors, 4 regions (HIP/TCx/PCx/FWM), merged with per-donor
  Braak and CERAD. FPKM table + donor metadata joined on donor_id.
- Used by script 10.

### Rexach et al. 2024 cross-dementia atlas (streamed, not a local file)
- **Source:** CELLxGENE Census dataset `ac0c6561-7a48-4185-af6f-af799f699172`
- **Paper:** Rexach et al. 2024, *Cell*, DOI 10.1016/j.cell.2024.08.019
- Script 11 streams RELN/RORB directly from Census (see its header). To run offline,
  pre-download that dataset's RELN/RORB columns to `data/rexach_reln_rorb.h5ad` and
  point the script at it.

### `data/rod_reln_rorb.parquet` — Rodriguez-Rodriguez et al. 2025 onset cohort
- **Source:** GEO **GSE287652**, supplementary file
  `GSE287652_Exc_neu_integrated_with_Leng.rds.gz` (~551 MB compressed Seurat object).
- **Paper:** Rodríguez-Rodríguez et al. 2025, *bioRxiv*, DOI 10.1101/2024.12.31.630881
  (PMID 39803521). Human entorhinal cortex at the onset of AD neuropathology,
  10 donors: C1–C5 (Braak 0) and AD1–AD5 (Braak II).
- **Two caveats that change the analysis** (both handled in script 18):
  1. The published excitatory object is **integrated with the Leng cohort**, which is
     our own primary discovery cohort. Subset to `set == "Rodriguez"` before any
     test, or the "independent replication" is partly circular.
  2. Samples **C5 and AD1** were excluded by the original authors because amyloid
     did not match their pathology group; the object encodes this as
     `BraakStage == "other"`. Script 18 honours that exclusion.
- **Prepare:** run the `EXTRACT_R` snippet in script 18's header (needs R + Seurat +
  arrow) to write `data/rod_reln_rorb.parquet`.
- Used by script 18.

### `data/seaad_merfish_hpf.h5ad` — SEA-AD MERFISH, hippocampal formation
- **Source:** as the MEC section above; specimen `1444211893` (HPF).
- **Donor provenance (confirmed by the data generators):** this section and the MEC
  section come from the **same donor, H24.30.005**, so the pair forms a *within-donor*
  regional control. Both were profiled to localize transcriptionally defined cell
  types in a **neurotypical** context; the donor is unlikely to carry significant
  pathological burden, and MEC/HIP spatial data across levels of pathological burden
  have not yet been generated. No Braak or disease label attaches to either section.
- Used by scripts 07 (MEC) and 19 (both).

## SEA-AD Multiregion 2026 release (ten regions) — produced by script 25

Public, no access agreement. Streamed from the SEA-AD S3 bucket by
`analysis/25_multiregion_extract.py`, which downloads one object, extracts two gene
columns plus per-nucleus QC metadata, writes a small parquet, and deletes the h5ad.

| File | Size | Contents |
|---|---|---|
| `data/seaad2026_<REGION>_final.parquet` | 4–27 MB each | RELN/RORB (normalised **and** raw UMI) + 17 obs columns, per region |
| `data/seaad2026_MEC_allnuclei.parquet` | ~19 MB | MEC "all-nuclei" object (see the warning below) |
| `data/Global_and_Local_CPS.csv` | ~55 KB | donor-level continuous pseudo-progression scores, tau and amyloid, global and local |

Ten regions: `MEC LEC HIP ITG MTG STG FI AnG PFC V1C`. 6,013,346 nuclei total, of
which 3,110,538 are glutamatergic; 34–84 donors per region. That total equals the
release's own `Used in analysis` count, which is a useful external check that the
extraction lost nothing.

**Bandwidth warning.** The source objects total **161 GB**. Script 25 holds one at a
time (peak 39 GB) and skips regions already extracted, so it is resumable — but a
cold run is a multi-hour download. The extracted parquets total ~135 MB, so extract
once and keep them.

**Two traps, both of which bit us:**

1. **`MTG/RNAseq/` contains TWO `*_final-nuclei.h5ad` objects** — the 2026 SEA-AD
   cohort (32.98 GB) and a 2022 neurotypical **Reference** (6.38 GB). Selecting the
   first match silently picks the Reference, which has no AD staging and would enter
   a regional comparison as a wrong-cohort ringer. Script 25 requires both `SEAAD_`
   and `2026-06-22` in every filename and refuses anything else.
2. **The prefrontal data is under `PFC/`**, and is labelled `DFC` in the per-nucleus
   metadata. The `DFC/` and `DLPFC/` prefixes exist but are **empty placeholders**;
   a plan built from those names silently drops the region.

Related: do **not** enumerate the bucket by paginating all keys to discover which
regions exist. It holds >500,000 objects, so a loop with any page cap under-reports
regions. Use a delimiter-based prefix query (script 25's `lsdir`).

**The `all-nuclei` objects do NOT contain the QC-failed nuclei.** Their cell counts
are identical to `final-nuclei` (MEC 970,118 both; `Used in analysis` 100% True; LEC
182,260 both), so expression for the ~2.6M QC-failed nuclei is **not distributed**.
The multiregion `*_all-nuclei_cell-annotation.csv` (8,641,342 rows, 2,627,996 with
`Used in analysis = False`) records which barcodes failed and under which gate, plus
their taxonomy — but carries no expression. Per-subclass QC failure rates are
therefore computable; a "does the convergence appear in discarded nuclei?" test is not.

**Positivity.** Both normalised and raw-UMI values are stored for the two genes so
either rule can be applied without re-downloading. The manuscript uses the
normalised rule (`>= 1`), which reproduces the published prevalences exactly
(MEC excitatory: RELN⁺ 13.339%, RORB⁺ 20.178%). The raw-UMI rule gives 36.6% / 44.3%
in the same cells — a different threshold, not a robustness check. Do not mix them.

## Braak stages present, per cohort

Verified against the delivered objects, not taken from the source publications:

| Cohort | Braak values present | Braak I present? |
|---|---|---|
| Leng 2021 EC (snRNA-seq) | 0, 2, 6 | no |
| SEA-AD, all ten regions (2026 release) | 0, II, III, IV, V, VI | no — verified per region, 0 donors at stage I |
| Allen ADTBI (bulk) | 0, 1, 2, 3, 4, 5, 6 | **yes — 40 of 377 samples** |

**AD group definition.** Where a binary contrast is used (script 01 only), the AD group is
**Braak >= 2**. Braak I is deliberately excluded from both groups: transentorhinal-only tau is
common in cognitively normal aging and is not treated as AD here. The Leng cohort contains no
Braak I donors, so this threshold is behaviourally identical to the `> 0` test used in the v3.0
release, but the script now asserts the absence of Braak I rather than relying on it, so
re-running on a cohort that does contain Braak I fails loudly instead of misclassifying it.

ADTBI, the only cohort containing Braak I, is analysed with continuous Spearman correlations
across all seven stages, so Braak I enters as stage 1 and no dichotomy applies.

## Derived files (produced by script 04, then consumed by others)

- `data/cellclass_all_regions.csv` — per-region × cell-class interaction ORs. **Produced by 04**, consumed by 02 and 05.
- `data/staged_ec_all_regions.json` — SEA-AD staged per-region results. **Produced by 04**, consumed by 14.

Run `04_staged_ec_replication.py` first if you need these; scripts 02, 05, 14
depend on its outputs (copy them into `data/`).

## Values inlined in scripts (no file needed)

- Script 06 (APOE) uses the isoform-group ORs computed in 04's APOE panel; the
  values are inlined in the script with their provenance noted. The upstream
  computation is in 04.

## Access terms

All primary datasets are public. SEA-AD and Allen ADTBI are released under the
Allen Institute terms of use; CELLxGENE Census data retain their original
dataset licenses. No controlled-access (dbGaP/Synapse) data are required to run
any script here — the proposed genetic-epistasis experiment that *would* need
controlled access is described in the manuscript but not implemented in this repo.
