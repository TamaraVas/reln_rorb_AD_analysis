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
RELN⁺, RORB⁺, log-depth, Braak, CERAD, APOE genotype, cell class/subclass.
- **Source:** SEA-AD multi-region 10x snRNA-seq (2026 release), streamed per-donor
  from the public bucket `s3://sea-ad-single-cell-profiling/` and reduced to the
  RELN/RORB/depth/pathology columns. The reduction script is `04`'s data-prep
  (documented in its docstring); the raw per-donor h5ad objects are the ABC Atlas
  SEA-AD-Multiregion-10X release.
- Regions/donors: MEC 81, LEC 34, HIP 40, V1C 43.
- Used by script 04 (which then writes the two derived files below).

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
