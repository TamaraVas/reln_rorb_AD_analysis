# Changelog

## v3.0 — first release since v1

**If you cloned v1, this is the release to read.** v2.0 and v2.1 were prepared but
never published, so v3.0 carries their contents plus later work. Cumulative diff
against v1: 9 files added, 4 modified, none removed; the other 19 are byte-identical.
The README's "Changes since the first public release" section is the narrative
version of this entry.

### Added
- `analysis/21_oligodendrocyte_specificity.py` — RELN⁺ oligodendrocytes rise with
  Braak stage (OR 1.21/step, p = 0.004; 95,556 nuclei, 81 donors), but the script's
  purpose is to bound that observation rather than advance it. Class-specific among
  the six non-neuronal subclasses (all others 0.88–0.96, which is the argument against
  ambient RNA); **brain-wide, not entorhinal** (MEC 1.21, LEC 1.19, HIP 1.19, V1C 1.30);
  RELN-only, with the within-oligodendrocyte interaction null; tracking general
  neuropathological burden rather than tau specifically, with dementia the largest
  single effect (OR 1.76). Reported in §3.6 as a specificity control.
- `analysis/22_directionality_replication.py` — documents why the manuscript makes no
  directional claim. Reports both conditional acquisition slopes and the convergence
  interaction side by side in both cohorts: the slopes' *ordering* reverses between
  cohorts while the convergence replicates in both.
- `CHANGELOG.md` (this file, added at v2.0).
- Scripts 16–20 (added at v2.0/v2.1 — see entries below).

### Changed
- `README.md` — the *Changes since the first public release* section is new and written
  as a single cumulative v1→v3 diff, since no intermediate release was published.
  Script table extended to 22 entries with the manuscript section each supports.
- `METHODS_PROVENANCE.md` — entries for scripts 21–22.
- `DATA.md` — the four SEA-AD regional extracts (MEC/LEC/HIP/V1C) are now listed
  individually, since scripts 21–22 read more than MEC.

### Fixed (found by validating the new scripts against known-good data)
- Script 21 coded dementia status with a substring test, which scored the negative
  level `"No dementia"` as positive because it *contains* "dementia" — every cell was
  marked positive, inflating the odds ratio from 1.76 to 2.26. Now an exact match with
  an assertion on the level set.
- Script 22 contained a malformed model formula (`bp ~ rp * 0 + bkn + logd`) in the
  RORB-acquisition fit. Corrected and verified to reproduce 1.545 (p = 5.6×10⁻⁵).
- Stripped a `__pycache__` directory that the v2.1 tarball had captured.

### Withdrawn claims (documented, not silently dropped)
- The directional refinement ("a stable RORB⁺ identity acquiring RELN"): not
  identifiable in cross-sectional autopsy tissue, and the direction reverses between
  cohorts. See script 22.
- A pooled cross-cohort effect magnitude: between-cohort heterogeneity too high.
  `14_cross_cohort_meta.py` reports per-cohort estimates instead.

### Added after QC review
- `analysis/24_merfish_detection_controls.py` — measures MERFISH detection error from the control
  channels and separates decoding error (attenuating; corrected OR 2.47–2.55 against observed 2.45)
  from segmentation spillover (tested by segmentation method and by cell area).
- `analysis/23_census_dataset_provenance.py` — enumerates the 37 datasets across 21 source
  collections behind the cross-disorder comparison, with collection DOIs, and re-derives the
  allocortex coverage claim from the pinned Census build. Closes a provenance gap: that row
  previously read "many (public atlases)" with no dataset identifiers recorded anywhere.
- Manuscript §3.8 now reports a measured MERFISH detection-error floor (all control channels
  = 0.112% of signal; genomic controls exactly zero) and shows the noise floor cannot generate
  the co-localization (0.69 expected double-positive cells from noise against 928 observed).
- Results sections renumbered 3.1–3.13 (a gap at 3.13 remained from the earlier §3.12 relocation).

### Corrected attribution (pre-push audit)
- `analysis/18_roussarie_reciprocal.py` renamed to `analysis/18_rodriguez_reciprocal.py`.
  The cohort is Rodriguez-Rodriguez et al. 2025 (GEO GSE287652); "Roussarie" was an
  early working label that survived into the filename. The manuscript always cited the
  correct authors, and `DATA.md` was already correct, but the repository was not.
  **If you cloned v1 this file did not exist, so no link breaks; anyone tracking the
  unpublished v2.x tarballs should note the path change.**
- `analysis/_data_paths.py` provenance pointers now quote `DATA.md`'s literal section
  headings. Four of nine previously named sections that did not exist, so the guidance
  sent a reader looking for a heading they would not find.

### Verified by cold-run smoke test
- All scripts with obtainable inputs were executed as subprocesses from a clean
  `data/` directory and reproduce the manuscript values (see README, *Verification
  state of this release*). `18_rodriguez_reciprocal.py` remains partially unverified:
  its input derives from a ~1.5 GB non-redistributable Seurat object.
- Added `analysis/_data_paths.py`. Before this, a missing input produced a bare
  pandas/h5py traceback naming no provenance; all seven scripts now exit with the
  filename, its size, and the `DATA.md` section that documents its source. Honours
  `RELN_RORB_DATA` for inputs stored outside the repository.

## v2.1

- Added `analysis/20_second_rorb_cluster.py`: the RORB-alone negative control reported in §3.2 and Supplementary Fig. S10. Tests whether RORB expression on its own predicts RELN co-expression, using the second RORB-high entorhinal cluster (EC:Exc.5) as an internal control, with a depth-confound check and SEA-AD replication.
- Documented the <10-double-positive-cell separation threshold (`MIN_DP`) used to mark subclass estimates unestimable rather than quoting them.

## v2.0 — revision round

Documentation corrections and four added analyses from peer review. **No previously
published result changed**; the corrections are to how methods are described and how
significance is reported.

### Fixed — positivity rule was misdocumented (important)

The README previously stated a single rule, "a gene is scored positive if its raw
count ≥ 1." That is true of the Leng and MERFISH analyses and **false of every
SEA-AD analysis**. The `reln`/`rorb` columns of the SEA-AD extracts hold
depth-normalised, truncated expression values, not raw UMI counts — they are small
integers (0–5) and so read like counts. Script 04's `> 0` test is therefore a
normalised-expression threshold, roughly ≥7 raw transcripts at that cohort's median
depth.

The code always did this; only the documentation was wrong. It is now stated in
three places: the README "Statistical approach" section, a warning block in DATA.md,
and an inline note at the thresholding lines in script 04. A reader who wants to
verify: under a genuine raw-count rule detection rises steeply with library size,
whereas in these columns RELN detection is flat across depth deciles and RORB
detection *falls* (30.8% → 14.2%).

### Changed — donor-clustered standard errors are the reporting standard

Cell-level p-values treat each nucleus as independent when the experimental unit is
the donor. Odds ratios are unaffected; p-values move (e.g. the SEA-AD MEC
interaction: OR 1.481 either way, p 1.8×10⁻¹⁶⁴ cell-level → 1.5×10⁻³
donor-clustered). Script 04 already used clustered SEs; the new scripts follow suit.

### Added — four analysis scripts

| Script | Purpose |
|---|---|
| `16_robustness_controls.py` | Depth thinning to common library size, ambient-floor and ≥2-transcript threshold, leave-one-donor-out, donor-label permutation |
| `17_clustering_circularity.py` | Subtype enrichment re-run with RELN/RORB removed from the clustering feature set |
| `18_rodriguez_reciprocal.py` | Same-cell test applied to an independent Braak 0–II onset cohort (GEO GSE287652) |
| `19_within_donor_regional.py` | MEC-versus-hippocampus spatial contrast within a single donor |

### Added — data documentation

- **GSE287652** (independent onset cohort) with the two handling caveats that affect
  the analysis: the published object is integrated with our own discovery cohort and
  must be subset; the original authors' sample exclusions are honoured.
- **MERFISH donor provenance.** The data generators confirmed that both spatial
  sections come from one donor (H24.30.005), profiled in a **neurotypical** context
  and unlikely to carry significant pathological burden, and that MEC/HIP spatial
  data across levels of pathological burden have not yet been generated. The spatial
  analysis therefore establishes anatomical specificity, not a disease gradient.

### Housekeeping

- Removed macOS AppleDouble (`._*`) files that were present in the previous archive.
- `METHODS_PROVENANCE.md` extended with entries for scripts 16–19.

---

## v1.0 — initial release

Fifteen analysis scripts covering the same-cell co-occupancy model, regional and
cell-class specificity, Braak trajectory, SEA-AD staged replication, APOE isoform
panel, spatial co-localization, tangle-bearing neurons, gene-regulatory-network
control, bulk pathology, cross-tauopathy comparison, marker-pair specificity,
reelin receptor arm, cross-cohort meta-analysis, and a squidpy cross-check.
