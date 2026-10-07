"""
Shared data-path resolution and preflight checks.
=================================================

Every analysis script reads its inputs from a local `data/` directory that the
user must populate first -- the files are large and several are access-controlled,
so they are deliberately not in the repository (see DATA.md).

Without a preflight check, a missing file surfaces as a 20-line pandas or h5py
traceback that says nothing about where to get the file. `require()` turns that
into one actionable line naming the file, its size, and the DATA.md section that
documents its source.

Usage:

    from _data_paths import DATA, require

    MEC, LENG = require("seaad_mec_extract.parquet", "leng2021_ec_full.h5ad")
    df = pd.read_parquet(MEC)

Set RELN_RORB_DATA to point somewhere other than ./data:

    RELN_RORB_DATA=/scratch/reln_rorb/data python analysis/16_robustness_controls.py
"""
import os
import sys

DATA = os.environ.get("RELN_RORB_DATA", "data")

# filename -> (approx size, DATA.md heading that documents how to obtain it)
_PROVENANCE = {
    "seaad_mec_extract.parquet":  ("~4 MB",   "SEA-AD multi-region snRNA-seq per-nucleus extracts"),
    "seaad_lec_extract.parquet":  ("~1 MB",   "SEA-AD multi-region snRNA-seq per-nucleus extracts"),
    "seaad_hip_extract.parquet":  ("~1 MB",   "SEA-AD multi-region snRNA-seq per-nucleus extracts"),
    "seaad_v1c_extract.parquet":  ("~3 MB",   "SEA-AD multi-region snRNA-seq per-nucleus extracts"),
    "leng2021_ec_full.h5ad":      ("~125 MB", "`data/leng2021_ec_full.h5ad`  — Leng et al. 2021 entorhinal cortex"),
    "seaad_merfish_mec.h5ad":     ("~76 MB",  "`data/seaad_merfish_mec.h5ad` — SEA-AD MERFISH, medial entorhinal"),
    "seaad_merfish_hpf.h5ad":     ("~108 MB", "`data/seaad_merfish_hpf.h5ad` — SEA-AD MERFISH, hippocampal formation"),
    "rod_reln_rorb.parquet":      ("~2 MB",   "`data/rod_reln_rorb.parquet` — Rodriguez-Rodriguez et al. 2025 onset cohort"),
    "oterogarcia_tangle_cells.parquet": ("~1 MB", "`data/oterogarcia_tangle_cells.parquet` — Otero-Garcia et al. tangle-sorted somata"),
    "adtbi_merged.parquet":       ("~2 MB",   "`data/adtbi_merged.parquet` — Allen Aging, Dementia & TBI bulk RNA-seq"),
    "Exc_neu_integrated_with_Leng.rds": ("~1.5 GB", "`data/rod_reln_rorb.parquet` — Rodriguez-Rodriguez et al. 2025 onset cohort"),
    "seaad2026_MEC_final.parquet": ("~19 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_LEC_final.parquet": ("~4 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_HIP_final.parquet": ("~4 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_ITG_final.parquet": ("~8 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_MTG_final.parquet": ("~24 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_STG_final.parquet": ("~11 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_FI_final.parquet": ("~5 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_AnG_final.parquet": ("~5 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_PFC_final.parquet": ("~27 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_V1C_final.parquet": ("~13 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "seaad2026_MEC_allnuclei.parquet": ("~19 MB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
    "Global_and_Local_CPS.csv": ("~55 KB", "SEA-AD Multiregion 2026 release (ten regions) — produced by script 25"),
}


def require(*names):
    """Resolve input filenames under DATA, exiting with guidance if any are absent.

    Returns one path per name, in order, so a single call can unpack:
        MEC, LENG = require("seaad_mec_extract.parquet", "leng2021_ec_full.h5ad")
    """
    paths = [os.path.join(DATA, n) for n in names]
    missing = [(n, p) for n, p in zip(names, paths) if not os.path.exists(p)]
    if missing:
        lines = [
            "",
            f"Missing {len(missing)} of {len(names)} required input file(s) in '{DATA}/'.",
            "",
        ]
        for n, p in missing:
            size, section = _PROVENANCE.get(n, ("unknown size", "DATA.md"))
            lines.append(f"  {n}")
            lines.append(f"      expected at : {p}")
            lines.append(f"      size        : {size}")
            lines.append(f"      how to get  : DATA.md, section '{section}'")
        lines += [
            "",
            "These files are large and several are access-controlled, so they are not",
            "in the repository. DATA.md documents the source and the reduction step for",
            "each one. If your copies live elsewhere, point RELN_RORB_DATA at them:",
            "",
            f"    RELN_RORB_DATA=/path/to/data python {os.path.basename(sys.argv[0] or 'analysis/<script>.py')}",
            "",
        ]
        sys.exit("\n".join(lines))
    return paths[0] if len(paths) == 1 else tuple(paths)
