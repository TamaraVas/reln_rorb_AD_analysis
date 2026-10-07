"""
SEA-AD Multiregion 2026 extraction: two genes + per-nucleus QC from ten regions
==============================================================================

The 2026-06-22 multi-region release is ~161 GB of h5ad across ten regions. We need
two gene columns and a dozen obs columns from each, so this script downloads one
object at a time, extracts, writes a small parquet, and DELETES the h5ad. Peak disk
is one object (max 39 GB), not the whole release.

Why not slice remotely: both X and layers/UMIs are CSR with ~5.4e9 nonzeros, so
pulling two GENE columns means reading essentially the whole index array. Streaming
the file and extracting locally is cheaper. (An fsspec/aiohttp-backed h5py reader
also bypasses proxied networks; a plain Range-request reader is used for the small
metadata reads instead.)

FILE SELECTION IS A TRAP. MTG/RNAseq/ contains TWO *_final-nuclei.h5ad objects: the
2026 SEA-AD cohort (32.98 GB) and a 2022 neurotypical Reference (6.38 GB). Taking
the first match silently picks the Reference, which has no AD staging -- it would
enter a regional comparison as a wrong-cohort ringer. Every key is therefore checked
against REQUIRE_IN_NAME and refused otherwise.

Both normalised and raw-UMI values are captured for the two genes so either
positivity rule can be applied downstream without re-downloading:
  - manuscript SEA-AD rule : normalised value >= 1   (RELN+ 13.3%, RORB+ 20.2% in MEC)
  - raw-count rule (Leng)  : UMI count >= 1          (36.6% / 44.3% in MEC -- NOT comparable)

Inputs : none on disk -- streams from the public SEA-AD S3 bucket
Outputs: data/seaad2026_<REGION>_final.parquet (x10), data/seaad2026_MEC_allnuclei.parquet,
         data/Global_and_Local_CPS.csv
Method : HTTP Range reads for metadata; full download + h5py extraction per object

Part of the RELN-RORB EC convergence analysis. See README.md and DATA.md.
"""
import io
import os
import re
import sys
import urllib.parse
import urllib.request

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/"
RELEASE = "2026-06-22"
# Both tokens must appear in every filename we accept. See the docstring.
REQUIRE_IN_NAME = ("SEAAD_", RELEASE)
GENES = ["RELN", "RORB"]
REGIONS = ["MEC", "LEC", "HIP", "ITG", "MTG", "STG", "FI", "AnG", "PFC", "V1C"]
# The prefrontal data sits under PFC/ and is labelled DFC in the per-nucleus
# metadata; the DFC/ and DLPFC/ prefixes exist but are empty placeholders.
OBS_KEEP = [
    "Donor ID", "Braak", "Thal", "CERAD score", "ADNC", "Cognitive status",
    "APOE Genotype", "Class", "Subclass", "Supertype", "Number of UMIs",
    "Genes detected", "Fraction mitochondrial UMIs", "Doublet score",
    "Used in analysis", "Severely Affected Donor", "LATE",
]
CPS_KEY = ("Multiregion_2026/model_outputs/continuous_pseudo-progression_score/"
           "Global_and_Local_CPS.20260105.csv")
OUT = os.environ.get("RELN_RORB_DATA", "data")


class S3Range(io.RawIOBase):
    """Minimal seekable reader over HTTP Range requests (proxy-safe)."""

    def __init__(self, url):
        self.url, self._pos = url, 0
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            self.size = int(r.headers["Content-Length"])

    def readable(self):
        return True

    def seekable(self):
        return True

    def seek(self, off, whence=0):
        self._pos = off if whence == 0 else (self._pos + off if whence == 1 else self.size + off)
        return self._pos

    def tell(self):
        return self._pos

    def read(self, n=-1):
        if n < 0 or self._pos + n > self.size:
            n = self.size - self._pos
        if n <= 0:
            return b""
        req = urllib.request.Request(
            self.url, headers={"Range": "bytes=%d-%d" % (self._pos, self._pos + n - 1)})
        with urllib.request.urlopen(req, timeout=120) as r:
            b = r.read()
        self._pos += len(b)
        return b


def lsdir(prefix):
    """One page of keys under a prefix. Delimiter-based: do NOT paginate the whole
    bucket to enumerate regions -- it holds >500,000 objects and a capped loop
    silently under-reports which regions exist."""
    url = BASE + "?list-type=2&delimiter=/&prefix=%s&max-keys=1000" % urllib.parse.quote(prefix)
    with urllib.request.urlopen(url, timeout=60) as r:
        x = r.read().decode()
    keys = re.findall(r"<Key>([^<]+)</Key>", x)
    sizes = [int(s) for s in re.findall(r"<Size>(\d+)</Size>", x)]
    return list(zip(keys, sizes))


def check_key(key):
    name = key.split("/")[-1]
    bad = [t for t in REQUIRE_IN_NAME if t not in name]
    if bad:
        raise SystemExit("refusing %s: missing %s (wrong release or Reference object)" % (name, bad))
    return key


def plan():
    """Resolve one final-nuclei key per region, plus MEC all-nuclei."""
    rows = []
    for rg in REGIONS:
        prefix = "PFC/RNAseq/" if rg == "PFC" else "%s/RNAseq/" % rg
        ks = lsdir(prefix)
        fin = [(k, s) for k, s in ks
               if all(t in k for t in REQUIRE_IN_NAME) and "final-nuclei" in k and k.endswith(".h5ad")]
        aln = [(k, s) for k, s in ks
               if all(t in k for t in REQUIRE_IN_NAME) and "all-nuclei" in k and k.endswith(".h5ad")]
        if len(fin) != 1:
            raise SystemExit("%s: expected exactly 1 final-nuclei object, found %d" % (rg, len(fin)))
        rows.append(dict(region=rg, key=fin[0][0], gb=fin[0][1] / 1e9,
                         all_key=aln[0][0] if aln else None,
                         all_gb=aln[0][1] / 1e9 if aln else 0.0))
    return pd.DataFrame(rows)


def obs_col(h, name):
    o = h["obs"][name]
    if isinstance(o, h5py.Group):
        cats = np.array([s.decode() if isinstance(s, bytes) else s for s in o["categories"][:]])
        return pd.Categorical.from_codes(o["codes"][:], cats).astype(str)
    v = o[:]
    return np.array([s.decode() if isinstance(s, bytes) else s for s in v]) if v.dtype.kind in "OS" else v


def gene_cols(h, group):
    """Two gene columns out of a CSR matrix, read once row-block-wise."""
    idx = np.array([str(s.decode() if isinstance(s, bytes) else s) for s in h["var"]["index"][:]])
    want = {g: int(np.where(idx == g)[0][0]) for g in GENES if (idx == g).any()}
    if len(want) != len(GENES):
        raise SystemExit("genes not found in var index: %s" % (set(GENES) - set(want)))
    g = h[group]
    n = int(dict(g.attrs)["shape"][0])
    data, indices, indptr = g["data"], g["indices"], g["indptr"]
    out = {k: np.zeros(n, dtype=np.float32) for k in want}
    step = 200_000
    for start in range(0, n, step):
        stop = min(start + step, n)
        lo, hi = int(indptr[start]), int(indptr[stop])
        if hi <= lo:
            continue
        block = sp.csr_matrix((data[lo:hi], indices[lo:hi], indptr[start:stop + 1] - lo),
                              shape=(stop - start, len(idx)))
        for k, j in want.items():
            out[k][start:stop] = np.asarray(block[:, j].todense()).ravel()
    return out


def extract(key, tag):
    dest = os.path.join(OUT, os.path.basename(key))
    out_pq = os.path.join(OUT, "seaad2026_%s.parquet" % tag)
    if os.path.exists(out_pq):
        print("  %s already extracted, skipping" % tag, flush=True)
        return
    check_key(key)
    print("  downloading %s ..." % key.split("/")[-1], flush=True)
    urllib.request.urlretrieve(BASE + urllib.parse.quote(key), dest)
    try:
        with h5py.File(dest, "r") as h:
            df = pd.DataFrame({c: obs_col(h, c) for c in OBS_KEEP if c in h["obs"]})
            norm = gene_cols(h, "X")
            umi = gene_cols(h, "layers/UMIs") if "layers" in h and "UMIs" in h["layers"] else {}
            for g in GENES:
                df["%s_norm" % g] = norm[g]
                if g in umi:
                    df["%s_umi" % g] = umi[g]
        df.to_parquet(out_pq, index=False)
        print("  %s -> %d rows, %d cols" % (tag, len(df), df.shape[1]), flush=True)
    finally:
        if os.path.exists(dest):
            os.remove(dest)


def main():
    os.makedirs(OUT, exist_ok=True)
    cps = os.path.join(OUT, "Global_and_Local_CPS.csv")
    if not os.path.exists(cps):
        urllib.request.urlretrieve(BASE + urllib.parse.quote(CPS_KEY), cps)
        print("CPS scores -> %s (%d bytes)" % (cps, os.path.getsize(cps)), flush=True)
    p = plan()
    p.to_csv(os.path.join(OUT, "seaad2026_download_plan.csv"), index=False)
    print(p[["region", "gb", "all_gb"]].to_string(index=False), flush=True)
    print("total final-nuclei: %.1f GB" % p.gb.sum(), flush=True)
    for _, r in p.iterrows():
        extract(r["key"], "%s_final" % r["region"])
    mec = p[p.region == "MEC"].iloc[0]
    if isinstance(mec["all_key"], str):
        # NOTE: the per-region all-nuclei objects have IDENTICAL cell counts to
        # final-nuclei (MEC 970,118 both; Used-in-analysis 100% True), so expression
        # for the ~2.6M QC-FAILED nuclei is NOT distributed. Kept for completeness.
        extract(mec["all_key"], "MEC_allnuclei")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
