"""
Cross-disorder dataset provenance (CELLxGENE Census)
====================================================

The cross-disorder comparison (Sec 3.11) draws on many contributed studies pooled
through CELLxGENE Census. Census is a standardisation and distribution layer: it
harmonises schema, ontology terms and cell-type labels, but it does NOT re-run
quality control on contributor data. QC provenance therefore devolves to each
constituent study, and a claim about QC for this row is only citable if the
constituent datasets can be named.

This script enumerates them. It resolves every dataset contributing nuclei to the
disease x region strata used in Sec 3.11, together with its dataset_id, source
collection and collection DOI, so each can be cited individually.

It also re-derives the coverage claim the section rests on: among all brain
disorders in Census, only Alzheimer's disease has entorhinal-cortex data and only
epilepsy has hippocampal-formation data. That claim is a property of the Census
build, so it is re-checked here rather than asserted from a previous run.

NOTE ON REPRODUCIBILITY -- Census is versioned. Pin CENSUS_BUILD to the build used
in the manuscript; a later build will contain more datasets and the coverage claim
may legitimately change.

Inputs : streamed from CELLxGENE Census (no local file)
Outputs: census_cross_disorder_provenance.csv  (one row per disease x region x dataset)
         census_source_collections.csv         (one row per source collection)
Method : Census dataset registry joined to obs-level disease/tissue/donor counts

Part of the RELN-RORB EC convergence analysis. See README.md.
"""
import numpy as np
import pandas as pd
import cellxgene_census as cxg

# Census build used in the manuscript. Do not float this to "latest".
CENSUS_BUILD = "2025-11-08"

# Diseases surveyed for allocortical coverage in Sec 3.11.
DISEASES = [
    "Alzheimer disease", "epilepsy", "schizophrenia",
    "progressive supranuclear palsy", "Lewy body dementia",
    "amyotrophic lateral sclerosis", "Parkinson disease",
    "frontotemporal dementia", "dementia", "tauopathy", "normal",
]

# Regions the section compares: allocortex (the region-matched axis) and
# neocortex (the disorder-matched axis).
REGION_PAT = "entorhinal|hippocamp|prefrontal"

# Strata below this many nuclei are noise in a coverage table, not cohorts.
MIN_NUCLEI = 200


def open_census(build=CENSUS_BUILD):
    """Open the pinned Census build.

    Behind a proxy, TileDB's S3 client needs to be told about it explicitly --
    it uses raw sockets and will otherwise fail DNS resolution rather than
    fall back to HTTP(S)_PROXY the way most Python HTTP clients do.
    """
    import os
    from urllib.parse import urlparse
    proxy = (os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
             or os.environ.get("HTTP_PROXY") or os.environ.get("http_proxy"))
    if not proxy:
        return cxg.open_soma(census_version=build)
    import tiledbsoma
    u = urlparse(proxy)
    ctx = tiledbsoma.SOMATileDBContext(tiledb_config={
        "vfs.s3.proxy_host": u.hostname,
        "vfs.s3.proxy_port": str(u.port or 80),
        "vfs.s3.proxy_scheme": u.scheme,
        "vfs.s3.region": "us-west-2",
        "vfs.s3.no_sign_request": "true",
    })
    return cxg.open_soma(census_version=build, context=ctx)


def load_strata(census):
    """Per-nucleus disease/tissue/donor/dataset labels for the surveyed diseases.

    is_primary_data == True is essential: Census carries the same nuclei under
    multiple collections, and without it every count is inflated by duplication.
    """
    flt = " or ".join(f"disease == '{d}'" for d in DISEASES)
    obs = (census["census_data"]["homo_sapiens"].obs
           .read(value_filter=f"({flt}) and is_primary_data == True",
                 column_names=["dataset_id", "disease", "tissue", "donor_id"])
           .concat().to_pandas())
    return obs[obs.tissue.str.contains(REGION_PAT, case=False, na=False)]


def coverage(strata):
    """disease x region coverage -- the table behind the Sec 3.11 coverage claim."""
    return (strata.groupby(["disease", "tissue"], observed=True)
            .agg(nuclei=("dataset_id", "size"),
                 datasets=("dataset_id", "nunique"),
                 donors=("donor_id", "nunique"))
            .reset_index().sort_values("nuclei", ascending=False))


def provenance(strata, registry):
    """One row per disease x region x dataset, carrying the citation."""
    rows = []
    for (dis, tis), grp in strata.groupby(["disease", "tissue"], observed=True):
        for did, g in grp.groupby("dataset_id", observed=True):
            if len(g) < MIN_NUCLEI:
                continue
            hit = registry[registry.dataset_id == did]
            if hit.empty:
                continue
            d = hit.iloc[0]
            rows.append(dict(
                disease=dis, tissue=tis, nuclei=len(g),
                donors=g.donor_id.nunique(), dataset_id=did,
                dataset_title=d.dataset_title, collection=d.collection_name,
                collection_doi=d.collection_doi,
                citation_short=d.collection_doi_label))
    return pd.DataFrame(rows).sort_values(
        ["disease", "tissue", "nuclei"], ascending=[True, True, False])


def check_coverage_claim(cov):
    """Re-derive the Sec 3.11 coverage claim from this build."""
    def diseases_in(pat):
        q = cov[cov.tissue.str.contains(pat, case=False, na=False)
                & (cov.disease != "normal") & (cov.nuclei >= MIN_NUCLEI)]
        return sorted(q.disease.unique())
    ec = diseases_in("entorhinal")
    hpf = diseases_in("hippocamp")
    print("\nCoverage claim (disorders only, 'normal' excluded):")
    print(f"  entorhinal cortex     : {ec}")
    print(f"  hippocampal formation : {hpf}")
    if ec != ["Alzheimer disease"]:
        print("  WARNING: entorhinal coverage differs from the manuscript build.")
    if hpf != ["epilepsy"]:
        print("  WARNING: hippocampal coverage differs from the manuscript build.")
    return ec, hpf


def main():
    census = open_census()
    try:
        registry = census["census_info"]["datasets"].read().concat().to_pandas()
        print(f"Census build {CENSUS_BUILD}: {len(registry)} datasets in registry")
        strata = load_strata(census)
        print(f"primary-data nuclei in surveyed regions: {len(strata):,}")

        cov = coverage(strata)
        print("\n=== disease x region coverage (>= 500 nuclei) ===")
        print(cov[cov.nuclei >= 500].to_string(index=False))

        check_coverage_claim(cov)

        prov = provenance(strata, registry)
        n_doi = prov.collection_doi.replace("", np.nan).notna().sum()
        print(f"\n=== provenance: {prov.dataset_id.nunique()} datasets, "
              f"{prov.collection_doi.nunique()} collections, "
              f"{n_doi}/{len(prov)} rows with a published DOI ===")
        print(prov[prov.tissue.str.contains("entorhinal|hippocamp", case=False)]
              [["disease", "tissue", "nuclei", "donors", "citation_short"]]
              .to_string(index=False))

        coll = (prov.groupby(["citation_short", "collection_doi"], observed=True)
                .agg(datasets=("dataset_id", "nunique"), nuclei=("nuclei", "sum"))
                .reset_index().sort_values("nuclei", ascending=False))

        prov.to_csv("census_cross_disorder_provenance.csv", index=False)
        coll.to_csv("census_source_collections.csv", index=False)
        print("\nwrote census_cross_disorder_provenance.csv, census_source_collections.csv")
    finally:
        census.close()


if __name__ == "__main__":
    main()
