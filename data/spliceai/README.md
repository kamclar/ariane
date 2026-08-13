# ARIANE SpliceAI Data

The binding classification profile is
`enigma_v1_2_spliceai_profile.json`. It implements ENIGMA BRCA1/2 VCEP v1.2
Appendix J with GRCh38, maximum distance 10,000, unmasked output, the reference
transcript, and the maximum of DS_AG, DS_AL, DS_DG, and DS_DL. Every accepted
record also retains the four REF and four ALT component scores.

## Current API-primary mode

Until both Appendix J caches have been completely rebuilt, ARIANE keeps
`SPLICEAI_USE_PRECOMPUTED_CACHE=0`. Runtime uses the profile-pinned Broad API
cache and then the configured Broad-compatible API. Legacy coding and intronic
files remain build inputs and audit material, but are not active classification
sources and do not create degraded-source warnings.

After a complete checksum-validated rebuild, the immutable caches can be
activated explicitly with `SPLICEAI_USE_PRECOMPUTED_CACHE=1`.

Both coding and intronic caches are built with
`scripts/build_spliceai_reference_caches.py` against local instances of the
container image pinned by digest in the profile.

Companion metadata:

- `spliceai_brca_snv_reference_cache.metadata.json`

Final report:

- `variant_space_scan/docs/final_spliceai_precompute_report.md`

## Lookup Order

For the current API-primary `reference_transcript` policy:

1. In-memory cache
2. `spliceai_api_cache.json`
3. A configured Broad-compatible API call using exactly the same profile

The immutable precomputed files are skipped while
`SPLICEAI_USE_PRECOMPUTED_CACHE=0`.

There is no classification mode that silently switches to a maximum across
other transcripts. An environment request for a conflicting transcript policy
stops startup.

Old runtime records use a different key and are ignored. A response that does
not echo GRCh38, distance 10,000 and mask 0, or lacks delta/REF/ALT fields, is
rejected. Missing scores remain unavailable and are never converted to zero.

## Reproducible Build

The builder writes resumable checkpoints under `data/spliceai/build/`. A
checkpoint is accepted only when its scoring profile, source checksum and own
checksum match. Production JSON and metadata are replaced only after all
expected records succeed. The existing production cache is never used as a
resume source.

```powershell
python scripts\build_spliceai_reference_caches.py all `
  --api-url http://127.0.0.1:8080/spliceai/ `
  --workers 3 --timeout 240 --delay 0
```

## Maintenance

The public Broad SpliceAI API can change and explicitly disallows batch use.
Reference builds therefore use only local instances of the pinned image.
Independent spot checks may use the public service at an interactive rate and
must record API failures separately from numeric differences.

Recommended cadence:

- Quarterly
- Before clinical/release updates
- After any known Broad API or SpliceAI model update
- After changing transcript policy

Recommended check:

- Random sample of 80 to 100 cached variants
- Slow requests with retry
- Compare local cache score against the public Broad API score filtered to the
  same reference transcript
- Record failures separately from numeric mismatches

The expected outcome is zero numeric mismatches among successful public API
responses. Transient public API failures should be retried before interpreting
them.

Example command from the repository root:

```powershell
python variant_space_scan\validate_spliceai_reference_pilot.py --cache data\spliceai\spliceai_brca_snv_reference_cache.json --report variant_space_scan\outputs\spliceai_app_cache_vs_broad_reference_validation.YYYYMMDD.json --sample-size 100 --sleep-seconds 5 --retries 2 --retry-sleep-seconds 20 --timeout 90
```

The same JSON report can be converted to CSV with:

```powershell
python variant_space_scan\export_spliceai_validation_csv.py variant_space_scan\outputs\spliceai_app_cache_vs_broad_reference_validation.YYYYMMDD.json variant_space_scan\outputs\spliceai_app_cache_vs_broad_reference_validation.YYYYMMDD.csv
```
