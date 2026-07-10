# ARIANE SpliceAI Data

## Active Offline Cache

ARIANE uses `spliceai_brca_snv_reference_cache.json` as the primary SpliceAI
source for BRCA1/BRCA2 coding SNVs under the `reference_transcript` policy.

The cache was produced by `variant_space_scan` using the patched local
Broad-compatible Docker endpoint and validated against the public Broad API.

Companion metadata:

- `spliceai_brca_snv_reference_cache.metadata.json`

Final report:

- `variant_space_scan/docs/final_spliceai_precompute_report.md`

## Lookup Order

For `SPLICEAI_TRANSCRIPT_POLICY=reference_transcript`:

1. In-memory cache
2. `spliceai_brca_snv_reference_cache.json`
3. `spliceai_api_cache.json`
4. Public or configured Broad-compatible API fallback

For `SPLICEAI_TRANSCRIPT_POLICY=max_any_transcript`, the reference-transcript
offline cache is skipped because it does not contain max-across-transcripts
scores.

## Maintenance

The public Broad SpliceAI API can change if model, annotation, transcript, or
deployment details change. The offline cache should therefore be spot-checked
periodically.

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
