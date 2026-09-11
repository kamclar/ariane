# ARIANE SpliceAI Data

The binding classification profile is
`enigma_v1_2_spliceai_profile.json`. It implements ENIGMA BRCA1/2 VCEP v1.2
Appendix J with GRCh38, maximum distance 10,000, unmasked output, the reference
transcript, and the maximum of DS_AG, DS_AL, DS_DG, and DS_DL. Every accepted
record also retains the four REF and four ALT component scores.

## Runtime mode

ARIANE computes SpliceAI on demand. Runtime uses the profile-pinned result cache
and then the local Broad-compatible service on `127.0.0.1:8082`. Precomputed coding and intronic
variant spaces are not classification sources and cannot be activated by an
environment variable. Existing files remain only as historical validation and
audit material.

Historical coding and intronic comparison datasets can be rebuilt with
`scripts/build_spliceai_reference_caches.py` against a local instance of the
container image pinned by digest in the profile. They are not required for
runtime classification.

Companion metadata:

- `spliceai_brca_snv_reference_cache.metadata.json`

Final report:

- `variant_space_scan/docs/final_spliceai_precompute_report.md`

## Lookup Order

For the current local on-demand `reference_transcript` policy:

1. In-memory cache
2. `${ARIANE_RUNTIME_CACHE_DIR}/spliceai_api_cache.json`, or
   `.runtime-cache/spliceai_api_cache.json` in local development
3. The configured local Broad-compatible service using exactly the same profile

The lookup does not read a precomputed gene-wide variant space.

There is no classification mode that silently switches to a maximum across
other transcripts. An environment request for a conflicting transcript policy
stops startup.

Old runtime records use a different key and are ignored. A response that does
not echo GRCh38, distance 10,000 and mask 0, or lacks delta/REF/ALT fields, is
rejected. Missing scores remain unavailable and are never converted to zero.

## Historical validation datasets

The builder can create comparison datasets under `data/spliceai/build/`. A
checkpoint is accepted only when its scoring profile, source checksum and own
checksum match. JSON and metadata are replaced only after all expected records
succeed. These datasets do not enter the runtime lookup order.

```powershell
python scripts\build_spliceai_reference_caches.py all `
  --api-url http://127.0.0.1:8080/spliceai/ `
  --workers 3 --timeout 240 --delay 0
```

## Maintenance

The public Broad SpliceAI API can change and explicitly disallows batch use.
ARIANE therefore uses a local image identified by its immutable registry digest.
`latest` is never used by the runtime or deployment scripts. Independent spot
checks may use the public service at an interactive rate and must record API
failures separately from numeric differences.

Install or update the local service on the Ubuntu host with:

```bash
sudo bash /home/ubuntu/ariane/scripts/server-ops/install-spliceai-service.sh
```

The installer reads the image from the active profile, pulls that exact digest,
starts it temporarily on a private port and runs
`scripts/validate_spliceai_service.py`. It changes the ARIANE environment only
after the response matches the versioned validation case. The permanent service
binds only to `127.0.0.1:8082` and does not use the optional server-side database.

An update starts with a candidate digest. It must pass the validation case,
classification regression suite and a representative BRCA comparison before the
profile is changed. Changing the profile checksum prevents reuse of results from
the previous engine.

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

The active image was built on 2026-09-03 and declares SpliceAI commit
`7f36ca847e1b1885167dab79681dbb75c09c6743`. It returns three-decimal scores
that match the versioned Broad API reference cases. The previous image from
2026-07-03 used an unpinned checkout of the model and returned older
two-decimal scores. It is not an approved runtime source.

Example command from the repository root:

```powershell
python variant_space_scan\validate_spliceai_reference_pilot.py --cache data\spliceai\spliceai_brca_snv_reference_cache.json --report variant_space_scan\outputs\spliceai_app_cache_vs_broad_reference_validation.YYYYMMDD.json --sample-size 100 --sleep-seconds 5 --retries 2 --retry-sleep-seconds 20 --timeout 90
```

The same JSON report can be converted to CSV with:

```powershell
python variant_space_scan\export_spliceai_validation_csv.py variant_space_scan\outputs\spliceai_app_cache_vs_broad_reference_validation.YYYYMMDD.json variant_space_scan\outputs\spliceai_app_cache_vs_broad_reference_validation.YYYYMMDD.csv
```
