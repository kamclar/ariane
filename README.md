# ARIANE

**Automated ACMG Rule-based Interpretation and Annotation ENgine**

BRCA1/2 variant classification following ENIGMA VCEP v1.2.

## Quick start (local)

```bash
# create virtual environment
python -m venv venv
source venv/bin/activate   # or venv\Scripts\activate on Windows

# install dependencies
pip install -r requirements.txt

# required classification datasets and the checksum-verified HGVS reference
# bundle are versioned in this repository

# run
uvicorn backend.main:app --reload --port 8000

# open http://localhost:8000
```

The supported production runtime is Linux with Python 3.12. On Ubuntu install
`python3-dev`, `build-essential`, and `libpq-dev` before Python dependencies.
The biocommons dependency stack does not currently provide all required wheels
for native Windows; use WSL for the complete local runtime on Windows.

## Deploy to Railway

```bash
# push to GitHub, then connect repo in Railway dashboard
# or use Railway CLI:
railway up
```

Attach a Railway Volume to the service for persistent runtime lookup caches.
The application automatically uses `RAILWAY_VOLUME_MOUNT_PATH` and stores its
mutable files below `ariane-runtime-cache/`. Outside Railway, set
`ARIANE_RUNTIME_CACHE_DIR` to a writable persistent directory. Precomputed
snapshots remain in the repository and are never modified at runtime.
Without either deployment setting, local development uses the ignored
`.runtime-cache/` directory. The mutable files are `bayesdel_api_cache.json`
and `spliceai_api_cache.json`. They are read before a network request and
remain available when an upstream API is temporarily unavailable. Genomic
coordinates come only from checksum-validated local sources registered in
`data/coordinates/coordinate_sources.manifest.json`.

## Project structure

This overview lists the main runtime boundaries and entry points. Tests,
scripts and individual data files are omitted here.

```
ariane/
├── backend/
│   ├── main.py                         # FastAPI composition root
│   ├── bootstrap.py                    # explicit startup validation and runtime wiring
│   ├── api/
│   │   ├── classification.py           # classification HTTP transport and quotas
│   │   ├── manual.py                   # manual evidence and normalization routes
│   │   ├── system.py                   # health, rules, resources and maintenance
│   │   └── auth.py / session.py        # API key and browser-session boundaries
│   ├── contracts/                      # typed classification, batch and review DTOs
│   ├── domain/                         # shared classification records and audit structures
│   ├── policy/                         # gene, combination and SpliceAI policies
│   ├── infrastructure/                 # runtime paths, health, caches and repositories
│   ├── classification_dag/
│   │   ├── engine.py                   # DAG validation and execution
│   │   ├── providers.py                # evidence provider nodes
│   │   ├── provider_wiring.py          # production provider adapters
│   │   ├── runtime.py                  # production graph assembly
│   │   └── nodes/                      # criterion-family rule nodes
│   ├── services/
│   │   ├── evidence_orchestration.py   # evidence lookup coordination
│   │   ├── variant_classification_service.py # single/batch workflow, cache and usage
│   │   ├── classification_presentation.py
│   │   └── ps1_reference_resolution.py
│   ├── population_frequency/
│   │   ├── service.py                  # gnomAD lookup service
│   │   ├── snapshot_repository.py      # validated snapshot loading
│   │   ├── criteria.py                 # BA1, BS1 and PM2 decisions
│   │   ├── coverage.py                 # coverage checks
│   │   ├── indel_size.py               # Appendix G indel sizing
│   │   ├── lookup.py                   # frequency record lookup
│   │   ├── policy.py                   # dataset and policy bindings
│   │   └── models.py                   # population evidence records
│   ├── criteria/                        # criterion evaluators and interactions
│   ├── reference_data/                  # validated ENIGMA tables and registries
│   ├── review/
│   │   ├── definitions.py              # manual form definitions and sources
│   │   ├── strength.py                 # policy-bound strength derivation
│   │   ├── validation.py               # form readiness and completeness
│   │   ├── service.py                  # amended working classification
│   │   └── *_review.py                 # variant-specific review builders
│   ├── presentation/                    # narratives, external comparison and display ordering
│   ├── variant_processing/              # HGVS normalization and variant typing
│   ├── lookups/
│   │   ├── spliceai.py                 # SpliceAI API lookup
│   │   ├── bayesdel.py                 # BayesDel lookup
│   │   ├── clinvar.py                  # ClinVar comparison lookup
│   │   ├── clingen.py                  # ClinGen ERepo comparison lookup
│   │   └── coordinates.py              # validated local GRCh37/38 resolution
│   └── data/                            # immutable reference datasets
├── frontend/
│   ├── index.html                    # page shell and explicit template includes
│   ├── templates/                    # feature-level Alpine templates
│   └── static/
│       ├── css/                      # base, forms, results, evidence, layout and modes
│       └── js/
│           ├── app.js                  # Alpine application assembly
│           ├── api.js                  # backend API calls
│           ├── batch.js                # batch input and results
│           ├── classification.js       # classification result state
│           ├── core.js                 # shared frontend state
│           ├── formatters.js           # display formatting
│           ├── graphs.js               # decision graph rendering
│           ├── manual-review.js        # manual evidence forms
│           └── rules.js                # rule explorer state
├── requirements.txt
├── railway.toml
└── README.md
```

## Evidence hierarchy

Classification follows this order (higher level overrides lower):

1. **BA1** - stand-alone benign (gnomAD non-cancer FAF95 > 0.1%), pouze po
   kontrole pokrytí, QC filtru a výjimky pro patogenní founder varianty
2. **Approved ERepo PVS1 RNA registry** - exact ENIGMA v1.2 RNA assertions
3. **Table 9** - calibrated PS3/BS3 functional evidence
4. **Table 4** - PVS1/PM5 structural rules
5. **gnomAD** - BS1, PM2
6. **SpliceAI/BayesDel** - PP3/BP4/BP7 according to the ENIGMA variant-type decision tree; PP3 is not stacked with PVS1
7. **BP1** - outside functional domain
8. **RNA evidence review recommendation** - for evidence not covered by an approved exact registry record
9. **External comparison** - live ClinVar + ClinGen ERepo, read-only

## Data sources

- ENIGMA VCEP v1.2 (2024-11-18): Table 4, Table 9
- ClinGen ERepo ENIGMA BRCA1/2 VCEP v1.2 PVS1 RNA assertions in a local
  checksum-validated registry
- gnomAD v2.1.1 exomes non-cancer
- gnomAD v3.1.2 genomes non-cancer
- SpliceAI: local Broad-compatible service using the digest-pinned image from
  `data/spliceai/enigma_v1_2_spliceai_profile.json`
- BayesDel: myvariant.info
- ClinVar: NCBI eutils

The gnomAD releases, official Hail Table identities and panel intervals are
pinned in `backend/data/gnomad/gnomad_panel_manifest.json`. Check for published
releases without changing the active classification data:

```bash
python scripts/refresh_gnomad_panel_snapshot.py check-updates
```

Refresh and validate the panel snapshots in a separate data-build environment:

```bash
pip install -r requirements-data.txt
python scripts/refresh_gnomad_panel_snapshot.py refresh
python scripts/refresh_gnomad_panel_snapshot.py validate
```

New releases are never activated automatically. The manifest is gene-extensible,
but an interval alone cannot activate a gene. Every target must reference an
explicit active gene-specific policy containing its VCEP provenance, transcript,
frequency datasets, population groups, thresholds, coverage rules and excluded
variant types. The BRCA policy is never inherited by another gene.
Runtime gene and VCEP configuration is held in
`backend/data/gene_policy_manifest.json` with checksum metadata. It is the
authoritative source for active genes, reference transcripts, VCEP policy IDs,
decision thresholds, functional domains and applicable rules. Source-specific
manifests are checked against it at startup and cannot silently override the
policy.

Each VCEP policy also declares an `implementation_profile`. The runtime fails
closed when that profile has no registered DAG implementation. Input gene
prefixes, VCEP links, HGVS startup checks, PVS1 decision assets and domain
descriptions are manifest-driven. BRCA-specific source names remain explicit
only for datasets whose documented scope is BRCA1/2.

After an approved policy or threshold change, increment `manifest_version` and
refresh the checksum metadata with:

```powershell
.\venv\Scripts\python.exe scripts\update_gene_policy_manifest_metadata.py --write
```

The update check also verifies that a newer release contains the equivalent
small-variant Hail Table. A release directory for another data type is not
reported as a usable frequency-data update.

Population scoring follows ENIGMA Appendix G. Only AFR, AMR, EAS, NFE and SAS
contribute to BA1/BS1 and outbred-population presence for PM2. Founder groups
ASJ, FIN and AMI, plus other non-scoring groups present in a release, are stored
and displayed as context but cannot change a criterion. Well-established
pathogenic founder variants are checked separately and cannot receive BA1/BS1.
- ClinGen: Evidence Repository API
- RNA evidence review recommendation: ARIANE review aid for RNA-dependent or
  predicted splice-effect situations; not an ACMG/AMP or ENIGMA criterion and
  not included in scoring
- Splice PS1 review recommendation: ARIANE review aid for possible
  same-splicing-impact PS1 scenarios; not scored automatically

## Disclaimer

ARIANE is under expert validation. Do not use it for clinical decisions without expert review.

## Automated scope and manual review

The automated score covers evidence that can be evaluated from the bundled
ENIGMA tables and local lookup data. Evidence requiring expert interpretation
or manual confirmation is intentionally excluded from Module 1.

For exon duplications, the form accepts a laboratory-supplied duplication
arrangement. The default is `Unknown`. Select `Confirmed tandem` only when the
laboratory data support tandem arrangement; the application never infers it
from copy number alone.

Case-control, Fanconi anemia, family co-segregation, RNA evidence outside an
exact approved ERepo assertion or an eligible curated ST2 qualitative branch,
curated functional evidence outside Table 9, curated initiation-codon PVS1
evidence, and curated splice PS1 evidence (`PS3`, `PS4`, `PM3`, `PP1`, `BS2`,
`BS3`, `BS4`, `PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`, `PS1_SPLICE`) are handled in
manual review. Exact published PVS1 RNA assertions in the versioned local ERepo
registry and unambiguous ENIGMA-curated ST2 patient-mRNA records can enter the
automatic Module 1 result. PP4 and BP5 outside the validated clinical LR
snapshot can also be reviewed manually. After a variant is classified, the
user can enter these evidence types in a separate manual-review panel. ARIANE
derives the permitted strength in the backend from ENIGMA VCEP v1.2 thresholds
or requires a complete documented VCEP calibration review. A manual strength
override is not accepted. BRCA1/2 VCEP v1.2 permits only Strong for PS3 and BS3.

RNA evidence is accepted only as a structured, reviewer-curated mRNA-only assay
record. `PVS1_RNA` supports PVS1 (RNA) at Supporting, Moderate, Strong, or Very
Strong strength. `BP7_RNA` supports BP7_Strong (RNA). Protein-only or combined
mRNA/protein assays should be reviewed as PS3/BS3 evidence, not as RNA-only
PVS1/BP7 evidence.

Initiation-codon variants are recognized automatically, but they are not scored
from Met1 loss alone. `PVS1_INIT` supports a structured reviewer-curated
initiation-codon PVS1 flowchart record with Supporting, Moderate, Strong, or
Very Strong strength.
ARIANE shows an initiation-codon PVS1 review recommendation for Met1/start-loss
variants and pre-fills the safest manual-review fields, but it does not infer
the downstream alternative start analysis or award automatic PVS1 points.
ARIANE shows an initiation-codon PVS1 review recommendation for Met1/start-loss
variants and pre-fills the safest manual-review fields, but it does not infer
the downstream alternative start analysis or award automatic PVS1 points.

The application always preserves the original Module 1 result and displays a
separate amended working result. The reviewer must provide their identifier,
assessment date, evidence notes, and references. The complete audit record can
be exported as JSON. It is not stored automatically by the server.

For RNA-dependent or predicted splice-effect situations, ARIANE may show an
RNA evidence review recommendation. This is a review aid only. It does not add
criteria, points, or change the Module 1 classification.

For splice-relevant variants, ARIANE may also show a Splice PS1 review
candidate notice. This is separate from Table 4/PVS1 and from the automated
protein-level PS1 implementation. It indicates that a reviewer should look for
a known P/LP reference variant with the same documented or confidently
predicted splice consequence.

### API status and batch responses

The versioned public API is available under `/api/v1`. It is a beta contract for
external integrations and uses the same production classification DAG as the
web interface. See [docs/public_api.md](docs/public_api.md) for requests,
responses, limits and a reference client.

The classification routes under `/api/v1` require an individual API key.
Only its SHA-256 digest is stored on the server. The public capabilities route
does not require a key. Unversioned `/api` classification and data routes also
require a key. The browser uses separate `/ui-api` routes with a short-lived,
server-signed HttpOnly session established when the page loads. No reusable API
key is present in the HTML or JavaScript. Interactive classification permits 30
requests per minute per client IP with a burst of 3.

Local development also requires a session signing secret before starting
Uvicorn. Use a development-only value of at least 32 bytes, for example
`ARIANE_UI_SESSION_SECRET=local-development-secret-change-me`. Server
installation scripts generate a cryptographically random value instead.

`POST /api/v1/classify/batch` accepts at most 10 items and preserves their input
order. Five items are recommended for uncached work. Each item is validated
separately. An invalid variant is returned as an item with `status: "error"`;
it does not prevent valid items in the same batch from being classified. HTTP
422 is reserved for an invalid top-level request, an empty list, or a batch
exceeding the item limit. Clients should not retry a 422 response without
correcting the request.

The v1 response places review aids in a separate `manual_review` object.
`affects_automatic_classification: false` means that a recommendation adds no
criterion, points, or automatic class change. A separate, documented
manual-evidence request is required to calculate an amended working result.

`PS1_SPLICE` can then be added manually as a structured curated record when a
reviewer confirms the reference variant, P/LP classification source, same splice
event, similar or stronger prediction evidence, and Appendix J/Table 17
strength. ARIANE does not infer this strength automatically.

The manual-review form can search factual P/LP splice candidates derived
directly from the complete official ENIGMA Supplementary Table 2 snapshot. A
selection prefills only source facts such as the reference variant, reported
splice event, assay context and multifactorial class. It does not confirm PS1
eligibility, same-event matching, prediction strength or criterion strength.

ClinVar review stars are displayed as the official review level of the
aggregate ClinVar assertion. Individual submitters are not assigned an
ARIANE-generated star score. ENIGMA submissions are identified separately as
`ClinGen/ENIGMA curated submitter`. Additional curated submitters should only be
added from an explicit, documented list.

See `docs/manual_evidence_review.md` for thresholds, sources, and limitations.

## Tests

Use Python 3.12. A virtual environment created with Python 3.13 or newer is not
supported by the pinned scientific dependencies. Install development
dependencies and run the offline regression suite without network access:

```bash
pip install -r requirements-dev.txt
python -m ruff check .
python -m mypy
python -m pytest tests -q
```

The same lint, type and test checks run automatically for every push and pull
request through `.github/workflows/quality.yml`. Tool configuration is kept in
`pyproject.toml`; exact development-tool versions are pinned in
`requirements-dev.txt`.

The VUS explanation layer and regression golden cases are documented in
`docs/vus_explanation_and_golden_cases.md`.
