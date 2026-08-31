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
`.runtime-cache/` directory. The mutable files are
`coordinates_api_cache.json`, `bayesdel_api_cache.json` and
`spliceai_api_cache.json`. They are read before a network request and remain
available when an upstream API is temporarily unavailable.

## Project structure

This overview lists the main runtime boundaries and entry points. Tests,
scripts and individual data files are omitted here.

```
ariane/
├── backend/
│   ├── main.py                         # FastAPI application and routes
│   ├── config.py                       # paths and runtime settings
│   ├── models.py                       # API request and response models
│   ├── classification_dag/
│   │   ├── domain.py                   # typed classification records
│   │   ├── engine.py                   # DAG validation and execution
│   │   ├── policy.py                   # ENIGMA point combination policy
│   │   ├── providers.py                # evidence provider nodes
│   │   ├── provider_wiring.py          # production provider adapters
│   │   ├── runtime.py                  # production graph assembly
│   │   └── nodes/                      # criterion-family rule nodes
│   ├── services/
│   │   ├── evidence_orchestration.py   # evidence lookup coordination
│   │   ├── variant_classification_service.py
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
│   ├── modules/
│   │   ├── pvs1.py                     # PVS1 and PM5 evaluation
│   │   ├── table4.py                   # ENIGMA Table 4 lookup
│   │   ├── table9.py                   # ENIGMA Table 9 lookup
│   │   ├── bp1.py                      # BP1 evaluation
│   │   ├── pp3_bp4.py                  # PP3 and BP4 evaluation
│   │   ├── bp7.py                      # BP7 evaluation
│   │   ├── manual_evidence.py          # validated manual evidence forms
│   │   └── external.py                 # external comparison formatting
│   ├── lookups/
│   │   ├── spliceai.py                 # SpliceAI API lookup
│   │   ├── bayesdel.py                 # BayesDel lookup
│   │   ├── clinvar.py                  # ClinVar comparison lookup
│   │   ├── clingen.py                  # ClinGen ERepo comparison lookup
│   │   └── coordinates.py              # HGVS to GRCh37/38 resolution
│   └── data/                            # immutable reference datasets
├── frontend/
│   ├── index.html
│   └── static/
│       ├── css/
│       │   └── style.css
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
2. **Table 9** - calibrated PS3/BS3 functional evidence
3. **Table 4** - PVS1/PM5 structural rules
4. **gnomAD** - BS1, PM2
5. **SpliceAI/BayesDel** - PP3/BP4/BP7 according to the ENIGMA variant-type decision tree; PP3 is not stacked with PVS1
6. **BP1** - outside functional domain
7. **RNA evidence review recommendation** - informational only, no scoring
8. **External comparison** - ClinVar + ClinGen ERepo, read-only

## Data sources

- ENIGMA VCEP v1.2 (2024-11-18): Table 4, Table 9
- gnomAD v2.1.1 exomes non-cancer
- gnomAD v3.1.2 genomes non-cancer
- SpliceAI: Broad API (spliceai-38-xwkwwwxdwq-uc.a.run.app)
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

This tool is a research prototype. Do not use for clinical decisions without expert review.

## Automated scope and manual review

The automated score covers evidence that can be evaluated from the bundled
ENIGMA tables and local lookup data. Evidence requiring expert interpretation
or manual confirmation is intentionally excluded from Module 1.

For exon duplications, the form accepts a laboratory-supplied duplication
arrangement. The default is `Unknown`. Select `Confirmed tandem` only when the
laboratory data support tandem arrangement; the application never infers it
from copy number alone.

Case-control, Fanconi anemia, family co-segregation, curated RNA evidence,
curated functional evidence outside Table 9, curated initiation-codon PVS1
evidence, and curated splice PS1 evidence (`PS3`, `PS4`, `PM3`, `PP1`, `BS2`,
`BS3`, `BS4`, `PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`, `PS1_SPLICE`) are not part of
the automatic Module 1 result. PP4 and BP5 outside the validated clinical LR
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

Install development dependencies and run the offline regression suite without
network access:

```bash
pip install -r requirements-dev.txt
python -m pytest tests -q
```

The VUS explanation layer and regression golden cases are documented in
`docs/vus_explanation_and_golden_cases.md`.
