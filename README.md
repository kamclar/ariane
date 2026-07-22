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

# place data files in backend/data/
#   enigma_table4.json
#   enigma_table9.json
#   gnomad_brca_cache/ (directory with regional cache files)

# run
uvicorn backend.main:app --reload --port 8000

# open http://localhost:8000
```

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

## Project structure

```
ariane/
├── backend/
│   ├── main.py              # FastAPI app, routes
│   ├── config.py            # settings, thresholds
│   ├── models.py            # Pydantic request/response models
│   ├── modules/
│   │   ├── classifier.py    # main evaluation with evidence hierarchy
│   │   ├── pvs1.py          # PVS1/PM5 - Table 4 decision tree
│   │   ├── table4.py        # Table 4 loading and lookup functions
│   │   ├── table9.py        # Table 9 PS3/BS3 lookup
│   │   ├── bp1.py           # BP1 - outside functional domain
│   │   ├── pp3_bp4.py       # PP3/BP4 - BayesDel + SpliceAI
│   │   ├── bp7.py           # BP7 - synonymous without splice effect
│   │   ├── frequency.py     # BA1/BS1/PM2 - gnomAD frequencies
│   │   ├── external.py      # external comparison logic
│   │   └── utils.py         # shared helpers (AA position, domains)
│   ├── lookups/
│   │   ├── spliceai.py      # Broad SpliceAI API + Drive cache
│   │   ├── bayesdel.py      # myvariant.info BayesDel lookup
│   │   ├── clinvar.py       # ClinVar eutils VCV parser
│   │   ├── clingen.py       # ClinGen Evidence Repository API
│   │   └── coordinates.py   # HGVS → GRCh37/38 resolution
│   └── data/
│       ├── enigma_table4.json
│       ├── enigma_table9.json
│       └── spliceai_api_cache.json
├── frontend/
│   ├── index.html
│   └── static/
│       ├── css/style.css
│       └── js/app.js
├── requirements.txt
├── railway.toml
└── README.md
```

## Evidence hierarchy

Classification follows this order (higher level overrides lower):

1. **BA1** - stand-alone benign (gnomAD FAF > 0.1%) → class 1, stop
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
- SpliceAI: Broad API (spliceai-38-xwkwwwxdwq-uc.a.run.app)
- BayesDel: myvariant.info
- ClinVar: NCBI eutils
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
curated initiation-codon PVS1 evidence, and curated splice PS1 evidence (`PS4`,
`PM3`, `PP1`, `BS2`, `BS4`, `PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`,
`PS1_SPLICE`) are not part of the automatic Module 1 result. After a variant is
classified, the user can enter these evidence types in a separate manual-review
panel. ARIANE proposes the strength from the ENIGMA VCEP v1.2 thresholds, and
the reviewer may select another strength allowed for non-RNA, non-initiation,
and non-splice-PS1 structured criteria.

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

A pilot, unreviewed seed set for this review is available at
`backend/data/splice_ps1_reference_set.json`; see
`docs/splice_ps1_reference_set.md`. The manual-review UI can use it to prefill
reference fields and a provisional strength suggestion. It is not used for
automatic scoring.

ClinVar review stars are displayed as the official review level of the
aggregate ClinVar assertion. Individual submitters are not assigned an
ARIANE-generated star score. ENIGMA submissions are identified separately as
`ClinGen/ENIGMA curated submitter`. Additional curated submitters should only be
added from an explicit, documented list.

See `docs/manual_evidence_review.md` for thresholds, sources, and limitations.

## Tests

Run the offline regression suite without network access:

```bash
python -m unittest discover -s tests -v
```

The VUS explanation layer and regression golden cases are documented in
`docs/vus_explanation_and_golden_cases.md`.
