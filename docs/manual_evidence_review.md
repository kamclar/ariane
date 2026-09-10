# Manual Evidence Review

## Purpose

ARIANE Module 1 produces an automatic result from automatable ACMG/AMP and
ENIGMA BRCA1/2 VCEP criteria. Evidence requiring review of patients, families,
study design, or literature is intentionally excluded from that result.

The manual-review panel supports `PS3`, `PS4`, `PM3`, `PP1`, `PP4`, `BS2`,
`BS3`, `BS4`, `BP5`, `PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`, `PS1_SPLICE`, and
`PS1_PROTEIN`. It creates a separate amended working result and never replaces
the original Module 1 classification.

## Form navigation and assistance

The interface groups the forms by evidence family: RNA and splicing, prior
variant evidence, functional evidence, clinical and case-control evidence,
family and segregation evidence, and co-occurrence evidence. Groups are used
only for navigation. They do not alter criterion eligibility or strength.

When the classification backend returns an explicit RNA, initiation-codon,
splice PS1, or protein PS1 review recommendation, the corresponding form is
shown under `Recommended reviews for this variant`. The recommendation opens
the existing expert-review form and never assigns the criterion.

ARIANE prefills objective facts that are already available from the classified
variant and pinned evidence sources. These include variant identity, transcript,
protein consequence, variant type, available SpliceAI context, exact Table 4,
Table 9, ST2 or ST7 records, and recorded source references where applicable.
An automatically filled field remains part of the submitted audit record.

Expert conclusions are not filled automatically. The reviewer must still
confirm matters such as cohort independence, VCEP applicability, assay
calibration, the same splice event, prediction comparability, and the relevance
of the tested tissue or biological system.

Each selected form displays one of three backend-derived states:

- `Not started`: the criterion has not been selected;
- `Needs information`: required evidence, notes, references, threshold data, or
  an ENIGMA stipulation is missing;
- `Ready`: the submitted fields meet the backend completeness rules and the
  backend can derive a criterion strength.

The browser requests these states from `/api/manual-evidence/status`. It does
not maintain a second list of required fields or calculate criterion strength.
The final amended result is still calculated separately by
`/api/manual-evidence/evaluate`.

## Criterion Thresholds

| Criterion | ARIANE strength suggestion |
| --- | --- |
| `PS4` | Strong when p-value is at most 0.05, odds ratio is at least 4, and the lower confidence limit is greater than 2 |
| `PM3` | Supporting at 1 evidence point, Moderate at 2-3 points, Strong at 4 or more points |
| `PP1` | Supporting at LR 2.08, Moderate at LR 4.3, Strong at LR 18.7, Very Strong at LR 350 |
| `BS2` | Supporting at 1 evidence point, Moderate at 2-3 points, Strong at 4 or more points |
| `BS4` | Supporting at LR 0.48, Moderate at LR 0.23, Strong at LR 0.05, Very Strong at LR 0.00285 |
| `PVS1_RNA` | Reviewer-curated mRNA-only damaging transcript record; Supporting, Moderate, Strong, or Very Strong |
| `BP7_RNA` | Reviewer-curated mRNA-only no-damaging-effect record; Strong only after variant type, domain and required BS3 checks |
| `PVS1_INIT` | Reviewer-curated initiation-codon PVS1 flowchart record; Supporting, Moderate, Strong, or Very Strong |
| `PS1_SPLICE` | Reviewer-curated same-splicing-impact PS1 record; Supporting, Moderate, or Strong |
| `PS1_PROTEIN` | Strong for a verified Pathogenic reference or Moderate for a verified Likely Pathogenic reference after complete protein-PS1 splice checks |

For `PP1`, the thresholds mean greater than or equal to the listed LR. For
`BS4`, they mean less than or equal to the listed LR.

A reported combined LR is sufficient to derive the BS4 strength. BS4 Strong
does not by itself satisfy the ENIGMA Table 3 route to Likely Benign unless at
least two independent LR components are documented. The BS4 form therefore
accepts structured component records containing the LR, source and a unique
independence-group identifier. The backend multiplies the component LRs and
requires the product to match any separately reported combined LR. One LR can
still produce BS4 Strong, but that single criterion alone remains insufficient
for Likely Benign. BS4 Very Strong alone is sufficient for Likely Benign under
the ENIGMA classification text.

The numeric threshold is not sufficient when ENIGMA defines an additional
stipulation:

- `PS4` requires country and ethnicity matching between case and control datasets.
- `PM3` requires a co-occurring P/LP variant classified using VCEP specifications
  and confirmation that the assessed variant does not meet benign population evidence.
- `BS2` requires a co-occurring P/LP variant classified using VCEP specifications.
- `PP1 Very Strong` requires a predicted or experimentally proven effect on
  protein or mRNA splicing. Without this evidence an LR of 350 or more is capped
  at `PP1 Strong`.
- `BP7_RNA` is permitted for intronic and synonymous variants and for
  missense or in-frame variants outside the ENIGMA clinically important
  functional domains. A missense variant inside such a domain must already
  meet BS3. An in-frame variant inside such a domain is not eligible for BP7
  Strong (RNA) under the stated BP7 rule.

The reviewer enters the evidence values and supporting provenance. ARIANE
derives the criterion strength from the ENIGMA v1.2 rule. A generic manual
strength override is not permitted. Evidence that does not meet a rule
threshold receives no criterion and no points.

Criterion eligibility, evidence completeness, strength, points, evidence
interactions and the amended classification are calculated only by the backend
manual-evidence DAG. The browser does not contain threshold tables, a second
strength calculator or evidence-completeness rules. It submits the raw form.
The backend selects thresholds and applicable rules from the checksum-validated
`backend/data/gene_policy_manifest.json` entry for the assessed gene.
The backend requires at least one enabled criterion, assessor, date, evidence
notes and at least one reference for every enabled criterion. Before submission
the browser states that the strength will be calculated by the backend. After
submission it displays `suggested_strength`, `selected_strength`, points and
validation errors returned by the backend.

`PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`, `PS1_SPLICE`, `PS1_PROTEIN`, `PS3`, `BS3`,
and `BP5` require structured curated records. Their allowed strength is derived
only after all required supporting fields have been validated.

The BP7 RNA form does not contain a general eligibility checkbox. The server
derives the variant type from the classified c. and p. notation, checks the
complete affected protein interval against the ENIGMA domains and reads BS3
from either the automated result or the submitted manual evidence. For a domain
missense variant, the prerequisite is satisfied by applied BS3 with ENIGMA
Table 9 provenance or by a complete manual BS3 record whose VCEP calibration,
controls, variant-specific result, source and reviewer pass backend validation.
Missing variant context, an unresolved protein position, a standalone checkbox
or incomplete BS3 fails closed and adds no BP7 RNA points.

For Met1/start-loss variants, ARIANE can show an initiation-codon review
recommendation and prefill `met1_loss_confirmed` in the `PVS1_INIT` manual
record. The reviewer must still complete the alternative start, upstream P/LP
evidence, functional-impact, strength, notes, and references fields.

`PS1_SPLICE` is separate from automated protein-level `PS1`. It requires a
known Pathogenic or Likely Pathogenic reference variant, the reference
classification source, confirmation that the variant under assessment has the
same splice event, similar or stronger prediction evidence, and a manual
Appendix J/Table 17 strength decision. ARIANE does not provide an active
splice-PS1 reference registry and does not prefill a criterion strength.
The interface can prefill factual source fields from P/LP records selected
directly from the complete official ST2 snapshot. This is candidate discovery
only. It does not confirm reference eligibility, same-event matching or the
Appendix J/Table 17 branch. The reviewer must document the complete eligibility
assessment from the primary source.

ST2 alone cannot establish splice PS1. ENIGMA additionally requires all of the
following checks:

- the reference P/LP classification was assigned using VCEP specifications;
- the assessed and reference variants produce precisely the same splice event;
- the splice prediction for the assessed variant is similar to or stronger than
  the reference prediction;
- the applicable Appendix J/Table 17 branch is selected from the positions of
  both variants within the donor or acceptor motif;
- the assessed variant's baseline PP3 or PVS1 result is included in the Table 17
  decision;
- any concurrent protein-level consequence is reviewed for an exonic variant.

The ST2 prefill therefore leaves same-event confirmation, prediction comparison
and criterion strength unset.

`PS1_PROTEIN` is prefilled when ARIANE finds a matching P/LP missense reference
but a runtime or splice condition remains unresolved. A reference marked `excluded` is
shown with its exclusion reason but cannot be manually confirmed as protein PS1.
The prefill includes the normalized reference variant and protein consequence,
the ST7 class and source, objective variant comparison, available SpliceAI
results, and recorded checks of the defined RNA/splice sources. ARIANE then
checks ClinVar and ClinGen ERepo for a current ENIGMA VCEP assertion in the
background. An ST7 P/LP record is an accepted classification basis. Fields that
cannot be established remain explicitly unresolved. It requires the same normalized
missense substitution, a different nucleotide change, SpliceAI at most 0.1 for
both variants, and a completed check of named RNA/splice sources for both
variants. The strength is derived from the reference class and cannot be freely
overridden. The reviewer must establish whether the reference classification
used PS1. An unresolved answer does not add points. If PS1 was used, the reviewer
must identify that dependency and exclude a direct reciprocal dependency.

The protein PS1 form also accepts a reference c. HGVS description. The backend
normalizes that reference against the configured transcript, derives and verifies
its canonical p. consequence, compares it with the assessed variant, and obtains
SpliceAI for both variants through the configured profile. It also checks the exact
reference variant in ClinVar and ClinGen ERepo. These facts are filled into the
form, but they do not by themselves add PS1 points.

ClinVar review stars are used only to describe review status. No star count
qualifies a reference for PS1. A separately identified assertion from the
applicable ENIGMA/ClinGen VCEP may qualify independently of the aggregate star
count. An ordinary ClinVar aggregate conclusion does not prefill the reference
classification, verification, classification source, or evidence references.
The named RNA/splice source review and reciprocal-dependency check remain required.
If normalization, SpliceAI, ClinVar, or ClinGen is unavailable, the corresponding
field is reported as unavailable and no missing value is interpreted as evidence.

## Audit Record

Each enabled criterion requires:

- evidence values used for the threshold calculation
- a reviewer note describing the evidence and limitations
- at least one PMID, DOI, URL, or internal evidence record
- assessor name or identifier
- assessment date

The browser can export the original Module 1 result, submitted evidence,
the backend-derived strengths, and amended working result as JSON. Raw submitted
manual evidence does not contain a separate frontend-derived strength.

An authenticated reviewer can also save the review as an immutable server-side
draft. ARIANE recomputes the amended result before saving and does not trust a
classification submitted by the browser. The record contains the complete
Module 1 result, enabled manual evidence, amended result, reviewer identity and
role, application and policy versions, and SHA-256 hashes of all three content
sections.

Approval creates a new immutable record version that references the saved draft.
The draft is not updated or deleted. New evidence must be evaluated and saved as
a new draft. Records are stored in `review_records.sqlite3` below the directory
configured by `ARIANE_RUNTIME_DATA_DIR`. The Railway default is
`${RAILWAY_VOLUME_MOUNT_PATH}/ariane-runtime-data`; local development uses
`.runtime-data/`. This storage is separate from replaceable API caches and is
excluded from Git.

The current beta uses the existing protected administrator account as the
authenticated account and separately records the stated reviewer name and
professional role. This provides authenticated, versioned beta records but is
not a substitute for institutional user management, role assignment or an
electronic-signature system. A production laboratory deployment should connect
the same record API to its identity provider and authorization policy.

Evidence notes must not contain direct patient identifiers unless the deployment
has an approved protected-data environment and an applicable retention policy.

## Development-only amended-result walkthrough

The following input tests the form workflow and must not be treated as clinical
evidence or saved as an approved review:

1. Classify `BRCA1 c.5366C>T p.(Ala1789Val)`.
2. Confirm that the Module 1 result is `Class 3`, 5 points, with `PS3 Strong`
   and `PP3 Supporting`. Stop the walkthrough if the current source bundle gives
   a different starting result.
3. Open `Family and segregation evidence` and select `PP1`.
4. Enter LR `2.08`, notes `Synthetic UI test only. No clinical evidence
   asserted.`, and reference `TEST:manual-review-demo`.
5. Enter assessor `Test reviewer` and role `Test only`.
6. Calculate the amended result. The backend should derive `PP1 Supporting`,
   add 1 point, and return `Class 4, Likely Pathogenic`, 6 points.

Do not use `Save immutable review draft` or `Approve saved version` for this
synthetic walkthrough. Persistent records are intended only for real evidence
reviewed under the applicable local governance process.

## ClinVar Display

ARIANE displays official ClinVar review stars for the aggregate assertion:

| Stars | ClinVar review level |
| ---: | --- |
| 0 | no assertion criteria or no classification |
| 1 | criteria provided by one submitter or conflicting submitters |
| 2 | criteria provided by multiple submitters with no conflicts |
| 3 | reviewed by an expert panel |
| 4 | practice guideline |

These stars describe the review status of the assertion. They are not a general
quality score for a laboratory.

Individual ENIGMA submissions are marked as
`ClinGen/ENIGMA curated submitter`. ARIANE does not currently assign a custom
credibility score to other laboratories. Any future list should have explicit,
documented inclusion criteria and versioning.

## Primary Sources

- [ACMG/AMP sequence variant interpretation guidelines](https://pubmed.ncbi.nlm.nih.gov/25741868/)
- [Tavtigian et al. point-based classification framework](https://pubmed.ncbi.nlm.nih.gov/32720330/)
- [ENIGMA BRCA1/2 VCEP v1.2 criteria registry](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2.0)
- [Specifications v1.2](https://cspec.genome.network/cspec/File/id/11e62fec-23b0-4a3e-b2df-751855301746/data)
- [Appendix v1.2](https://cspec.genome.network/cspec/File/id/5a75d1a0-1222-46a2-8802-68a4f2251a3a/data)
- [Supplementary tables v1.2](https://cspec.genome.network/cspec/File/id/3dadda2f-94a3-497f-aa35-3bb6e828ddd5/data)
- [Specifications Table 4](https://cspec.genome.network/cspec/File/id/10301df8-45e0-4309-adba-c121eb057d3e/data)
- [Specifications Table 9](https://cspec.genome.network/cspec/File/id/c540f11d-0be2-45d6-a0bf-ae5327a04885/data)
- [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
