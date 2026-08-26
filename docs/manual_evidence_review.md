# Manual Evidence Review

## Purpose

ARIANE Module 1 produces an automatic result from automatable ACMG/AMP and
ENIGMA BRCA1/2 VCEP criteria. Evidence requiring review of patients, families,
study design, or literature is intentionally excluded from that result.

The manual-review panel supports `PS4`, `PM3`, `PP1`, `BS2`, `BS4`,
`PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`, `PS1_SPLICE`, and `PS1_PROTEIN`. It creates a separate
amended working result and never replaces the original Module 1 classification.

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

`PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`, `PS1_SPLICE`, and `PS1_PROTEIN` require structured
curated records. Their allowed strength is derived only after all required
supporting fields have been validated.

The BP7 RNA form does not contain a general eligibility checkbox. The server
derives the variant type from the classified c. and p. notation, checks the
complete affected protein interval against the ENIGMA domains and reads BS3
from the original automated result. For a domain missense variant, only an
applied BS3 with ENIGMA Table 9 provenance satisfies the prerequisite. Missing
variant context, an unresolved protein position or missing BS3 fails closed and
adds no BP7 RNA points.

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
whose registry status is `review_required`. A reference marked `excluded` is
shown with its exclusion reason but cannot be manually confirmed as protein PS1.
It requires a qualifying VCEP classification verification, the same normalized
missense substitution, a different nucleotide change, SpliceAI at most 0.1 for
both variants, and a completed check of named RNA/splice sources for both
variants. The strength is derived from the reference class and cannot be freely
overridden. If the reference classification is known to use PS1, the reviewer
must identify that dependency and exclude a direct reciprocal dependency.

The protein PS1 form also accepts a reference c. HGVS description. The backend
normalizes that reference against the configured transcript, derives and verifies
its canonical p. consequence, compares it with the assessed variant, and obtains
SpliceAI for both variants through the configured profile. It also checks the exact
reference variant in ClinVar and ClinGen ERepo. These facts are filled into the
form, but they do not by themselves add PS1 points.

ClinVar review stars are used only to describe the candidate source. A two-star
aggregate assertion is not an ENIGMA VCEP assertion and therefore remains a
manual-review candidate. A three-star record qualifies as VCEP verification only
when the underlying assertion is from the applicable ENIGMA/ClinGen expert panel.
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
The server does not currently persist these records.

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
