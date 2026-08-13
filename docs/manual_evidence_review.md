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
| `BP7_RNA` | Reviewer-curated mRNA-only no-damaging-effect record; Strong only |
| `PVS1_INIT` | Reviewer-curated initiation-codon PVS1 flowchart record; Supporting, Moderate, Strong, or Very Strong |
| `PS1_SPLICE` | Reviewer-curated same-splicing-impact PS1 record; Supporting, Moderate, or Strong |
| `PS1_PROTEIN` | Strong for a verified Pathogenic reference or Moderate for a verified Likely Pathogenic reference after complete protein-PS1 splice checks |

For `PP1`, the thresholds mean greater than or equal to the listed LR. For
`BS4`, they mean less than or equal to the listed LR.

The reviewer can replace the ARIANE suggestion with another strength permitted
for non-RNA and non-initiation criteria. Both the suggestion and the selected
strength are retained in the audit output. A reviewer override does not
establish that the evidence meets the ENIGMA threshold.

`PVS1_RNA`, `BP7_RNA`, `PVS1_INIT`, `PS1_SPLICE`, and `PS1_PROTEIN` require structured
curated records and do not accept a free-text strength override without the
required supporting fields.

For Met1/start-loss variants, ARIANE can show an initiation-codon review
recommendation and prefill `met1_loss_confirmed` in the `PVS1_INIT` manual
record. The reviewer must still complete the alternative start, upstream P/LP
evidence, functional-impact, strength, notes, and references fields.

`PS1_SPLICE` is separate from automated protein-level `PS1`. It requires a
known Pathogenic or Likely Pathogenic reference variant, the reference
classification source, confirmation that the variant under assessment has the
same splice event, similar or stronger prediction evidence, and a manual
Appendix J/Table 17 strength decision. The pilot
`backend/data/splice_ps1_reference_set.json` can help identify candidates, but
it is unreviewed seed material and is not used for automatic scoring.

In the UI, the `PS1_SPLICE` candidate selector can prefill reference-variant
fields from the pilot set. The reviewer must still verify the match, choose the
Appendix J/Table 17 strength, and document the rationale.
When a candidate is selected, ARIANE may prefill a provisional strength from
the reference classification (`Pathogenic` -> Strong, `Likely Pathogenic` ->
Moderate). This is a convenience default only and must be confirmed or changed
by the reviewer.

`PS1_PROTEIN` is prefilled when ARIANE finds a matching P/LP missense reference
whose registry status is `review_required`. A reference marked `excluded` is
shown with its exclusion reason but cannot be manually confirmed as protein PS1.
It requires a qualifying VCEP classification verification, the same normalized
missense substitution, a different nucleotide change, SpliceAI at most 0.1 for
both variants, and a completed check of named RNA/splice sources for both
variants. The strength is derived from the reference class and cannot be freely
overridden. If the reference classification is known to use PS1, the reviewer
must identify that dependency and exclude a direct reciprocal dependency.

## Audit Record

Each enabled criterion requires:

- evidence values used for the threshold calculation
- a reviewer note describing the evidence and limitations
- at least one PMID, DOI, URL, or internal evidence record
- assessor name or identifier
- assessment date

The browser can export the original Module 1 result, submitted evidence,
ARIANE suggestions, reviewer overrides, and amended working result as JSON.
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
- [ENIGMA BRCA1/2 VCEP v1.2 criteria registry](https://cspec.genome.network/cspec/ui/svi/doc/GN097)
- [Specifications v1.2](https://cspec.genome.network/cspec/File/id/02537f62-66a3-497f-aa35-3bb6e828ddd5/data)
- [Appendix v1.2](https://cspec.genome.network/cspec/File/id/5a75d1a0-1222-46a2-8802-68a4f2251a3a/data)
- [Supplementary tables v1.2](https://cspec.genome.network/cspec/File/id/3dadda2f-94a3-497f-aa35-3bb6e828ddd5/data)
- [Specifications Table 4](https://cspec.genome.network/cspec/File/id/10301df8-45e0-4309-adba-c121eb057d3e/data)
- [Specifications Table 9](https://cspec.genome.network/cspec/File/id/c540f11d-0be2-45d6-a0bf-ae5327a04885/data)
- [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
