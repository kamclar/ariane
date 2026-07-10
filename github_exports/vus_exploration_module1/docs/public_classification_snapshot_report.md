# Public Classification Snapshot

Generated: 2026-06-22 14:21

## Purpose

This analysis adds a public assertion snapshot to the VUS review worklist. It
does not use public assertions as a validation target. ClinVar, ClinGen, and
ENIGMA-linked assertions are treated as public clinical/laboratory context that can help
prioritize manual review, identify conflicts, and avoid missing already
curated interpretations.

Queried variants: `80` selected from the highest-priority VUS action
plan.

## Snapshot Categories

| public_snapshot_category | count |
| --- | --- |
| multi-submitter assertion | 28 |
| single-submitter assertion | 28 |
| conflicting public assertions | 11 |
| no public assertion | 7 |
| panel-level public assertion | 4 |
| ClinGen/ENIGMA assertion | 2 |

## Discordance Worklist Labels

| discordance_label | count |
| --- | --- |
| module1_vus_public_vus | 40 |
| module1_vus_public_pathogenic_direction | 13 |
| public_conflict | 11 |
| no_public_classification | 8 |
| module1_vus_public_benign_direction | 7 |
| other_or_unclear | 1 |

`no public assertion` is a source category, while `no_public_classification` is
a worklist label. These counts can differ. In this snapshot, 7 variants had no
ClinVar record, but one additional variant had a ClinVar record with aggregate
classification `not provided`, giving 8 variants without a usable public
classification label.

## Lookup Status

ClinVar:

| clinvar_status | count |
| --- | --- |
| ok | 73 |
| not_found | 7 |

ClinGen Evidence Repository:

| clingen_status | count |
| --- | --- |
| not_found | 78 |
| ok | 2 |

## ClinGen ERepo BRCA1/2 Coverage Context

ClinGen Evidence Repository is useful when a BRCA1/2 record is present, but it
is not comprehensive for these genes. A manual check of the public ClinGen
ERepo BRCA1/2 table showed the following small record counts:

| Gene | Total ERepo records | Pathogenic | Likely pathogenic | VUS | Likely benign | Benign |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCA1 | 75 | 9 | 1 | 16 | 33 | 16 |
| BRCA2 | 68 | 12 | 6 | 14 | 26 | 10 |

Public ClinGen ERepo entry point: https://erepo.clinicalgenome.org/evrepo/

Therefore, `clingen_status=not_found` should be interpreted as "not found in
this small ERepo subset", not as absence of public interpretation. ClinVar
expert-panel or ENIGMA-linked assertions may still exist for the same variant.

## Category By Evidence Action Group

| public_snapshot_category | action_group | count |
| --- | --- | --- |
| single-submitter assertion | computational_signal_needs_noncomputational_support | 28 |
| multi-submitter assertion | computational_signal_needs_noncomputational_support | 26 |
| conflicting public assertions | computational_signal_needs_noncomputational_support | 9 |
| no public assertion | computational_signal_needs_noncomputational_support | 7 |
| ClinGen/ENIGMA assertion | computational_signal_needs_noncomputational_support | 2 |
| panel-level public assertion | computational_signal_needs_noncomputational_support | 2 |
| conflicting public assertions | manual_conflict_adjudication | 1 |
| conflicting public assertions | near_pathogenic_threshold | 1 |
| multi-submitter assertion | manual_conflict_adjudication | 1 |
| multi-submitter assertion | near_pathogenic_threshold | 1 |
| panel-level public assertion | manual_conflict_adjudication | 1 |
| panel-level public assertion | near_pathogenic_threshold | 1 |

## Highest Priority Public Context Rows

| priority_score | priority_tier | public_snapshot_category | discordance_label | gene | c_notation | p_notation | clinvar_aggregate_classification | clinvar_review_status | clinvar_n_submitters | preferred_submitters | clingen_classification |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 90 | tier1_urgent | multi-submitter assertion | module1_vus_public_pathogenic_direction | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | Pathogenic | criteria provided, multiple submitters, no conflicts | 15 | Laboratory for Molecular Medicine, Mass General Brigham Personalized Medicine: Pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Pathogenic; Ambry Genetics: Pathogenic; Quest Diagnostics Nichols Institute San Juan Capistrano: Pathogenic; GeneDx: Pathogenic; Baylor Genetics: Pathogenic; Color Diagnostics, LLC DBA Color Health: Pathogenic |  |
| 90 | tier1_urgent | panel-level public assertion | module1_vus_public_benign_direction | BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | Benign | reviewed by expert panel | 6 | Evidence-based Network for the Interpretation of Germline Mutant Alleles (ENIGMA): Benign; GeneDx: Likely benign; Labcorp Genetics (formerly Invitae), Labcorp: Likely benign; Ambry Genetics: Uncertain significance |  |
| 90 | tier1_urgent | panel-level public assertion | module1_vus_public_pathogenic_direction | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | Pathogenic | reviewed by expert panel | 39 | Baylor Genetics: Pathogenic; GeneDx: Pathogenic; Baylor Genetics: Pathogenic; Color Diagnostics, LLC DBA Color Health: Pathogenic; Laboratory for Molecular Medicine, Mass General Brigham Personalized Medicine: Pathogenic; Ambry Genetics: Pathogenic; Quest Diagnostics Nichols Institute San Juan Capistrano: pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Pathogenic |  |
| 90 | tier1_urgent | conflicting public assertions | public_conflict | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | Conflicting classifications of pathogenicity | criteria provided, conflicting classifications | 3 | Baylor Genetics: Uncertain significance; Ambry Genetics: Likely pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 90 | tier1_urgent | multi-submitter assertion | module1_vus_public_pathogenic_direction | BRCA1 | c.5453A&gt;G | p.(Asp1818Gly) | Pathogenic/Likely pathogenic | criteria provided, multiple submitters, no conflicts | 11 | Ambry Genetics: Pathogenic; Color Diagnostics, LLC DBA Color Health: Pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Pathogenic |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_vus | BRCA1 | c.299A&gt;G | p.(Glu100Gly) | Uncertain significance | criteria provided, multiple submitters, no conflicts | 2 | Baylor Genetics: Uncertain significance; Ambry Genetics: Uncertain significance |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_vus | BRCA1 | c.299A&gt;T | p.(Glu100Val) | Uncertain significance | criteria provided, multiple submitters, no conflicts | 2 | Color Diagnostics, LLC DBA Color Health: Uncertain significance; Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 85 | tier1_urgent | panel-level public assertion | module1_vus_public_pathogenic_direction | BRCA1 | c.4094T&gt;C | p.(Leu1365Ser) | Pathogenic | reviewed by expert panel | 5 | Evidence-based Network for the Interpretation of Germline Mutant Alleles (ENIGMA): Pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Pathogenic |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4095A&gt;C | p.(Leu1365Phe) | Uncertain significance | criteria provided, single submitter | 1 | Ambry Genetics: Uncertain significance |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_benign_direction | BRCA1 | c.4095A&gt;G | p.(Leu1365=) | Likely benign | criteria provided, single submitter | 1 | Labcorp Genetics (formerly Invitae), Labcorp: Likely benign |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4095A&gt;T | p.(Leu1365Phe) | Uncertain significance | criteria provided, single submitter | 1 | Ambry Genetics: Uncertain significance |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_vus | BRCA1 | c.4096G&gt;A | p.(Gly1366Ser) | Uncertain significance | criteria provided, multiple submitters, no conflicts | 6 | Color Diagnostics, LLC DBA Color Health: Uncertain significance; Ambry Genetics: Uncertain significance; GeneDx: Uncertain significance; Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_vus | BRCA1 | c.4096G&gt;C | p.(Gly1366Arg) | Uncertain significance | criteria provided, multiple submitters, no conflicts | 2 | Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance; Ambry Genetics: Uncertain significance |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4096G&gt;T | p.(Gly1366Cys) | Uncertain significance | criteria provided, single submitter | 1 | Ambry Genetics: Uncertain significance |  |
| 85 | tier1_urgent | panel-level public assertion | module1_vus_public_pathogenic_direction | BRCA1 | c.4097G&gt;C | p.(Gly1366Ala) | Pathogenic | reviewed by expert panel | 5 | GeneDx: Pathogenic; Evidence-based Network for the Interpretation of Germline Mutant Alleles (ENIGMA): Pathogenic; Color Diagnostics, LLC DBA Color Health: Pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Pathogenic |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_pathogenic_direction | BRCA1 | c.4097G&gt;T | p.(Gly1366Val) | Pathogenic | criteria provided, multiple submitters, no conflicts | 30 | Baylor Genetics: Pathogenic; Baylor Genetics: Pathogenic; Color Diagnostics, LLC DBA Color Health: Pathogenic; GeneDx: Pathogenic; Ambry Genetics: Pathogenic; Quest Diagnostics Nichols Institute San Juan Capistrano: pathogenic; Laboratory for Molecular Medicine, Mass General Brigham Personalized Medicine: Pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Pathogenic |  |
| 85 | tier1_urgent | conflicting public assertions | public_conflict | BRCA1 | c.4184A&gt;C | p.(Gln1395Pro) | Conflicting classifications of pathogenicity | criteria provided, conflicting classifications | 10 | GeneDx: Uncertain significance; Baylor Genetics: Uncertain significance; Color Diagnostics, LLC DBA Color Health: Uncertain significance; Ambry Genetics: Uncertain significance; Labcorp Genetics (formerly Invitae), Labcorp: Likely benign |  |
| 85 | tier1_urgent | conflicting public assertions | public_conflict | BRCA1 | c.4184A&gt;G | p.(Gln1395Arg) | Conflicting classifications of pathogenicity | criteria provided, conflicting classifications | 10 | GeneDx: Uncertain significance; Baylor Genetics: Uncertain significance; Color Diagnostics, LLC DBA Color Health: Uncertain significance; Ambry Genetics: Uncertain significance; Labcorp Genetics (formerly Invitae), Labcorp: Likely benign |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4184A&gt;T | p.(Gln1395Leu) | Uncertain significance | criteria provided, single submitter | 1 | Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 85 | tier1_urgent | ClinGen/ENIGMA assertion | module1_vus_public_pathogenic_direction | BRCA1 | c.4185G&gt;C | p.(Gln1395His) | Likely pathogenic | reviewed by expert panel | 4 | Ambry Genetics: Likely pathogenic; Labcorp Genetics (formerly Invitae), Labcorp: Likely pathogenic; ClinGen ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel, ClinGen: Likely Pathogenic | Likely Pathogenic |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4185G&gt;T | p.(Gln1395His) | Uncertain significance | criteria provided, single submitter | 1 | Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_vus | BRCA1 | c.4186C&gt;A | p.(Gln1396Lys) | Uncertain significance | criteria provided, multiple submitters, no conflicts | 4 | GeneDx: Uncertain significance; Ambry Genetics: Uncertain significance; Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4186C&gt;G | p.(Gln1396Glu) | Uncertain significance | criteria provided, single submitter | 1 | Ambry Genetics: Uncertain significance |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_vus | BRCA1 | c.4187A&gt;C | p.(Gln1396Pro) | Uncertain significance | criteria provided, multiple submitters, no conflicts | 2 | Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance; Ambry Genetics: Uncertain significance |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4187A&gt;G | p.(Gln1396Arg) | Uncertain significance | criteria provided, single submitter | 1 | Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4187A&gt;T | p.(Gln1396Leu) | Uncertain significance | criteria provided, single submitter | 1 | Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_benign_direction | BRCA1 | c.4188G&gt;A | p.(Gln1396=) | Likely benign | criteria provided, multiple submitters, no conflicts | 2 | Ambry Genetics: Likely benign; Labcorp Genetics (formerly Invitae), Labcorp: Likely benign |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_vus | BRCA1 | c.4188G&gt;C | p.(Gln1396His) | Uncertain significance | criteria provided, single submitter | 1 | Quest Diagnostics Nichols Institute San Juan Capistrano: Uncertain significance |  |
| 85 | tier1_urgent | single-submitter assertion | module1_vus_public_benign_direction | BRCA1 | c.4188G&gt;T | p.(Gln1396His) | Likely benign | criteria provided, single submitter | 1 | Labcorp Genetics (formerly Invitae), Labcorp: Likely benign |  |
| 85 | tier1_urgent | multi-submitter assertion | module1_vus_public_vus | BRCA1 | c.4357G&gt;A | p.(Ala1453Thr) | Uncertain significance | criteria provided, multiple submitters, no conflicts | 2 | Labcorp Genetics (formerly Invitae), Labcorp: Uncertain significance; Ambry Genetics: Uncertain significance |  |

## Interpretation Boundary

This is a snapshot of public assertions, not a benchmark. A public pathogenic
or benign assertion is a reason to inspect the variant, submitter history,
evidence, disease context, and ClinGen/ENIGMA status. It should not be counted
as automatic ACMG/ENIGMA evidence unless the underlying criterion and source
requirements are met.

## Outputs

- `tables/public_classification_snapshot/public_classification_snapshot_variants.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_category.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_discordance.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_clinvar_status.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_clingen_status.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_category_action_group.csv`
- `external_sources/public_classification_snapshot_sources/public_classification_snapshot_cache.json`
