# VUS Evidence Action Plan

Generated: 2026-06-22 12:10

## Purpose

This analysis converts unresolved ARIANE Module 1 VUS bottlenecks into a
practical evidence worklist. It does not add ACMG/ENIGMA points and it does
not treat local context as evidence. The aim is to ask: for this VUS group,
which kind of additional variant-level evidence would be most useful?

## Action Groups

| action_group | count | high_priority_count | high_spliceai_count | median_points | top_candidate_evidence |
| --- | --- | --- | --- | --- | --- |
| background_absence_only | 5048 | 0 | 10 | 1.0 | any independent variant-level evidence; otherwise low priority |
| conflicting_computational_splice_context | 1276 | 0 | 0 | -1.0 | manual RNA/splicing review and transcript-specific interpretation |
| computational_signal_needs_noncomputational_support | 842 | 812 | 842 | 2.0 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| near_pathogenic_threshold | 446 | 38 | 5 | 5.0 | supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion |
| module1_no_signal | 331 | 0 | 1 | 0 | ClinVar/ClinGen/ENIGMA assertion context, literature, segregation, functional data |
| near_benign_threshold | 191 | 2 | 0 | -4 | BS2, BS4, population evidence, or curated functional strength review |
| manual_conflict_adjudication | 20 | 7 | 8 | -0.5 | manual review of mixed benign and pathogenic evidence |

## Bottleneck To Action Mapping

| bottleneck_category | action_group | count | high_priority_count | high_spliceai_count | top_candidate_evidence |
| --- | --- | --- | --- | --- | --- |
| PM2_only | background_absence_only | 5046 | 0 | 10 | any independent variant-level evidence; otherwise low priority |
| mixed_splice_benign_and_pathogenic | conflicting_computational_splice_context | 1276 | 0 | 0 | manual RNA/splicing review and transcript-specific interpretation |
| computational_pathogenic_evidence_not_enough | computational_signal_needs_noncomputational_support | 842 | 812 | 842 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| strong_pathogenic_combo_one_step_short | near_pathogenic_threshold | 446 | 38 | 5 | supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion |
| no_automated_evidence | module1_no_signal | 331 | 0 | 1 | ClinVar/ClinGen/ENIGMA assertion context, literature, segregation, functional data |
| benign_functional_evidence_not_enough | near_benign_threshold | 170 | 2 | 0 | BS2, BS4, population evidence, or curated functional strength review |
| benign_evidence_not_enough | near_benign_threshold | 21 | 0 | 0 | BS2, BS4, population evidence, or ClinGen/ENIGMA assertion context |
| other_conflicting_evidence | manual_conflict_adjudication | 19 | 7 | 8 | manual review of mixed benign and pathogenic evidence |
| PM2_plus_weak_context | background_absence_only | 2 | 0 | 0 | stronger independent evidence; do not classify from context alone |
| conflicting_PVS1_BS3 | manual_conflict_adjudication | 1 | 0 | 0 | functional assay metadata, NMD/splicing review, ClinGen/ENIGMA assertion context |

## Per Gene

| gene | action_group | count | high_priority_count | high_spliceai_count | top_candidate_evidence |
| --- | --- | --- | --- | --- | --- |
| BRCA2 | background_absence_only | 4322 | 0 | 1 | any independent variant-level evidence; otherwise low priority |
| BRCA2 | conflicting_computational_splice_context | 1186 | 0 | 0 | manual RNA/splicing review and transcript-specific interpretation |
| BRCA1 | background_absence_only | 726 | 0 | 9 | any independent variant-level evidence; otherwise low priority |
| BRCA2 | computational_signal_needs_noncomputational_support | 556 | 526 | 556 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| BRCA1 | near_pathogenic_threshold | 342 | 33 | 4 | supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion |
| BRCA1 | computational_signal_needs_noncomputational_support | 286 | 286 | 286 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| BRCA2 | module1_no_signal | 280 | 0 | 0 | ClinVar/ClinGen/ENIGMA assertion context, literature, segregation, functional data |
| BRCA1 | near_benign_threshold | 126 | 2 | 0 | BS2, BS4, population evidence, or curated functional strength review |
| BRCA2 | near_pathogenic_threshold | 104 | 5 | 1 | supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion |
| BRCA1 | conflicting_computational_splice_context | 90 | 0 | 0 | manual RNA/splicing review and transcript-specific interpretation |
| BRCA2 | near_benign_threshold | 65 | 0 | 0 | BS2, BS4, population evidence, or curated functional strength review |
| BRCA1 | module1_no_signal | 51 | 0 | 1 | ClinVar/ClinGen/ENIGMA assertion context, literature, segregation, functional data |
| BRCA2 | manual_conflict_adjudication | 16 | 4 | 5 | manual review of mixed benign and pathogenic evidence |
| BRCA1 | manual_conflict_adjudication | 4 | 3 | 3 | manual review of mixed benign and pathogenic evidence |

## Top Candidates

| priority_score | priority_tier | action_group | bottleneck_category | gene | c_notation | p_notation | criteria_combo | spliceai_score | candidate_evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 90 | tier1_urgent | near_pathogenic_threshold | strong_pathogenic_combo_one_step_short | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | PP3+PS3 | 0.65 | supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion |
| 90 | tier1_urgent | manual_conflict_adjudication | other_conflicting_evidence | BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | BP5+PM2_Supporting+PP3+PS3 | 0.57 | manual review of mixed benign and pathogenic evidence |
| 90 | tier1_urgent | near_pathogenic_threshold | strong_pathogenic_combo_one_step_short | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | PP3+PS3 | 0.22 | supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion |
| 90 | tier1_urgent | near_pathogenic_threshold | strong_pathogenic_combo_one_step_short | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | PP3+PS3 | 0.4 | supporting pathogenic evidence such as PP1, PS4, PM3, RNA evidence, or curated clinical assertion |
| 90 | tier1_urgent | manual_conflict_adjudication | other_conflicting_evidence | BRCA1 | c.5453A&gt;G | p.(Asp1818Gly) | BP5+PM2_Supporting+PP3+PS3 | 0.35 | manual review of mixed benign and pathogenic evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.299A&gt;G | p.(Glu100Gly) | PM2_Supporting+PP3 | 0.56 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.299A&gt;T | p.(Glu100Val) | PM2_Supporting+PP3 | 0.28 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4094T&gt;C | p.(Leu1365Ser) | PM2_Supporting+PP3 | 0.27 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4095A&gt;C | p.(Leu1365Phe) | PM2_Supporting+PP3 | 0.22 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4095A&gt;G | p.(Leu1365=) | PM2_Supporting+PP3 | 0.37 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4095A&gt;T | p.(Leu1365Phe) | PM2_Supporting+PP3 | 0.3 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4096G&gt;A | p.(Gly1366Ser) | PM2_Supporting+PP3 | 0.6 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4096G&gt;C | p.(Gly1366Arg) | PM2_Supporting+PP3 | 0.69 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4096G&gt;T | p.(Gly1366Cys) | PM2_Supporting+PP3 | 0.68 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4097G&gt;C | p.(Gly1366Ala) | PM2_Supporting+PP3 | 0.21 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4097G&gt;T | p.(Gly1366Val) | PM2_Supporting+PP3 | 0.36 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4184A&gt;C | p.(Gln1395Pro) | PM2_Supporting+PP3 | 0.28 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | manual_conflict_adjudication | other_conflicting_evidence | BRCA1 | c.4184A&gt;G | p.(Gln1395Arg) | BP5+PM2_Supporting+PP3 | 0.36 | manual review of mixed benign and pathogenic evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4184A&gt;T | p.(Gln1395Leu) | PM2_Supporting+PP3 | 0.46 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4185G&gt;C | p.(Gln1395His) | PM2_Supporting+PP3 | 0.95 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4185G&gt;T | p.(Gln1395His) | PM2_Supporting+PP3 | 0.96 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4186C&gt;A | p.(Gln1396Lys) | PM2_Supporting+PP3 | 0.48 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4186C&gt;G | p.(Gln1396Glu) | PM2_Supporting+PP3 | 0.87 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4187A&gt;C | p.(Gln1396Pro) | PM2_Supporting+PP3 | 0.87 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4187A&gt;G | p.(Gln1396Arg) | PP3 | 0.87 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4187A&gt;T | p.(Gln1396Leu) | PM2_Supporting+PP3 | 0.87 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4188G&gt;A | p.(Gln1396=) | PM2_Supporting+PP3 | 0.87 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4188G&gt;C | p.(Gln1396His) | PM2_Supporting+PP3 | 0.87 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4188G&gt;T | p.(Gln1396His) | PM2_Supporting+PP3 | 0.87 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |
| 85 | tier1_urgent | computational_signal_needs_noncomputational_support | computational_pathogenic_evidence_not_enough | BRCA1 | c.4357G&gt;A | p.(Ala1453Thr) | PM2_Supporting+PP3 | 0.27 | RNA, functional, PS4, PP1, PM3, or curated assertion evidence |

## Interpretation

The largest group is usually `background_absence_only`, meaning PM2 alone or
PM2 with weak context. These variants are not the most efficient first manual
review target unless external public assertions, literature, or laboratory data
already exist.

The most actionable groups are:

1. `near_pathogenic_threshold`, where one independent pathogenic evidence item
   may be enough to move the case from VUS toward likely pathogenic.
2. `near_benign_threshold`, where stronger benign evidence or an additional
   benign criterion may resolve the case.
3. `computational_signal_needs_noncomputational_support`, where SpliceAI/PP3
   points to a hypothesis but cannot classify the variant by itself.

## Outputs

- `tables/vus_evidence_action_plan/vus_evidence_action_plan_variants.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_action_group.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_bottleneck.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_gene.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_by_gene_exon.csv`
- `tables/vus_evidence_action_plan/vus_evidence_action_plan_top_candidates.csv`
- `plots/29_vus_evidence_action_plan/vus_evidence_action_groups.svg`
