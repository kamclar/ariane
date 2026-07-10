# Precomputed BRCA1/2 Classification Exploration

This directory contains exploratory analyses of the precomputed ARIANE Module 1
coding SNV variant map for BRCA1 and BRCA2.

The analyses describe what the current automatable Module 1 rule set does
across the precomputed variant space. They do not replace manual classification
and they do not provide an independent clinical reference set.

## Writing And Interpretation Note

When these analyses are later turned into continuous text, two layers must be
kept separate:

1. Analysis of the generated classification itself.
   This asks how ARIANE Module 1 criteria behave: which criteria are applied,
   which combinations create class 1-5 outputs, which criteria hold a class in
   place, and which VUS groups are one evidence item away from moving.

2. Analysis of patterns in the generated data landscape.
   This asks what seems to happen across positions, domains, exons, splice
   boundaries, pathogenic-density regions, or mechanism buckets. These patterns
   are useful for understanding the model output and for generating hypotheses,
   but they are not the same as independent biological or clinical evidence.

The current dataset is a synthetic, precomputed coding SNV map classified by
the current ARIANE automatable Module 1 criteria. It should be described as a
map of precomputed variants classified by automated rules, not as a map of
biological reality. Therefore, statements such as "pathogenic density",
"pathogenic mechanism", "hotspot", "domain enrichment", or "VUS priority"
should be written as patterns in the generated ARIANE classification snapshot
unless they are independently validated against external curated evidence.

Local pathogenic enrichment is a context signal, not an ACMG/ENIGMA evidence
criterion.

## Current Input

Primary input:

- `variant_space_scan/outputs/brca_module1_full_snv_classification.csv`

Scope:

- BRCA1 and BRCA2
- coding SNVs from the current precomputed manifest
- ARIANE Module 1 automated ENIGMA/ACMG criteria
- local SpliceAI cache integrated into the classification snapshot

Observed automated criteria in this snapshot:

- `PM2_Supporting`, `BA1`, `BS1_Supporting`, `BS1_Strong`
- `PP3`, `BP4`, `BP7`
- `PVS1`, `PM5_PTC`, `BP1`
- `PS3`, `BS3`, `PS1`, `PP4`, `BP5`

Important limitation:

- The current exploration does not yet cover indels, multi-exon CNVs, deep
  intronic variants, promoter/UTR variants, or other noncoding regions outside
  this SNV snapshot.

## What Is Already Done

| Area | Script | Report | Main output |
| --- | --- | --- | --- |
| Basic summary | `analyze_precomputed_classification.py` | `precomputed_classification_analysis.md` | class, gene, criterion, variant-type summaries |
| Visual overview | `make_visualizations.py` | `visualization_report.md` | overview plots, position plots, clustered profiles |
| Hotspots and coldspots | `find_hotspots.py` | `hotspot_coldspot_report.md` | genomic regions enriched or depleted for selected signals |
| Exon-level VUS and conflict summary | `exon_vus_conflict_analysis.py` | `exon_vus_conflict_report.md` | exon summaries, VUS priorities, conflict checks |
| BRCA1/BRCA2 comparison | `gene_comparison_analysis.py` | `gene_comparison_report.md` | normalized profiles and gene-level composition |
| VUS prioritization | `vus_prioritization_analysis.py` | `vus_prioritization_report.md` | ranked VUS candidates and priority reasons |
| Boundary enrichment | `boundary_enrichment_analysis.py` | `boundary_enrichment_report.md` | splice-boundary enrichment and threshold sensitivity |
| Criteria co-occurrence | `criteria_cooccurrence_analysis.py` | `criteria_cooccurrence_report.md` | criterion combinations and pairwise co-occurrence |
| Criteria sanity audit | `criteria_sanity_audit.py` | `criteria_sanity_audit_report.md` | expected edge cases and possible review points |
| PVS1 but VUS audit | manual/report | `pvs1_but_vus_case_audit.md` | focused review of PVS1-plus-conflict behavior |
| VUS bottleneck | `vus_bottleneck_analysis.py` | `vus_bottleneck_report.md` | why variants remain VUS and what evidence is missing |
| VUS bottleneck deep dive | `vus_bottleneck_deep_dive.py` | `vus_bottleneck_deep_dive_report.md` | higher-resolution VUS buckets and action queue |
| Curation review queue | `curation_review_queue.py` | `curation_review_queue.md` | practical review tiers for curator follow-up |
| Criterion sensitivity | `criterion_sensitivity_analysis.py` | `criterion_sensitivity_report.md` | criteria that hold classifications in place and VUS evidence-impact simulation |
| Frequency baseline audit | `frequency_baseline_audit.py` | `frequency_baseline_audit_report.md` | current PM2/BS1/BA1 dependence before future gnomAD v4.1 calibration |
| Pathogenic mechanism exploration | `pathogenic_mechanism_analysis.py` | `pathogenic_mechanism_report.md` | class 4/5 drivers by truncation, functional evidence, domains, and pathogenic residue context |
| Structure/function mapping | `structure_function_mapping.py` | `structure_function_mapping_report.md` | non-truncating no-splice pathogenic variants mapped to broad protein structure/function features |
| VUS pathogenic-region proximity | `vus_pathogenic_region_analysis.py` | `vus_pathogenic_region_report.md` | VUS in pathogenic regions and VUS closer to pathogenic than benign local neighborhoods |
| BRCA1 BRCT1 mixed cluster | `brca1_brct1_mixed_cluster_analysis.py` | `brca1_brct1_mixed_cluster_report.md` | focused explanation of the BRCA1 c.5044-c.5143 mixed VUS cluster seen in the +/-20 bp neighborhood plot |
| Regional driver decomposition | `regional_driver_decomposition_analysis.py` | `regional_driver_decomposition_report.md` | decomposes regional pathogenic enrichment into truncation, splice, functional, benign, and residual signals |
| Benign structure/function mapping | `benign_structure_function_analysis.py` | `benign_structure_function_report.md` | benign/likely benign variants mapped to the same domains and UniProt active/interface annotations |
| Benign counterexamples | `benign_counterexamples_analysis.py` | `benign_counterexamples_report.md` | benign variants inside important domains or active/interface contexts as negative controls |
| BS3 domain conflicts | `bs3_domain_conflict_analysis.py` | `bs3_domain_conflict_report.md` | BS3 benign evidence inside functional domains and conflicts with PP3 or high SpliceAI |
| BS3 functional source audit | `audit_bs3_functional_sources.py` | `bs3_functional_source_audit_report.md` | audits BRCA1 BS3 RING/BRCT evidence against downloaded Findlay 2018 SGE and Dace/Findlay 2023 update sources |
| Position/codon class conflicts | `position_class_conflict_analysis.py` | `position_class_conflict_report.md` | positions and codons where different possible SNVs receive different grouped classes |
| Positional context interpretation | manual/report | `positional_context_interpretation.md` | conceptual interpretation of what local pathogenic enrichment can and cannot mean |
| VUS prioritization exploratory report | manual/report | `vus_prioritization_module1_exploratory_report.md` | working Module 1 summary of VUS review candidate prioritization |
| Budapest 2026 notes | manual/report | `budapest2026_notes.md` | conference notes for gnomAD v4.1 calibration, MAVE functional evidence, VRS, and VA-Spec |
| Positional context follow-up | `positional_context_followup_analyses.py` | `positional_context_followup_report.md` | tests whether variant type, splice signal, functional evidence, and mixed positions/codons explain local context patterns |
| VUS manuscript critique response | `vus_manuscript_critique_response_analysis.py` | `vus_manuscript_critique_response_report.md` | quantifies how much the VUS priority tiers depend on local-neighborhood score components |
| Findlay SGE VUS tier analysis | `findlay_sge_vus_tier_analysis.py` | `findlay_sge_vus_tier_report.md` | compares BRCA1 VUS priority tiers against continuous Findlay 2018 SGE function scores where available |
| VUS evidence action plan | `vus_evidence_action_plan_analysis.py` | `vus_evidence_action_plan_report.md` | converts VUS bottlenecks into practical next-evidence review groups |
| Public classification snapshot | `public_classification_snapshot_analysis.py` | `public_classification_snapshot_report.md` | queries ClinVar and ClinGen/ENIGMA context for highest-priority VUS as a public assertion worklist |

Plots are grouped under:

- `plots/01_overview`
- `plots/02_position`
- `plots/03_boundary_spliceai`
- `plots/04_clusters`
- `plots/05_hotspots`
- `plots/06_exon_vus_conflict`
- `plots/07_gene_comparison`
- `plots/08_vus_prioritization`
- `plots/09_boundary_enrichment`
- `plots/10_criteria_cooccurrence`
- `plots/11_criteria_sanity_audit`
- `plots/12_vus_bottleneck`
- `plots/13_vus_bottleneck_deep_dive`
- `plots/15_frequency_baseline`
- `plots/18_vus_pathogenic_regions`
- `plots/19_benign_structure_function`
- `plots/20_benign_counterexamples`
- `plots/21_bs3_domain_conflicts`
- `plots/22_position_class_conflicts`
- `plots/23_brca1_brct1_mixed_cluster`
- `plots/24_position_context_concept`
- `plots/25_regional_driver_decomposition`
- `plots/26_positional_context_followup`
- `plots/27_vus_manuscript_critique_response`
- `plots/28_findlay_sge_vus_tier`
- `plots/29_vus_evidence_action_plan`
- `plots/16_pathogenic_mechanisms`
- `plots/17_structure_function_mapping`
- `plots/14_criterion_sensitivity`

Tables are grouped under `tables/` with matching subdirectories where the
analysis has multiple outputs.

## Main Findings So Far

- The strongest current claim is not that ARIANE can predict VUS
  pathogenicity. The stronger and safer claim is that a precomputed Module 1
  map can be turned into an auditable review queue that separates variant-level
  signals, neighborhood-driven context, and evidence-conflict cases without
  treating those signals as classification criteria.
- Most analyses confirm that the landscape is largely driven by the automated
  Module 1 rules and especially by SpliceAI/PP3, PVS1/PM5_PTC, BP1, BS3, and
  PM2 behavior.
- The largest practical VUS bottleneck is not a small number of dramatic
  outliers, but many low-information variants where PM2 is the only evidence.
- The most useful curation queue is therefore not simply the largest VUS group.
  It is the group where one credible additional evidence item could change the
  interpretation.
- The clearest high-priority VUS group is `PS3+PP3` one-step-short pathogenic
  VUS.
- The main conflict group is benign functional evidence versus computational or
  truncating pathogenic evidence. These are review targets, not automatic
  software errors.
- Boundary analyses show that SpliceAI signal is strongly tied to splice
  proximity, but distance alone is not enough to separate benign and pathogenic
  grouped classes.
- Structure/function analyses show that domain or interface context is useful
  for prioritization and interpretation, but it is not evidence by itself.
  Benign/likely benign counterexamples occur inside RING, BRCT, BRCA2 DBD,
  PALB2, RAD51/BRC, and SEM1/DSS1 contexts. Therefore, a VUS in one of these
  regions should be reviewed earlier only when the local neighborhood, criteria,
  and residue-level evidence point in the same direction.
- BS3 explains many benign variants inside functional domains, especially BRCA1
  BRCT and RING in this snapshot. The important review subset is smaller:
  variants with `BS3` plus `PP3` or high SpliceAI, where the functional assay
  metadata matters because the assay may or may not capture splicing, NMD, or
  the relevant protein mechanism.
- External-source audit shows that BRCA1 `BS3` in RING/BRCT is not just broad
  domain-level annotation. In ENIGMA Table 9, 1861/1865 BRCA1 BRCT `BS3`
  records and 868/871 N-terminal RING-region `BS3` records match variant-level
  rows in Findlay et al. 2018 saturation genome editing Supplementary Table 1.
  A small RING subset has older Findlay `LOF` labels while local ENIGMA Table 9
  cites a Dace/Findlay 2023 interim update that reinterprets those variants as
  benign-like. That interim source is recorded as a Table 9 citation here, not
  independently verified as a standalone publication; these remain manual audit
  targets.
- Position/codon conflict analysis confirms that a position is not equivalent
  to a concrete variant. Many coding positions and codons contain different
  possible SNVs with different generated grouped classes. The dominant expected
  mechanism is often truncating versus missense/synonymous consequence, but VUS
  in fully mixed contexts are useful manual-review candidates.
- A visible VUS cluster in the +/-20 bp local-neighborhood plot maps to BRCA1
  c.5044-c.5143, approximately p.1682-p.1715, inside the BRCT1 subdomain.
  It is benign-enriched but mixed: the region contains many class 1/2 variants,
  many VUS, and a smaller pathogenic group. Its behavior is driven by a mixture
  of BRCT functional evidence, BS3/PS3, local PP3/SpliceAI signal near a
  donor-like boundary, and nearby truncating substitutions.
- Regional driver decomposition shows that local pathogenic enrichment is
  largely explained by truncation-related criteria, splice/PP3 signal, and
  functional evidence. After excluding truncation, splice, and functional
  evidence drivers together, the residual pathogenic fraction is 0 in the
  current 100 bp bin analysis. This supports the interpretation that positional
  context is useful for triage, not as a standalone evidence criterion.
- Positional context follow-up confirms the same point from another angle.
  Nonsense/PTC variants are overwhelmingly pathogenic in the generated map,
  while missense and synonymous variants are mostly benign or VUS. Mixed
  position and codon contexts almost disappear after removing truncation,
  splice, and functional-evidence drivers: positions with both benign and
  pathogenic grouped classes drop from 1654 to 0, and mixed codons with both
  benign and pathogenic grouped classes drop from 1737 to 0.
- Manuscript critique response analysis shows that the neighborhood component
  has a large effect on the current high-priority queue. There are 859
  tier1/tier2 VUS in the original heuristic score, but only 286 remain
  tier1/tier2 after removing local-neighborhood points. This means
  neighborhood-driven VUS must be labelled as context-driven review candidates,
  not as independently high-priority variants.
- Findlay 2018 SGE provides a continuous functional axis for 673 BRCA1 VUS in
  RING/BRCT-covered regions, but it is not fully independent because local
  BRCA1 `PS3`/`BS3` evidence is partly derived from the same source. The
  supported interpretation is coarse: tier4 is clearly more function-preserved
  or benign-leaning (median `function.score.mean` -0.828; 14 LOF, 139 INT, 100
  FUNC), but the data do not validate a monotonic tier1 > tier2 > tier3 order.
  Tier3 contains 303 LOF-class variants, 302 of them with `PS3`, usually
  `PM2_Supporting+PS3`. This suggests that the heuristic may under-prioritize
  functionally damaging missense VUS that have PS3 but lack splice, boundary,
  PS1, PP4, or stronger neighborhood signals. Mechanically, `PS3` has 25
  review-priority points and tier3 starts at 25, so a PS3-driven VUS lands at
  the floor of tier3 by construction. Of the 303 tier3 LOF variants, 302 overlap
  with `near_pathogenic_threshold` / `strong_pathogenic_combo_one_step_short`.
  This shows that priority tier and distance to likely pathogenic are related
  but not identical.
- The evidence action plan turns bottlenecks into review questions. The largest
  group remains low-information/background absence only (5048 VUS), but the
  most actionable groups are `computational_signal_needs_noncomputational_support`
  (842 VUS, 812 high priority), `near_pathogenic_threshold` (446 VUS), and
  `near_benign_threshold` (191 VUS). These groups point to different evidence
  searches instead of one generic VUS queue. The practical workflow should not
  start with tier1/tier2 alone; it should also pull near-pathogenic-threshold
  variants from tier3, especially `PM2_Supporting+PS3` variants.
- The public classification snapshot for 80 highest-priority VUS found 28
  multi-submitter assertions, 28 single-submitter assertions, 11 conflicting
  public assertions, 7 variants with no public assertion, 4 panel-level public
  assertions, and 2 ClinGen/ENIGMA assertions. ClinVar lookup succeeded for
  73/80 variants and ClinGen Evidence Repository lookup found 2/80. These are
  public clinical/laboratory context labels, not a clinical reference set. They
  identify a discordance worklist: 13 Module 1 VUS have public
  pathogenic-direction assertions, 7 have public benign-direction assertions,
  and 11 have public conflicts. The `no public assertion` source category has
  7 variants, while the `no_public_classification` worklist label has 8 because
  one ClinVar record has aggregate classification `not provided`.
- ClinGen ERepo has low BRCA1/2 coverage by design in the current public table:
  BRCA1 has 75 records (9 pathogenic, 1 likely pathogenic, 16 VUS, 33 likely
  benign, 16 benign) and BRCA2 has 68 records (12 pathogenic, 6 likely
  pathogenic, 14 VUS, 26 likely benign, 10 benign). Therefore, ERepo `not_found`
  is a low-coverage-source result, not evidence that no public interpretation
  exists. Entry point: https://erepo.clinicalgenome.org/evrepo/

## What Is Still Missing

### 1. Independent Clinical Comparison

The current exploration mostly studies ARIANE output against ARIANE-derived
features. This tells us how the automated system behaves, but not whether each
prediction matches independent clinical or laboratory assertions.

Useful next validation sets:

- ENIGMA/ClinGen assertions where available
- ClinVar variants with strong review status
- known benchmark examples from EQA material
- local laboratory curated variants, if available

Goal:

- compare ARIANE Module 1 class against an external curated class
- flag systematic mismatches
- separate expected Module 1 limitations from possible implementation problems

### 2. Reproducibility Metadata

The exploration needs a small machine-readable manifest recording:

- input snapshot path
- input row count
- generation date
- relevant cache/data versions
- SpliceAI source used
- script names and output paths

This will make it clear which reports belong to which precomputed snapshot.

### 3. One-Command Rebuild

At the moment the analyses are separate scripts. A small runner would help:

- execute the analyses in the intended order
- stop on errors
- write a run log
- optionally skip expensive plots

### 4. Non-SNV And Noncoding Extension

The current exploration is coding SNV-centered. Future analyses should be
separate, because the evidence model and variant generation differ.

Missing spaces:

- coding indels
- splice-region indels
- canonical splice-site variants not represented as coding SNVs
- UTR variants
- deep intronic variants
- exon-level and multi-exon CNVs

### 5. Application Integration Checklist

Some exploration results already affected the application through the VUS
explanation layer. We should keep a short checklist of which findings are:

- already implemented in ARIANE
- suitable for UI display only
- suitable for curation queue export
- not suitable for automated classification
- blocked until external validation exists

### 6. Compact Executive Summary

There are many reports now. It would help to have one short human-readable
summary with:

- what was analyzed
- the strongest findings
- what changed in the app
- what remains research-only
- what should be validated before clinical-facing use

### 7. Priority Score Sensitivity And Distance-To-Threshold

The Findlay SGE analysis showed that direct functional `PS3` evidence can be
under-prioritized by the current review score. A useful next iteration would:

- add a separate `distance_to_LP` or `one_step_short` field next to priority
  tier
- test priority-score sensitivity to increasing `PS3` weight from 25 to at
  least the high-SpliceAI weight of 30
- report variants where review priority and distance to likely pathogenic
  disagree
- keep any score change limited to review triage, not automatic classification

### 8. Review-Yield Pilot

To support a stronger practical claim, run a small manual pilot:

- 30 tier1/tier2 variants
- 30 tier3 variants
- 30 tier4 or PM2-only variants

For each variant, record whether manual review finds actionable external
evidence, such as RNA, functional assay details, segregation, case-control
data, trans observations, public ClinVar/ClinGen/ENIGMA assertions, or local
laboratory observations. The desired output is a review-yield estimate, for
example whether high-priority tiers are enriched for variants with actionable
external evidence compared with low-priority PM2-only VUS.

## Recommended Next Step

The next most useful technical step is reproducibility metadata plus an
exploration manifest. The next most useful scientific step is a small
review-yield pilot plus priority-score sensitivity and independent clinical
comparison against ENIGMA/ClinVar/local curated variants.

If the goal is application development, prioritize the application integration
checklist. If the goal is biology/variant interpretation, prioritize external
validation and non-SNV/noncoding extension.
