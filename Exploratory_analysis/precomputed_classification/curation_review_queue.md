# Curator Review Queue

Generated: 2026-06-20 09:50

Inputs:

- `tables/vus_bottleneck_deep_dive/deep_dive_action_queue.csv`
- `tables/criteria_sanity_audit/sanity_variant_level_audit.csv`

## Purpose

This checklist converts the exploratory VUS analyses into a practical curator review queue. It does not change classifications. It selects variants where manual evidence review is most likely to be useful, and separates them from large low-information groups.

## Tier Summary

| tier | tier_label | count | top_source_bucket | max_spliceai |
| --- | --- | --- | --- | --- |
| A | PS3+PP3 one-step-short pathogenic VUS | 4 | PS3_PP3_splice_or_functional_followup | 0.65 |
| B | Strong pathogenic evidence one step short | 442 | PS3_PM2_needs_one_more_pathogenic_support | 0.23 |
| C | Mixed benign and pathogenic evidence | 11 | pathogenic_and_benign_evidence | 0.83 |
| C | PP3 plus BS3 conflict | 128 | pathogenic_and_benign_evidence;PP3_with_BS3 | 0.99 |
| C | PVS1 plus benign functional conflict | 1 | pathogenic_and_benign_evidence;PVS1_but_VUS | 0.01 |
| D | Benign-leaning BP4/BP7 plus PM2 VUS | 1276 | BP4_BP7_PM2_near_pathogenic_density | 0.10 |
| E | PM2-only background VUS | 5046 | PM2_only_pathogenic_neighborhood | 0.68 |

## How To Read The Tiers

| tier | meaning |
| --- | --- |
| A | highest priority pathogenic-side VUS: PS3 and PP3 point in the same direction, but the ENIGMA combination is incomplete |
| B | strong pathogenic evidence one step short, often PS3 plus PM2 |
| C | evidence conflicts, such as PVS1 plus BS3 or PP3 plus BS3 |
| D | benign-leaning splice VUS, usually BP4/BP7 plus PM2 |
| E | PM2-only background VUS, lowest priority unless external evidence appears |

## Tier A: PS3+PP3 One-Step-Short

Why: functional and computational/splice evidence point toward pathogenicity, but another independent accepted evidence item is still needed.

| tier | gene | c_notation | p_notation | criteria_combo | total_points | spliceai_score | what_to_check |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | BRCA1 | c.190T&gt;G | p.(Cys64Gly) | PP3+PS3 | 5 | 0.65 | Look for RNA confirmation, PP1, PM3, PS4, PP4, curated PS1, or other accepted independent pathogenic evidence. |
| A | BRCA1 | c.5123C&gt;G | p.(Ala1708Gly) | PP3+PS3 | 5 | 0.4 | Look for RNA confirmation, PP1, PM3, PS4, PP4, curated PS1, or other accepted independent pathogenic evidence. |
| A | BRCA1 | c.5123C&gt;A | p.(Ala1708Glu) | PP3+PS3 | 5 | 0.22 | Look for RNA confirmation, PP1, PM3, PS4, PP4, curated PS1, or other accepted independent pathogenic evidence. |
| A | BRCA2 | c.9218A&gt;G | p.(Asp3073Gly) | PP3+PS3 | 5 | 0.22 | Look for RNA confirmation, PP1, PM3, PS4, PP4, curated PS1, or other accepted independent pathogenic evidence. |

## Tier B: Strong Pathogenic Evidence One Step Short

Why: these are often close to likely pathogenic, but still need one supporting or moderate pathogenic item under the ENIGMA combination rules.

| tier | gene | c_notation | p_notation | criteria_combo | total_points | spliceai_score | what_to_check |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B | BRCA1 | c.5077G&gt;T | p.(Ala1693Ser) | PM2_Supporting+PS3 | 5 | 0.1 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.3G&gt;A | p.(Met1Ile) | PM2_Supporting+PS3 | 5 | 0.23 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5072C&gt;A | p.(Thr1691Lys) | PM2_Supporting+PS3 | 5 | 0.05 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5145C&gt;A | p.(Ser1715Arg) | PM2_Supporting+PS1+PS3 | 9 | 0.03 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5060T&gt;A | p.(Val1687Asp) | PM2_Supporting+PS3 | 5 | 0.17 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5062G&gt;T | p.(Val1688Phe) | PM2_Supporting+PS3 | 5 | 0.16 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5060T&gt;G | p.(Val1687Gly) | PM2_Supporting+PS3 | 5 | 0.15 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5063T&gt;A | p.(Val1688Asp) | PM2_Supporting+PS3 | 5 | 0.15 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA2 | c.7007G&gt;C | p.(Arg2336Pro) | PS3 | 4 | 0.18 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5059G&gt;T | p.(Val1687Phe) | PM2_Supporting+PS3 | 5 | 0.16 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5058T&gt;A | p.(His1686Gln) | PM2_Supporting+PS3 | 5 | 0.15 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5058T&gt;G | p.(His1686Gln) | PM2_Supporting+PS3 | 5 | 0.15 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.131G&gt;C | p.(Cys44Ser) | PM2_Supporting+PS1+PS3 | 9 | 0.0 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5056C&gt;A | p.(His1686Asn) | PM2_Supporting+PS3 | 5 | 0.19 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5057A&gt;C | p.(His1686Pro) | PM2_Supporting+PS3 | 5 | 0.18 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5056C&gt;G | p.(His1686Asp) | PM2_Supporting+PS3 | 5 | 0.16 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5057A&gt;G | p.(His1686Arg) | PM2_Supporting+PS3 | 5 | 0.15 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5145C&gt;G | p.(Ser1715Arg) | PM2_Supporting+PS3 | 5 | 0.16 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA2 | c.7978T&gt;G | p.(Tyr2660Asp) | PM2_Supporting+PS3 | 5 | 0.03 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5066T&gt;A | p.(Met1689Lys) | PM2_Supporting+PS3 | 5 | 0.04 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5153G&gt;T | p.(Trp1718Leu) | PM2_Supporting+PS3 | 5 | 0.07 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5154G&gt;C | p.(Trp1718Cys) | PM2_Supporting+PS3 | 5 | 0.01 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.5154G&gt;T | p.(Trp1718Cys) | PM2_Supporting+PS3 | 5 | 0.01 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.181T&gt;A | p.(Cys61Ser) | PM2_Supporting+PP4+PS3 | 9 | 0.08 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |
| B | BRCA1 | c.218T&gt;G | p.(Leu73Arg) | PM2_Supporting+PS3 | 5 | 0.16 | Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion. |

## Tier C: Evidence Conflicts

Why: these are not automatic upgrades or downgrades. They need curator adjudication of conflicting evidence directions.

| tier | gene | c_notation | p_notation | criteria_combo | total_points | spliceai_score | what_to_check |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C | BRCA1 | c.5196T&gt;G | p.(His1732Gln) | BS3+PM2_Supporting+PP3 | -2 | 0.99 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.5470A&gt;G | p.(Ile1824Val) | BS3+PM2_Supporting+PP3 | -2 | 0.93 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA2 | c.6842G&gt;T | p.(Gly2281Val) | BS3+PM2_Supporting+PP3 | -2 | 0.92 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.192T&gt;C | p.(Cys64=) | BS3+PM2_Supporting+PP3 | -2 | 0.87 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.193A&gt;C | p.(Lys65Gln) | BS3+PM2_Supporting+PP3 | -2 | 0.87 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.195G&gt;C | p.(Lys65Asn) | BS3+PM2_Supporting+PP3 | -2 | 0.87 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.195G&gt;A | p.(Lys65=) | BS3+PM2_Supporting+PP3 | -2 | 0.86 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.195G&gt;T | p.(Lys65Asn) | BS3+PM2_Supporting+PP3 | -2 | 0.86 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.194A&gt;G | p.(Lys65Arg) | BS3+PM2_Supporting+PP3 | -2 | 0.84 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA2 | c.1763A&gt;G | p.(Asn588Ser) | BS1_Supporting+PP3 | 0 | 0.83 | Curator evidence-adjudication review before any class movement. |
| C | BRCA1 | c.194A&gt;T | p.(Lys65Met) | BS3+PM2_Supporting+PP3 | -2 | 0.8 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.5073A&gt;G | p.(Thr1691=) | BS3+PM2_Supporting+PP3 | -2 | 0.71 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.5196T&gt;A | p.(His1732Gln) | BS3+PM2_Supporting+PP3 | -2 | 0.71 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.162G&gt;T | p.(Gln54His) | BS3+PM2_Supporting+PP3 | -2 | 0.7 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.187T&gt;G | p.(Leu63Val) | BS3+PM2_Supporting+PP3 | -2 | 0.7 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.5270A&gt;T | p.(Asp1757Val) | BS3+PM2_Supporting+PP3 | -2 | 0.69 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.194A&gt;C | p.(Lys65Thr) | BS3+PM2_Supporting+PP3 | -2 | 0.67 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.256C&gt;G | p.(Leu86Val) | BS3+PM2_Supporting+PP3 | -2 | 0.67 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.84G&gt;A | p.(Leu28=) | BS3+PM2_Supporting+PP3 | -2 | 0.65 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.5044G&gt;A | p.(Glu1682Lys) | BP5+PM2_Supporting+PP3+PS3 | 2 | 0.57 | Curator evidence-adjudication review before any class movement. |
| C | BRCA1 | c.83T&gt;A | p.(Leu28Gln) | BS3+PM2_Supporting+PP3 | -2 | 0.57 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.165G&gt;T | p.(Lys55Asn) | BS3+PM2_Supporting+PP3 | -2 | 0.55 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.231G&gt;T | p.(Thr77=) | BS3+PM2_Supporting+PP3 | -2 | 0.54 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.4931A&gt;T | p.(Glu1644Val) | BS3+PM2_Supporting+PP3 | -2 | 0.53 | Review functional assay details and consider RNA or independent clinical evidence. |
| C | BRCA1 | c.147G&gt;T | p.(Leu49=) | BS3+PM2_Supporting+PP3 | -2 | 0.51 | Review functional assay details and consider RNA or independent clinical evidence. |

## Tier D: Benign-Leaning BP4/BP7 Plus PM2

Why: benign splice prediction is present, but PM2 keeps the automated result in VUS. RNA or stronger benign evidence would be most useful.

| tier | gene | c_notation | p_notation | criteria_combo | total_points | spliceai_score | what_to_check |
| --- | --- | --- | --- | --- | --- | --- | --- |
| D | BRCA1 | c.222A&gt;G | p.(Gln74=) | BP4+BP7+PM2_Supporting | -1 | 0.1 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.7971A&gt;G | p.(Lys2657=) | BP4+BP7+PM2_Supporting | -1 | 0.1 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.9120T&gt;A | p.(Val3040=) | BP4+BP7+PM2_Supporting | -1 | 0.02 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8757T&gt;A | p.(Gly2919=) | BP4+BP7+PM2_Supporting | -1 | 0.02 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8757T&gt;C | p.(Gly2919=) | BP4+BP7+PM2_Supporting | -1 | 0.02 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8757T&gt;G | p.(Gly2919=) | BP4+BP7+PM2_Supporting | -1 | 0.02 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8952A&gt;C | p.(Ser2984=) | BP4+BP7+PM2_Supporting | -1 | 0.01 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8952A&gt;G | p.(Ser2984=) | BP4+BP7+PM2_Supporting | -1 | 0.01 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8952A&gt;T | p.(Ser2984=) | BP4+BP7+PM2_Supporting | -1 | 0.01 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8955T&gt;C | p.(Val2985=) | BP4+BP7+PM2_Supporting | -1 | 0.07 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8955T&gt;G | p.(Val2985=) | BP4+BP7+PM2_Supporting | -1 | 0.06 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.8955T&gt;A | p.(Val2985=) | BP4+BP7+PM2_Supporting | -1 | 0.04 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.7620G&gt;C | p.(Leu2540=) | BP4+BP7+PM2_Supporting | -1 | 0.02 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.7620G&gt;T | p.(Leu2540=) | BP4+BP7+PM2_Supporting | -1 | 0.01 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA1 | c.5079T&gt;A | p.(Ala1693=) | BP4+BP7+PM2_Supporting | -1 | 0.04 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.7804A&gt;C | p.(Arg2602=) | BP4+BP7+PM2_Supporting | -1 | 0.06 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.7618C&gt;T | p.(Leu2540=) | BP4+BP7+PM2_Supporting | -1 | 0.04 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.7806G&gt;A | p.(Arg2602=) | BP4+BP7+PM2_Supporting | -1 | 0.03 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.7803T&gt;C | p.(Tyr2601=) | BP4+BP7+PM2_Supporting | -1 | 0.02 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |
| D | BRCA2 | c.9258A&gt;C | p.(Gly3086=) | BP4+BP7+PM2_Supporting | -1 | 0.0 | Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously. |

## Tier E: PM2-Only Background

Why: PM2 alone is not enough. Most of this tier should stay background unless another independent evidence source appears.

| tier | gene | c_notation | p_notation | criteria_combo | total_points | spliceai_score | what_to_check |
| --- | --- | --- | --- | --- | --- | --- | --- |
| E | BRCA1 | c.604C&gt;T | p.(Gln202Ter) | PM2_Supporting | 1 | 0.63 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA1 | c.667A&gt;T | p.(Lys223Ter) | PM2_Supporting | 1 | 0.4 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA1 | c.548G&gt;A | p.(Gly183Glu) | PM2_Supporting | 1 | 0.19 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.318A&gt;C | p.(Gly106=) | PM2_Supporting | 1 | 0.19 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.9120T&gt;G | p.(Val3040=) | PM2_Supporting | 1 | 0.19 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.8631A&gt;T | p.(Glu2877Asp) | PM2_Supporting | 1 | 0.18 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.8752G&gt;A | p.(Glu2918Lys) | PM2_Supporting | 1 | 0.18 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.9256G&gt;A | p.(Gly3086Arg) | PM2_Supporting | 1 | 0.18 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.7436A&gt;T | p.(Asp2479Val) | PM2_Supporting | 1 | 0.17 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.7804A&gt;G | p.(Arg2602Gly) | PM2_Supporting | 1 | 0.17 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.8752G&gt;C | p.(Glu2918Gln) | PM2_Supporting | 1 | 0.17 | Do not prioritize unless another independent signal or external evidence appears. |
| E | BRCA2 | c.9120T&gt;C | p.(Val3040=) | PM2_Supporting | 1 | 0.17 | Do not prioritize unless another independent signal or external evidence appears. |

## Practical Use

Start with tiers A and B for pathogenic-side review. Use tier C as a quality-control and evidence-adjudication queue. Tier D is useful for benign resolution, especially if RNA or additional benign clinical evidence is available. Tier E should mostly stay as background.

## Outputs

- `tables/curation_review_queue/curation_review_queue.csv`
- `tables/curation_review_queue/curation_review_queue_summary.csv`
