# Pathogenic Mechanism Exploration

Generated: 2026-06-20 14:37

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Why We Did This

This analysis looks only at generated class 4 and class 5 variants and asks what
kind of biological or rule-based mechanism appears to drive pathogenicity in the
current ARIANE Module 1 coding SNV snapshot.

It separates:

- truncating or loss-of-function-like outcomes, mostly `PVS1` and `PM5_PTC`
- functional evidence, mostly `PS3` and `PP4`
- same amino-acid evidence, `PS1`
- computational splice or protein evidence, mostly `PP3`
- location in ENIGMA-defined functional domains
- overlap with ENIGMA clinically important pathogenic residue positions

This is descriptive. It does not prove molecular mechanism for each variant.

When using this analysis in continuous text, keep two layers separate. One layer
is the analysis of ARIANE's generated classification and its criteria. The other
layer is the biological-looking pattern in the generated data landscape. The
mechanism buckets below describe what the current automated classification
appears to be driven by; they should not be presented as independently proven
biological mechanisms unless validated against external curated evidence.

## Important Limitation

ARIANE currently has ENIGMA-defined clinically important domains and selected
pathogenic residue positions. It does not yet have a full protein-function map
with all active sites, interaction interfaces, phosphorylation sites, DNA
contact residues, PALB2/RAD51/BARD1 mechanistic annotations, or structural
features. Therefore, "active site" in a strict biochemical sense is not
available from the current local data.

## Dataset

- Class 4/5 variants analyzed: 2441
- Missense class 4/5 variants: 103
- Truncating or loss-of-function-like bucket: 2326
- Variants inside an ENIGMA functional domain: 613
- Variants overlapping a known ENIGMA pathogenic residue position: 58

## Mechanism Buckets

| mechanism_bucket | count |
| --- | --- |
| truncating_pvs1_pm5_ptc | 2241 |
| truncating_pvs1_only | 85 |
| functional_plus_splice_or_computational | 61 |
| functional_multifactorial_pathogenic | 40 |
| multifactorial_pp4 | 12 |
| functional_ps3_only_or_with_pm2 | 2 |

## Aberrant Protein Buckets

| aberrant_protein_bucket | count |
| --- | --- |
| premature_stop_or_truncating_protein | 2326 |
| missense_protein_function_or_domain | 103 |
| synonymous_splice_or_functional_evidence | 9 |
| translation_initiation_effect | 3 |

## Variant Types

| variant_type | count |
| --- | --- |
| nonsense | 2326 |
| missense | 103 |
| synonymous | 9 |
| initiation_codon | 3 |

## ENIGMA Domain Context

| gene | enigma_domain | count |
| --- | --- | --- |
| BRCA2 | outside_enigma_functional_domain | 1188 |
| BRCA1 | outside_enigma_functional_domain | 640 |
| BRCA2 | DBD | 362 |
| BRCA1 | BRCT | 140 |
| BRCA1 | RING | 80 |
| BRCA1 | coiled_coil | 17 |
| BRCA2 | PALB2_binding | 14 |

## Missense Domain Context

This table is the most relevant one for questions about functional domains or
active-site-like protein effects. It excludes nonsense variants and asks where
class 4/5 missense variants fall relative to the ENIGMA domain/residue map.

| gene | enigma_domain | known_pathogenic_residue_position | count |
| --- | --- | --- | --- |
| BRCA1 | BRCT | False | 29 |
| BRCA1 | BRCT | True | 23 |
| BRCA1 | RING | False | 18 |
| BRCA2 | DBD | False | 13 |
| BRCA1 | RING | True | 12 |
| BRCA2 | DBD | True | 4 |
| BRCA1 | outside_enigma_functional_domain | False | 3 |
| BRCA2 | outside_enigma_functional_domain | False | 1 |

## Known Pathogenic Residue Context

These are positions where the generated class 4/5 variant falls at a residue
position listed in the local ENIGMA clinically important residue table. This is
position context, not proof that every alternate amino acid at that position has
the same evidence.

| gene | p_notation | known_context | count |
| --- | --- | --- | --- |
| BRCA1 | p.(Leu22Ter) | RING:Leu22Ser:c.65T&gt;C | 2 |
| BRCA1 | p.(Trp1837Ter) | BRCT:Trp1837Arg:c.5509T&gt;C | 2 |
| BRCA1 | p.(Tyr1853Ter) | BRCT:Tyr1853Cys:c.5558A&gt;G | 2 |
| BRCA2 | p.(Trp2626Ter) | DBD:Trp2626Cys:c.7878G&gt;C | 2 |
| BRCA1 | p.(Met18Thr) | RING:Met18Thr:c.53T&gt;C | 1 |
| BRCA1 | p.(Cys39Ter) | RING:Cys39Arg:c.115T&gt;C;RING:Cys39Tyr:c.116G&gt;A | 1 |
| BRCA1 | p.(His41Arg) | RING:His41Arg:c.122A&gt;G | 1 |
| BRCA1 | p.(Cys44Ser) | RING:Cys44Ser:c.130T&gt;A;RING:Cys44Tyr:c.131G&gt;A;RING:Cys44Phe:c.131G&gt;T | 1 |
| BRCA1 | p.(Cys44Phe) | RING:Cys44Ser:c.130T&gt;A;RING:Cys44Tyr:c.131G&gt;A;RING:Cys44Phe:c.131G&gt;T | 1 |
| BRCA1 | p.(Cys44Ter) | RING:Cys44Ser:c.130T&gt;A;RING:Cys44Tyr:c.131G&gt;A;RING:Cys44Phe:c.131G&gt;T | 1 |
| BRCA1 | p.(Cys44=) | RING:Cys44Ser:c.130T&gt;A;RING:Cys44Tyr:c.131G&gt;A;RING:Cys44Phe:c.131G&gt;T | 1 |
| BRCA1 | p.(Cys47Tyr) | RING:Cys47Tyr:c.140G&gt;A | 1 |
| BRCA1 | p.(Cys47Ter) | RING:Cys47Tyr:c.140G&gt;A | 1 |
| BRCA1 | p.(Cys61Gly) | RING:Cys61Gly:c.181T&gt;G;RING:Cys61Ser:c.181T&gt;A | 1 |
| BRCA1 | p.(Cys61Tyr) | RING:Cys61Gly:c.181T&gt;G;RING:Cys61Ser:c.181T&gt;A | 1 |
| BRCA1 | p.(Cys61Ser) | RING:Cys61Gly:c.181T&gt;G;RING:Cys61Ser:c.181T&gt;A | 1 |
| BRCA1 | p.(Cys61Ter) | RING:Cys61Gly:c.181T&gt;G;RING:Cys61Ser:c.181T&gt;A | 1 |
| BRCA1 | p.(Cys64Tyr) | RING:Cys64Tyr:c.191G&gt;A | 1 |
| BRCA1 | p.(Cys64Ser) | RING:Cys64Tyr:c.191G&gt;A | 1 |
| BRCA1 | p.(Cys64Phe) | RING:Cys64Tyr:c.191G&gt;A | 1 |
| BRCA1 | p.(Cys64Ter) | RING:Cys64Tyr:c.191G&gt;A | 1 |
| BRCA1 | p.(Cys64Trp) | RING:Cys64Tyr:c.191G&gt;A | 1 |
| BRCA1 | p.(Thr1685Ala) | BRCT:Thr1685Ala:c.5053A&gt;G;BRCT:Thr1685Ile:c.5054C&gt;T | 1 |
| BRCA1 | p.(Thr1685Asn) | BRCT:Thr1685Ala:c.5053A&gt;G;BRCT:Thr1685Ile:c.5054C&gt;T | 1 |
| BRCA1 | p.(Thr1685Ser) | BRCT:Thr1685Ala:c.5053A&gt;G;BRCT:Thr1685Ile:c.5054C&gt;T | 1 |
| BRCA1 | p.(Thr1685Ile) | BRCT:Thr1685Ala:c.5053A&gt;G;BRCT:Thr1685Ile:c.5054C&gt;T | 1 |
| BRCA1 | p.(Cys1697Arg) | BRCT:Cys1697Arg:c.5089T&gt;C | 1 |
| BRCA1 | p.(Cys1697Ter) | BRCT:Cys1697Arg:c.5089T&gt;C | 1 |
| BRCA1 | p.(Arg1699Gly) | BRCT:Arg1699Trp:c.5095C&gt;T | 1 |
| BRCA1 | p.(Gly1706Ter) | BRCT:Gly1706Glu:c.5117G&gt;A | 1 |

## Most Common Criteria Combinations

| criteria_combo | count |
| --- | --- |
| PM2_Supporting+PM5_PTC+PVS1 | 2010 |
| PM2_Supporting+PM5_PTC+PS3+PVS1 | 119 |
| PM5_PTC+PVS1 | 100 |
| PM2_Supporting+PVS1 | 68 |
| PM2_Supporting+PP3+PS3 | 61 |
| PM2_Supporting+PP4+PS3 | 27 |
| PM2_Supporting+PS3+PVS1 | 15 |
| PM5_PTC+PS3+PVS1 | 11 |
| PM2_Supporting+PP3+PP4+PS3 | 5 |
| PP4+PS3 | 4 |
| PP3+PP4+PS3 | 4 |
| PM2_Supporting+PP3+PP4 | 4 |
| PP3+PP4 | 4 |
| PM2_Supporting+PP4 | 4 |
| PM2_Supporting+PS1+PS3 | 2 |
| PS3+PVS1 | 2 |
| BS1_Supporting+PM5_PTC+PVS1 | 1 |

## Mechanism By Class

| predicted_class | mechanism_bucket | count |
| --- | --- | --- |
| 5 | truncating_pvs1_pm5_ptc | 2240 |
| 4 | truncating_pvs1_only | 68 |
| 4 | functional_plus_splice_or_computational | 61 |
| 5 | functional_multifactorial_pathogenic | 40 |
| 5 | truncating_pvs1_only | 17 |
| 4 | multifactorial_pp4 | 9 |
| 5 | multifactorial_pp4 | 3 |
| 4 | functional_ps3_only_or_with_pm2 | 2 |
| 4 | truncating_pvs1_pm5_ptc | 1 |

## Interpretation

- In this coding SNV snapshot, class 4/5 pathogenicity is dominated by
  truncating mechanisms: nonsense variants with `PVS1` and often `PM5_PTC`.
- Missense class 4/5 variants are much fewer and are mostly driven by curated
  functional or multifactorial evidence (`PS3`, `PP4`, `PS1`) and by location in
  ENIGMA functional domains or known pathogenic residue positions.
- SpliceAI/PP3 can contribute to pathogenic combinations, but by itself it does
  not explain most class 4/5 variants.
- We can describe "functional domain context" with current local data. We
  cannot yet describe all protein active sites or detailed structural mechanism
  without adding a curated protein feature source.

## Outputs

- `tables/pathogenic_mechanisms/pathogenic_mechanism_annotated_variants.csv`
- `tables/pathogenic_mechanisms/mechanism_summary.csv`
- `tables/pathogenic_mechanisms/aberrant_protein_summary.csv`
- `tables/pathogenic_mechanisms/variant_type_summary.csv`
- `tables/pathogenic_mechanisms/domain_summary.csv`
- `tables/pathogenic_mechanisms/gene_domain_summary.csv`
- `tables/pathogenic_mechanisms/criteria_combo_summary.csv`
- `tables/pathogenic_mechanisms/known_pathogenic_residue_context.csv`
- `tables/pathogenic_mechanisms/missense_domain_summary.csv`
- `tables/pathogenic_mechanisms/class_mechanism_summary.csv`
- `tables/pathogenic_mechanisms/gene_mechanism_summary.csv`
- `tables/pathogenic_mechanisms/missense_domain_gene_summary.csv`
- `plots/16_pathogenic_mechanisms/mechanism_summary.svg`
- `plots/16_pathogenic_mechanisms/aberrant_protein_summary.svg`
- `plots/16_pathogenic_mechanisms/gene_domain_summary.svg`
- `plots/16_pathogenic_mechanisms/gene_by_mechanism_heatmap.svg`
- `plots/16_pathogenic_mechanisms/gene_by_domain_heatmap.svg`
- `plots/16_pathogenic_mechanisms/class_by_mechanism_heatmap.svg`
- `plots/16_pathogenic_mechanisms/missense_domain_heatmap.svg`
