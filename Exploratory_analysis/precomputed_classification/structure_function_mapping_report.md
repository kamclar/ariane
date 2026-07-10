# Structure Function Mapping

Generated: 2026-06-20 16:19

Input: `variant_space_scan\outputs\brca_module1_full_snv_classification.csv`

## Why We Did This

This analysis focuses on generated class 4/5 variants that are pathogenic
without an obvious aberrant-protein or splice driver in the current snapshot.
It excludes:

- nonsense/frameshift/splice-site variants
- variants with `PVS1` or `PM5_PTC`
- variants with SpliceAI >= 0.20
- variants with `PP3`

The goal is to ask: among the remaining pathogenic variants, do they map to
protein regions where a missense or initiation-codon change could plausibly
affect protein function?

## Dataset

- Non-truncating, no-strong-splice class 4/5 variants: 37
- Missense variants in this subset: 34
- Variants at a known ENIGMA pathogenic residue position: 28

## Important Interpretation Boundary

This is not full 3D structural modeling. It maps variants to broad
structure/function features using local ENIGMA domains plus a small curated
feature scaffold from public protein annotations. A real 3D analysis would need
explicit PDB or AlphaFold residue coordinates, structure confidence, partner
interfaces, ligand/DNA contacts, and probably separate models for BRCA1-BARD1
and BRCA2-DSS1/RAD51/PALB2 complexes.

## Feature Summary

| feature | count |
| --- | --- |
| BRCT phosphopeptide-binding region | 20 |
| DNA-binding domain / helical-OB-DSS1 region | 8 |
| RING zinc-binding/E3 ligase region | 6 |
| outside_mapped_structure_function_features | 3 |

## Gene Feature Summary

| gene | feature | count |
| --- | --- | --- |
| BRCA1 | BRCT phosphopeptide-binding region | 20 |
| BRCA2 | DNA-binding domain / helical-OB-DSS1 region | 8 |
| BRCA1 | RING zinc-binding/E3 ligase region | 6 |
| BRCA1 | outside_mapped_structure_function_features | 2 |
| BRCA2 | outside_mapped_structure_function_features | 1 |

## Criteria Driving This Subset

| criteria_combo | count |
| --- | --- |
| PM2_Supporting+PP4+PS3 | 27 |
| PP4+PS3 | 4 |
| PM2_Supporting+PP4 | 4 |
| PM2_Supporting+PS1+PS3 | 2 |

## Variant Types

| variant_type | count |
| --- | --- |
| missense | 34 |
| initiation_codon | 3 |

## Known Pathogenic Residue Context

| known_pathogenic_residue_context | count |
| --- | --- |
| known_pathogenic_residue_position | 28 |
| no_known_pathogenic_residue_position | 9 |

## Curated Active Site Or Interface Status

The analysis now reads residue-level active-site/interface annotations from
`tables/structure_function_mapping/curated_active_site_interface_annotations.csv`.
It also reads UniProt region-level interface annotations from
`tables/structure_function_mapping/uniprot_brca1_interface_regions.csv` and
`tables/structure_function_mapping/uniprot_brca2_interface_regions.csv`.
This is deliberately separate from broad domain context. A variant can sit in a
BRCT or DBD domain without being annotated as an exact binding/interface residue.
Rows prefixed with `region_level` should be interpreted as regional context, not
as proof that the exact residue directly contacts a partner.

| curated_active_site_status | count |
| --- | --- |
| region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | 8 |
| not_curated_in_current_dataset | 7 |
| region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 6 |
| region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | 5 |
| zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys44 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 2 |
| UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | 2 |
| UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | 2 |
| zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue His41 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 1 |
| zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys47 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 1 |
| UniProt:Natural variant:in BC and OC; pathogenic; no interaction with BAP1; dbSNP:rs28897672 [source_extracted];zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys61 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | 1 |
| UniProt:Natural variant:in BC and FANCS; pathogenic; decreased localization to DNA damage sites and reduced interaction with UIMC1/RAP80; dbSNP:rs45553935 [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | 1 |
| region_level:UniProt:Region:Interaction with PALB2 (1-40) [source_extracted_region] | 1 |

## Top Variants

| gene | c_notation | p_notation | variant_type | predicted_class | structure_features | structure_domain_context | curated_active_site_status | known_pathogenic_residue_context | criteria_combo |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BRCA1 | c.1A&gt;G | p.(Met1Val) | initiation_codon | 5 | outside_mapped_structure_function_features | outside_mapped_structure_function_features | not_curated_in_current_dataset |  | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.3G&gt;T | p.(Met1Ile) | initiation_codon | 5 | outside_mapped_structure_function_features | outside_mapped_structure_function_features | not_curated_in_current_dataset |  | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.53T&gt;C | p.(Met18Thr) | missense | 5 | RING zinc-binding/E3 ligase region | RING zinc-binding/E3 ligase region | not_curated_in_current_dataset | RING:Met18Thr:c.53T&gt;C | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.122A&gt;G | p.(His41Arg) | missense | 5 | RING zinc-binding/E3 ligase region | RING zinc-binding/E3 ligase region | zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue His41 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | RING:His41Arg:c.122A&gt;G | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.130T&gt;A | p.(Cys44Ser) | missense | 5 | RING zinc-binding/E3 ligase region | RING zinc-binding/E3 ligase region | zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys44 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | RING:Cys44Ser:c.130T&gt;A;RING:Cys44Tyr:c.131G&gt;A;RING:Cys44Phe:c.131G&gt;T | PP4+PS3 |
| BRCA1 | c.131G&gt;T | p.(Cys44Phe) | missense | 5 | RING zinc-binding/E3 ligase region | RING zinc-binding/E3 ligase region | zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys44 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | RING:Cys44Ser:c.130T&gt;A;RING:Cys44Tyr:c.131G&gt;A;RING:Cys44Phe:c.131G&gt;T | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.140G&gt;A | p.(Cys47Tyr) | missense | 5 | RING zinc-binding/E3 ligase region | RING zinc-binding/E3 ligase region | zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys47 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | RING:Cys47Tyr:c.140G&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.182G&gt;C | p.(Cys61Ser) | missense | 4 | RING zinc-binding/E3 ligase region | RING zinc-binding/E3 ligase region | UniProt:Natural variant:in BC and OC; pathogenic; no interaction with BAP1; dbSNP:rs28897672 [source_extracted];zinc_coordination_core:BRCA1 RING C3HC4 zinc-coordinating residue Cys61 [seed_needs_source_line_review];region_level:UniProt:Zinc finger:RING-type (24-65) [source_extracted_region] | RING:Cys61Gly:c.181T&gt;G;RING:Cys61Ser:c.181T&gt;A | PM2_Supporting+PS1+PS3 |
| BRCA1 | c.5053A&gt;G | p.(Thr1685Ala) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BRCT:Thr1685Ala:c.5053A&gt;G;BRCT:Thr1685Ile:c.5054C&gt;T | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5054C&gt;T | p.(Thr1685Ile) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BRCT:Thr1685Ala:c.5053A&gt;G;BRCT:Thr1685Ile:c.5054C&gt;T | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5089T&gt;C | p.(Cys1697Arg) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BRCT:Cys1697Arg:c.5089T&gt;C | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5117G&gt;A | p.(Gly1706Glu) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BRCT:Gly1706Glu:c.5117G&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5143A&gt;C | p.(Ser1715Arg) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BRCT:Ser1715Arg:c.5143A&gt;C;BRCT:Ser1715Asn:c.5144G&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5144G&gt;A | p.(Ser1715Asn) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BRCT:Ser1715Arg:c.5143A&gt;C;BRCT:Ser1715Asn:c.5144G&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5207T&gt;C | p.(Val1736Ala) | missense | 4 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC and FANCS; pathogenic; decreased localization to DNA damage sites and reduced interaction with UIMC1/RAP80; dbSNP:rs45553935 [source_extracted];region_level:UniProt:Domain:BRCT 1 (1642-1736) [source_extracted_region] | BRCT:Val1736Ala:c.5207T&gt;C;BRCT:Val1736Gly:c.5207T&gt;G | PM2_Supporting+PP4 |
| BRCA1 | c.5212G&gt;A | p.(Gly1738Arg) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | BRCT:Gly1738Arg:c.5212G&gt;A;BRCT:Gly1738Glu:c.5213G&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5213G&gt;A | p.(Gly1738Glu) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | UniProt:Mutagenesis:Abolishes interaction with BRIP1. [source_extracted] | BRCT:Gly1738Arg:c.5212G&gt;A;BRCT:Gly1738Glu:c.5213G&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5243G&gt;A | p.(Gly1748Asp) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | not_curated_in_current_dataset | BRCT:Gly1748Asp:c.5243G&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5291T&gt;C | p.(Leu1764Pro) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Leu1764Pro:c.5291T&gt;C | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5297T&gt;G | p.(Ile1766Ser) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Ile1766Ser:c.5297T&gt;G | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5309G&gt;T | p.(Gly1770Val) | missense | 4 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Gly1770Val:c.5309G&gt;T | PM2_Supporting+PP4 |
| BRCA1 | c.5324T&gt;A | p.(Met1775Lys) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Met1775Lys:c.5324T&gt;A;BRCT:Met1775Arg:c.5324T&gt;G | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5324T&gt;G | p.(Met1775Arg) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | UniProt:Natural variant:in BC; pathogenic; strongly reduced transcription transactivation; abolishes interaction with BRIP1 and RBBP8; dbSNP:rs41293463 [source_extracted];UniProt:Natural variant:in BC; pathogenic; alters protein stability and abolishes ACACA and BRIP1 binding; dbSNP:rs41293463 [source_extracted];region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Met1775Lys:c.5324T&gt;A;BRCT:Met1775Arg:c.5324T&gt;G | PP4+PS3 |
| BRCA1 | c.5359T&gt;A | p.(Cys1787Ser) | missense | 4 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] |  | PM2_Supporting+PP4 |
| BRCA1 | c.5509T&gt;C | p.(Trp1837Arg) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Trp1837Arg:c.5509T&gt;C | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5513T&gt;A | p.(Val1838Glu) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Val1838Glu:c.5513T&gt;A | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5516T&gt;C | p.(Leu1839Ser) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Leu1839Ser:c.5516T&gt;C | PM2_Supporting+PP4+PS3 |
| BRCA1 | c.5558A&gt;G | p.(Tyr1853Cys) | missense | 5 | BRCT phosphopeptide-binding region | BRCT phosphopeptide-binding region | region_level:UniProt:Domain:BRCT 2 (1756-1855) [source_extracted_region] | BRCT:Tyr1853Cys:c.5558A&gt;G | PM2_Supporting+PP4+PS3 |
| BRCA2 | c.3G&gt;A | p.(Met1Ile) | initiation_codon | 5 | outside_mapped_structure_function_features | outside_mapped_structure_function_features | region_level:UniProt:Region:Interaction with PALB2 (1-40) [source_extracted_region] |  | PM2_Supporting+PP4+PS3 |
| BRCA2 | c.7879A&gt;T | p.(Ile2627Phe) | missense | 5 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | DBD:Ile2627Phe:c.7879A&gt;T | PM2_Supporting+PP4+PS3 |
| BRCA2 | c.7958T&gt;C | p.(Leu2653Pro) | missense | 5 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] |  | PM2_Supporting+PP4+PS3 |
| BRCA2 | c.7988A&gt;T | p.(Glu2663Val) | missense | 4 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] | DBD:Glu2663Val:c.7988A&gt;T | PM2_Supporting+PP4 |
| BRCA2 | c.8009C&gt;T | p.(Ser2670Leu) | missense | 5 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] |  | PM2_Supporting+PP4+PS3 |
| BRCA2 | c.8165C&gt;G | p.(Thr2722Arg) | missense | 5 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | region_level:UniProt:Region:Interaction with SEM1 (2481-2832) [source_extracted_region] |  | PM2_Supporting+PP4+PS3 |
| BRCA2 | c.9285C&gt;A | p.(Asp3095Glu) | missense | 4 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | not_curated_in_current_dataset |  | PM2_Supporting+PS1+PS3 |
| BRCA2 | c.9285C&gt;G | p.(Asp3095Glu) | missense | 5 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | not_curated_in_current_dataset |  | PP4+PS3 |
| BRCA2 | c.9371A&gt;T | p.(Asn3124Ile) | missense | 5 | DNA-binding domain / helical-OB-DSS1 region | DNA-binding domain / helical-OB-DSS1 region | not_curated_in_current_dataset | DBD:Asn3124Ile:c.9371A&gt;T | PP4+PS3 |

## Interpretation

- The class 4/5 variants without truncation and without strong splice/PP3 signal
  form a small subset.
- Most of this subset is driven by curated functional/multifactorial evidence,
  especially `PS3` and `PP4`, rather than by rules inferred from structure.
- The strongest structure/function concentration is in BRCA1 RING and BRCT
  regions, with some BRCA2 DNA-binding-domain variants.
- For truly structural interpretation, the next step is to add residue-level
  coordinates from selected experimental PDB structures and/or AlphaFold models
  and annotate partner interfaces.

## Active Site Annotation Curation

The first residue-level table was generated from downloaded UniProt BRCA1/BRCA2
snapshots plus a source-tracked seed for the BRCA1 RING zinc-coordinating core.
UniProt provides strong domain and interaction context, but not a complete
ready-made residue-level interface table for all BRCA1/2 functional surfaces.
Some BRCT and BRCA2 positions are now extracted from UniProt `Mutagenesis` or
`Natural variant` records with interaction effects. Remaining BRCT
phosphopeptide-binding surface details and BRCA2 DNA/DSS1/RAD51/PALB2 contact
residues still need explicit PDB/literature curation before we should treat them
as exact active-site/interface annotations.

Minimum useful fields for that table:

- `gene`
- `protein_position`
- `annotation_type`
- `annotation_label`
- `evidence_source`
- `source_id`
- `structure_id`
- `notes`

Candidate sources to curate from:

- UniProt feature annotations for BRCA1 and BRCA2
- selected experimental PDB structures
- AlphaFold models as spatial context, with confidence filtering
- InterPro/Pfam domain annotations
- literature for BRCA1-BARD1 RING, BRCA1 BRCT phosphopeptide binding, BRCA2
  DNA/DSS1/RAD51/PALB2 interfaces

## Suggested Next 3D Step

Start with the small subset in
`tables/structure_function_mapping/nontruncating_no_splice_pathogenic_variants.csv`.
For each variant, map the amino-acid position to a structure source:

- BRCA1 RING: experimental BRCA1-BARD1 N-terminal structures
- BRCA1 BRCT: experimental BRCT phosphopeptide-binding structures
- BRCA2 DBD: experimental DNA-binding/DSS1 region structures

Then annotate distance to zinc-coordinating residues, phosphopeptide-binding
surface, DNA/DSS1/RAD51/PALB2 interface, or other curated functional surfaces.

## Sources

- UniProt BRCA1 P38398: https://rest.uniprot.org/uniprotkb/P38398.txt
- UniProt BRCA2 P51587: https://rest.uniprot.org/uniprotkb/P51587.txt
- AlphaFold BRCA1 P38398: https://alphafold.ebi.ac.uk/entry/P38398
- AlphaFold BRCA2 P51587: https://alphafold.ebi.ac.uk/entry/P51587

## Outputs

- `tables/structure_function_mapping/nontruncating_no_splice_pathogenic_variants.csv`
- `tables/structure_function_mapping/feature_summary.csv`
- `tables/structure_function_mapping/gene_feature_summary.csv`
- `tables/structure_function_mapping/criteria_combo_summary.csv`
- `tables/structure_function_mapping/variant_type_summary.csv`
- `tables/structure_function_mapping/known_residue_summary.csv`
- `tables/structure_function_mapping/curated_active_site_status_summary.csv`
- `tables/structure_function_mapping/curated_active_site_interface_annotations.csv`
- `tables/structure_function_mapping/uniprot_brca1_interface_regions.csv`
- `tables/structure_function_mapping/uniprot_brca2_interface_regions.csv`
- `plots/17_structure_function_mapping/feature_summary.svg`
- `plots/17_structure_function_mapping/gene_feature_summary.svg`
- `plots/17_structure_function_mapping/curated_active_site_status_summary.svg`
- `plots/17_structure_function_mapping/brca1_lollipop.svg`
- `plots/17_structure_function_mapping/brca2_lollipop.svg`
