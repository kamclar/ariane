# Budapest 2026 Notes For ARIANE Exploration

Notes from local conference photos in:

- `F:\UOCHB\Enigma\prezentace_Budapest2026\01.jpg`
- `F:\UOCHB\Enigma\prezentace_Budapest2026\02.jpg`
- `F:\UOCHB\Enigma\prezentace_Budapest2026\03.jpg`
- `F:\UOCHB\Enigma\prezentace_Budapest2026\04.jpg`
- `F:\UOCHB\Enigma\prezentace_Budapest2026\05.jpg`

This file records what looks relevant for the exploratory analysis of the
precomputed ARIANE Module 1 classification snapshot.

## 1. gnomAD v4.1 Frequency Calibration

Conference note:

- frequency-data calibration is moving toward gnomAD v4.1.

Why this matters for ARIANE:

- The current exploratory snapshot is strongly influenced by `PM2_Supporting`.
- `PM2_Supporting` is currently one of the biggest drivers of VUS background.
- Any change in population-frequency data can change whether PM2, BS1, or BA1
  are applied.
- Therefore, gnomAD v4.1 calibration should be treated as a priority
  sensitivity analysis before changing the production classifier.

Relevant gnomAD updates:

- gnomAD v4.1 fixed an allele-number issue from v4.0 and added allele-number
  information across all callable sites.
- v4.1 reports joint exome and genome allele numbers when a variant is observed
  in either data type.
- v4.1 adds a warning flag for variants with highly discordant exome and genome
  frequencies.
- gnomAD v4.1.1 updates gene constraint metrics, LOFTEE flags, VRS IDs, gene and
  transcript quality flags, and download annotations.

Implications for exploration:

- Recompute `PM2`, `BS1`, and `BA1` using gnomAD v4.1/v4.1.1 where possible.
- Compare old versus new frequency criterion calls:
  - PM2 lost
  - PM2 gained
  - BS1 gained or lost
  - BA1 gained or lost
  - class changed because of frequency criterion change
- Track whether changes come from:
  - joint exome plus genome allele number
  - improved callable-site AN
  - exome/genome discordance flag
  - low coverage or poor mappability flags
- Keep v4.1 calibration separate from the current Module 1 snapshot until the
  rule interpretation is reviewed.

Suggested output files for a future analysis:

- `tables/frequency_v41_calibration/frequency_criterion_diff.csv`
- `tables/frequency_v41_calibration/class_transition_summary.csv`
- `tables/frequency_v41_calibration/pm2_changed_variants.csv`
- `tables/frequency_v41_calibration/bs1_ba1_changed_variants.csv`
- `plots/15_frequency_v41_calibration/class_transitions.svg`
- `plots/15_frequency_v41_calibration/pm2_change_by_gene.svg`

## 2. MaveMD And Functional Evidence

Photo `01.jpg` shows MaveMD summarizing experimental design details for a MAVE
study. The visible example is an MSH2 cell-fitness assay with clinical
performance values mapped to `BS3_Strong` and `PS3_Strong`.

Important caution from the slide:

- the example assay does not detect splicing variants
- the example assay does not detect NMD variants

Why this matters for ARIANE:

- ARIANE already uses Table 9 functional evidence as `PS3` or `BS3`.
- We have observed conflicts such as `PP3` plus `BS3` and the single
  `PVS1` plus `BS3` edge case.
- These conflicts may be biologically meaningful if the functional assay does
  not assess the same mechanism predicted by SpliceAI or PVS1.

Exploration implication:

- Future functional-evidence tables should carry assay metadata:
  - assay type
  - mechanism assessed
  - whether splicing is detected
  - whether NMD is detected
  - odds path thresholds or calibration source
- `PS3` and `BS3` should not be treated as mechanism-free evidence in conflict
  interpretation.

## 3. GA4GH VRS For Variant Identity

Photo `02.jpg` shows GA4GH VRS as a verbose but computable representation of
variants. Photo `04.jpg` shows Cat-VRS for amino-acid variants and their
relevant nucleotide variants.

Why this matters for ARIANE:

- ARIANE currently works mainly with HGVS c. and p. notation.
- Cross-source matching will become easier if variants also get stable
  computed identifiers.
- This is especially relevant for:
  - PS1 same amino-acid changes
  - mapping multiple nucleotide variants to the same protein consequence
  - matching ENIGMA, ClinVar, MaveDB, gnomAD, and local lab records

Exploration implication:

- Add future exploratory columns for normalized variant identity:
  - genomic allele
  - transcript allele
  - protein consequence
  - optional VRS identifier
  - optional Cat-VRS amino-acid grouping

This should be treated as infrastructure and matching support, not as a new
classification rule.

## 4. GA4GH VA-Spec For Evidence Representation

Photo `03.jpg` shows GA4GH VA-Spec representing granular functional evidence as:

- Variant Functional Score
- Variant Functional Class
- Variant Functional Evidence

Why this matters for ARIANE:

- This aligns well with how ARIANE needs to represent functional evidence:
  raw assay score, interpreted functional class, and final ACMG/ENIGMA evidence
  line such as `PS3_Strong` or `BS3_Strong`.
- It could become a cleaner model for importing external MAVE/functional data
  without losing the distinction between raw score and evidence strength.

Exploration implication:

- Add a future functional evidence import schema that separates:
  - study result
  - functional class statement
  - evidence line
  - assay limitation flags
  - source and calibration version

## 5. Practical Next Analyses

Recommended additions to the exploratory analysis:

1. `15_frequency_v41_calibration`
   - Compare current frequency calls against gnomAD v4.1/v4.1.1.
   - Focus on PM2, BS1, BA1 and final class transitions.

2. `16_functional_evidence_mechanism_audit`
   - Revisit `PS3/BS3` cases and annotate whether the source assay can detect
     splicing or NMD.
   - Prioritize conflicts where functional evidence and SpliceAI/PVS1 point in
     different biological directions.

3. `17_variant_identity_normalization`
   - Prototype a matching table between HGVS c., HGVS p., genomic coordinate,
     same amino-acid grouping, and optional VRS/Cat-VRS identifiers.

## Sources Checked

- gnomAD v4.1 release notes:
  https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/
- gnomAD v4.1.1 release notes:
  https://gnomad.broadinstitute.org/news/2026-03-gnomad-v4-1-1/
- GA4GH VRS documentation:
  https://vrs.ga4gh.org/en/2.0/
- GA4GH VA-Spec documentation:
  https://va-spec.ga4gh.org/en/1.0/
- MaveMD entry point:
  https://mavedb.org/mavemd
