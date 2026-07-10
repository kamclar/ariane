# Exploratory Analysis TODO

These are candidate next analyses for the generated ARIANE Module 1 BRCA1/2
coding SNV classification landscape.

Important boundary: these analyses describe patterns in the generated ARIANE
classification snapshot. They do not create independent clinical truth and they
do not add ACMG/ENIGMA criteria.

## Priority Queue

1. Benign counterexamples in pathogenic domains
   - Find benign/likely benign variants inside RING, BRCT, BRCA2 DBD, PALB2,
     RAD51/BRC, SEM1/DSS1, or other structurally/functionally important regions.
   - Purpose: protect against over-interpreting domain location as pathogenic
     evidence by itself.
   - Status: done.

2. Criterion-driver map by protein region
   - For each protein region, summarize which criteria drive class assignment:
     `BP1`, `BS3`, `BP4`, `BP7`, `PP3`, `PS3`, `PM5_PTC`, `PVS1`, `PP4`.
   - Purpose: separate biological-looking regional effects from rule behavior.

2a. BS3 and benign evidence in functional domains
   - Analyze where BS3 appears inside RING, BRCT, BRCA2 DBD, PALB2, RAD51/BRC,
     SEM1/DSS1, and other functional contexts.
   - Separate missense from synonymous variants and flag `BS3+PP3` or BS3 plus
     high SpliceAI conflicts.
   - Status: done.

3. VUS rescue map
   - For VUS in strong pathogenic regions, identify what evidence would be
     needed to move toward class 4 under existing ENIGMA/ACMG combination rules.
   - Examples: RNA evidence, functional evidence, PS1, PP1, PM3, PS4, PP4.

4. Discordant neighborhood scan
   - Find small local windows where benign and pathogenic generated variants
     occur side by side.
   - Purpose: sanity-check regions where mechanism may differ between nearby
     variants, or where automated criteria create mixed local landscapes.
   - Status: partly done as `position_class_conflict_analysis.py`, which checks
     same-CDS-position and same-codon mixed grouped classes.
   - Follow-up: focus on the smaller, more biologically interesting subset
     where the benign/pathogenic split is not simply nonsense versus
     missense/synonymous, especially non-truncating mixed codons with VUS.

## Later Ideas

5. Domain tolerance / contrast map
   - For each domain, compute benign/VUS/pathogenic ratios normalized by domain
     length and possible SNV count.

6. Pathogenic without known domain
   - Find class 4/5 variants outside mapped domains and without strong
     truncation/splice signal.
   - Purpose: identify missing annotations or criteria-driven outliers.

7. Known pathogenic residue saturation
   - For ENIGMA known pathogenic residues, inspect all possible SNVs in the same
     codon.
   - Purpose: support PS1/PM5 reasoning and detect same-amino-acid consequences.
