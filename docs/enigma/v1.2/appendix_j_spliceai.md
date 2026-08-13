# Appendix J: SpliceAI implementation checkpoint

This note is an implementation index, not a replacement for the official
source in `source/Appendix_V1.2.docx`.

ARIANE must use the ENIGMA v1.2 Appendix J profile:

- maximum of `DS_AG`, `DS_AL`, `DS_DG`, and `DS_DL`
- maximum distance 10,000 bases
- unmasked scores (`mask=0`)
- BP4 at a maximum delta score less than or equal to 0.10
- no computational criterion above 0.10 and below 0.20
- PP3 at a maximum delta score greater than or equal to 0.20
- REF and ALT component scores retained and displayed for audit

The machine-readable profile is
`data/spliceai/enigma_v1_2_spliceai_profile.json`. Classification code and
cache builders use that same file. Cache metadata with a different profile,
distance, masking mode, assembly, transcript policy, aggregation, or score
field set is rejected. It is never reinterpreted as ENIGMA-compatible data.
