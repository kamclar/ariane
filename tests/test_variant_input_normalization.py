import pytest
from pydantic import ValidationError

from backend.models import VariantRequest
from backend.modules.variant_input import normalize_variant_input


@pytest.mark.parametrize(
    ("gene", "notation", "expected_c", "expected_p", "transcript"),
    [
        ("BRCA1", "c.303T>G", "c.303T>G", "p.(Tyr101Ter)", "NM_007294.4"),
        ("BRCA1", "NM_007294.4:c.509G>A", "c.509G>A", "p.(Arg170Gln)", "NM_007294.4"),
        ("BRCA1", "c.5266dup", "c.5266dup", "p.(Gln1756ProfsTer74)", "NM_007294.4"),
        ("BRCA1", "c.548-9A>G", "c.548-9A>G", "p.(?)", "NM_007294.4"),
        ("BRCA2", "c.36dup", "c.36dup", "p.(Glu13Ter)", "NM_000059.4"),
    ],
)
def test_reference_transcript_input_derives_canonical_protein(
    gene, notation, expected_c, expected_p, transcript
):
    result = normalize_variant_input(gene, notation)
    assert result.c_notation == expected_c
    assert result.p_notation == expected_p
    assert result.reference_transcript == transcript
    assert result.normalization_source


@pytest.mark.parametrize(
    ("gene", "notation", "expected_c", "expected_p"),
    [
        ("BRCA2", "c.3703c>t", "c.3703C>T", "p.(Gln1235Ter)"),
        (
            "BRCA1",
            "nm_007294.4 : c.3668_3671DUP / p.cys1225FS",
            "c.3668_3671dup",
            "p.(Cys1225SerfsTer10)",
        ),
        (
            "BRCA1",
            "NM_007294.4\t:\tC.3668_3671dup\tP.(Cys1225fs)",
            "c.3668_3671dup",
            "p.(Cys1225SerfsTer10)",
        ),
    ],
)
def test_hgvs_input_accepts_case_slash_and_flexible_separator_whitespace(
    gene, notation, expected_c, expected_p
):
    result = normalize_variant_input(gene, notation)
    assert result.c_notation == expected_c
    assert result.p_notation == expected_p


@pytest.mark.parametrize(
    ("assembly", "notation"),
    [
        ("GRCh37", "17:41251830:C>T"),
        ("GRCh38", "chr17:43099813:C>T"),
        ("hg38", "17-43099813-C-T"),
    ],
)
def test_genomic_input_maps_to_reference_transcript_hgvs(assembly, notation):
    result = normalize_variant_input("BRCA1", notation, assembly=assembly)
    assert result.c_notation == "c.509G>A"
    assert result.p_notation == "p.(Arg170Gln)"
    assert result.reference_transcript == "NM_007294.4"


def test_genomic_input_requires_assembly():
    with pytest.raises(ValueError, match="assembly is required"):
        normalize_variant_input("BRCA1", "chr17:43099813:C>T")


def test_unknown_intronic_protein_consequence_is_explained():
    result = normalize_variant_input("BRCA1", "c.548-9A>G")
    assert result.p_notation == "p.(?)"
    assert "RNA evidence" in result.protein_consequence_explanation
    assert "DNA notation alone" in result.protein_consequence_explanation


def test_known_protein_consequence_has_no_unknown_explanation():
    result = normalize_variant_input("BRCA1", "c.303T>G")
    assert result.p_notation == "p.(Tyr101Ter)"
    assert result.protein_consequence_explanation == ""


def test_wrong_reference_allele_is_rejected_before_classification():
    with pytest.raises(ValueError, match="Reference allele mismatch"):
        normalize_variant_input("BRCA1", "c.181A>C", p_notation="p.(Cys61Gly)")


def test_contradictory_protein_description_is_rejected():
    with pytest.raises(ValueError, match="Protein consequence mismatch"):
        normalize_variant_input("BRCA1", "c.303T>G", p_notation="p.(Arg170Gln)")


@pytest.mark.parametrize(
    ("gene", "c_notation", "abbreviated_p", "canonical_p"),
    [
        (
            "BRCA1", "c.3668_3671dup", "p.(Cys1225fs)",
            "p.(Cys1225SerfsTer10)",
        ),
        (
            "BRCA2", "c.9097del", "p.(Thr3033fs)",
            "p.(Thr3033LeufsTer29)",
        ),
    ],
)
def test_abbreviated_frameshift_is_accepted_but_output_is_canonical(
    gene, c_notation, abbreviated_p, canonical_p
):
    result = normalize_variant_input(gene, c_notation, p_notation=abbreviated_p)
    assert result.p_notation == canonical_p


@pytest.mark.parametrize("wrong_p", ["p.(Arg1225fs)", "p.(Cys1226fs)"])
def test_abbreviated_frameshift_must_match_reference_amino_acid_and_position(wrong_p):
    with pytest.raises(ValueError, match="Protein consequence mismatch"):
        normalize_variant_input("BRCA1", "c.3668_3671dup", p_notation=wrong_p)


def test_wrong_transcript_is_rejected():
    with pytest.raises(ValueError, match="does not match BRCA1"):
        normalize_variant_input("BRCA1", "NM_000059.4:c.303T>G")


def test_unmapped_genomic_input_fails_closed():
    with pytest.raises(ValueError, match="No validated local BRCA1 mapping"):
        normalize_variant_input(
            "BRCA1", "chr17:43000000:C>G", assembly="GRCh38"
        )


def test_variant_request_retains_submitted_and_canonical_descriptions():
    request = VariantRequest(
        gene="BRCA1",
        c_notation="chr17:43099813:C>T",
        assembly="GRCh38",
    )
    assert request.submitted_notation == "chr17:43099813:C>T"
    assert request.c_notation == "c.509G>A"
    assert request.p_notation == "p.(Arg170Gln)"
    assert request.reference_transcript == "NM_007294.4"


def test_variant_request_reports_missing_assembly_as_validation_error():
    with pytest.raises(ValidationError, match="assembly is required"):
        VariantRequest(gene="BRCA1", c_notation="chr17:43099813:C>T")
