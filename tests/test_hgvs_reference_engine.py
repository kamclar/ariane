import json
import shutil

import pytest
from fastapi.testclient import TestClient

from backend.config import PANEL_REFERENCE_DIR
from backend.models import VariantRequest
from backend.modules.hgvs import protein_notations_compatible
from backend.modules.hgvs_engine import (
    VariantNormalizationError,
    derive_protein_consequence,
    load_hgvs_engine,
)
from backend.modules.hgvs_provider import load_panel_provider
from backend.modules.variant_input import normalize_variant_input


@pytest.mark.parametrize(
    ("gene", "submitted", "canonical_c", "expected_p", "status"),
    [
        ("BRCA1", "c.303T>G", "c.303T>G", "p.(Tyr101Ter)", "sequence_derived"),
        ("BRCA1", "c.4185G>A", "c.4185G>A", "p.(Gln1395=)", "sequence_derived_synonymous"),
        ("BRCA1", "c.5266dup", "c.5266dup", "p.(Gln1756ProfsTer74)", "sequence_derived"),
        ("BRCA1", "c.2102delA", "c.2102del", "p.(Lys701SerfsTer2)", "sequence_derived"),
        ("BRCA1", "c.2102del", "c.2102del", "p.(Lys701SerfsTer2)", "sequence_derived"),
        ("BRCA1", "c.5542C>T", "c.5542C>T", "p.(Gln1848Ter)", "sequence_derived"),
        ("BRCA1", "c.1A>G", "c.1A>G", "p.(Met1?)", "sequence_derived"),
        ("BRCA1", "c.5590T>A", "c.5590T>A", "p.(Ter1864ArgextTer39)", "sequence_derived"),
        ("BRCA2", "c.9097del", "c.9097del", "p.(Thr3033LeufsTer29)", "sequence_derived"),
        ("BRCA2", "c.36dup", "c.36dup", "p.(Glu13Ter)", "sequence_derived"),
        ("BRCA1", "c.-10A>G", "c.-10A>G", "p.?", "protein_consequence_unknown"),
        ("BRCA1", "c.*10G>A", "c.*10G>A", "p.?", "protein_consequence_unknown"),
        ("BRCA1", "c.548-9A>G", "c.548-9A>G", "p.?", "protein_consequence_unknown"),
    ],
)
def test_reference_engine_derives_and_normalizes(
    gene, submitted, canonical_c, expected_p, status
):
    result = derive_protein_consequence(gene, submitted)
    assert result.canonical_c_notation == canonical_c
    assert result.p_notation == expected_p
    assert result.consequence_status == status
    assert result.provenance["normalization_engine"] == "biocommons.hgvs"
    assert result.provenance["reference_bundle"] == "ariane-brca12-reference-v1"


def test_normalization_is_idempotent():
    engine = load_hgvs_engine()
    for gene, submitted in (
        ("BRCA1", "c.2102delA"),
        ("BRCA1", "c.5266dupC"),
        ("BRCA1", "c.3668_3671dup"),
        ("BRCA2", "c.9097del"),
    ):
        first = engine.c_to_p(gene, submitted)
        second = engine.c_to_p(gene, first.canonical_c_notation)
        assert second.canonical_c_notation == first.canonical_c_notation
        assert second.p_notation == first.p_notation


def test_wrong_reference_allele_fails_closed():
    with pytest.raises(VariantNormalizationError) as raised:
        derive_protein_consequence("BRCA1", "c.303A>G")
    assert raised.value.code == "reference_allele_mismatch"


@pytest.mark.parametrize("notation", ["c.2102delC", "c.5266dupA"])
def test_wrong_indel_sequence_suffix_fails_closed(notation):
    with pytest.raises(VariantNormalizationError) as raised:
        derive_protein_consequence("BRCA1", notation)
    assert raised.value.code == "reference_allele_mismatch"


def test_random_supplied_protein_is_rejected_for_any_coding_variant():
    with pytest.raises(ValueError, match="Protein consequence mismatch"):
        normalize_variant_input("BRCA1", "c.2102del/p.(Gln1Ter)")
    with pytest.raises(ValueError, match="Protein consequence mismatch"):
        normalize_variant_input("BRCA2", "c.9097del p.(Met1Val)")


def test_legacy_abbreviated_frameshift_is_accepted_but_canonical_output_is_full():
    result = normalize_variant_input("BRCA1", "c.2102delA/p.Lys701fs")
    assert result.c_notation == "c.2102del"
    assert result.p_notation == "p.(Lys701SerfsTer2)"


def test_uncertain_exon_cnv_is_preserved_with_unknown_protein():
    notation = "c.(793+1_794-1)_(1909+1_1910-1)del"
    result = normalize_variant_input("BRCA2", notation, p_notation="p.(?)")
    assert result.c_notation == notation
    assert result.p_notation == "p.?"
    assert result.consequence_status == "protein_consequence_unknown"


def test_panel_bundle_accessions_and_independent_sequence_lengths():
    panel = load_panel_provider()
    assert panel.gene_to_transcript == {
        "BRCA1": "NM_007294.4",
        "BRCA2": "NM_000059.4",
    }
    assert panel.transcript_to_protein == {
        "NM_007294.4": "NP_009225.1",
        "NM_000059.4": "NP_000050.3",
    }
    assert panel.seqfetcher.sequence_length("NM_007294.4") == 7088
    assert panel.seqfetcher.sequence_length("NM_000059.4") == 11954
    assert panel.seqfetcher.sequence_length("NP_009225.1") == 1863
    assert panel.seqfetcher.sequence_length("NP_000050.3") == 3418


def test_missing_bundle_metadata_stops_provider(tmp_path):
    copied = tmp_path / "panel"
    shutil.copytree(PANEL_REFERENCE_DIR, copied)
    (copied / "metadata.json").unlink()
    with pytest.raises(RuntimeError, match="metadata is missing"):
        load_panel_provider(copied)


def test_corrupt_bundle_file_stops_provider(tmp_path):
    copied = tmp_path / "panel"
    shutil.copytree(PANEL_REFERENCE_DIR, copied)
    transcript_fasta = copied / "fasta" / "transcripts.fa"
    transcript_fasta.write_text(
        transcript_fasta.read_text(encoding="ascii") + "A", encoding="ascii"
    )
    with pytest.raises(RuntimeError, match="checksum mismatch"):
        load_panel_provider(copied)


def test_manifest_does_not_contain_precomputed_protein_consequences():
    manifest = json.loads((PANEL_REFERENCE_DIR / "panel_manifest.json").read_text())
    assert all("p_notation" not in record for record in manifest["transcripts"])


def test_engine_and_snapshot_notation_compatibility_helper_handles_abbreviation():
    assert protein_notations_compatible(
        "p.(Gln1756fs)", "p.(Gln1756ProfsTer74)"
    )


@pytest.mark.parametrize(
    ("c_notation", "expected_p"),
    [
        ("c.1A>G", "p.(Met1?)"),
        ("c.5590T>A", "p.(Ter1864ArgextTer39)"),
    ],
)
def test_request_model_accepts_canonical_start_and_stop_extension_forms(
    c_notation, expected_p
):
    request = VariantRequest(gene="BRCA1", c_notation=c_notation)
    assert request.p_notation == expected_p


def test_normalize_api_uses_same_engine_and_returns_provenance():
    from backend.main import app

    response = TestClient(app).post(
        "/api/normalize", json={"gene": "BRCA1", "c_notation": "c.2102delA"}
    )
    assert response.status_code == 200
    result = response.json()
    assert result["c_notation"] == "c.2102del"
    assert result["p_notation"] == "p.(Lys701SerfsTer2)"
    assert result["consequence_status"] == "sequence_derived"
    assert result["normalization_provenance"]["reference_bundle"] == (
        "ariane-brca12-reference-v1"
    )


def test_normalize_api_rejects_random_protein_consequence():
    from backend.main import app

    response = TestClient(app).post(
        "/api/normalize",
        json={
            "gene": "BRCA1",
            "c_notation": "c.2102del",
            "p_notation": "p.(Gln1Ter)",
        },
    )
    assert response.status_code == 422
    assert "Protein consequence mismatch" in response.text
