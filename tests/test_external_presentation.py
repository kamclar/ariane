from types import SimpleNamespace

from backend.services.classification_presentation import _external_model


def evidence(*, clinvar, clingen):
    return SimpleNamespace(
        clinvar=clinvar,
        clingen=clingen,
        result={"predicted_class": 3},
        variant=SimpleNamespace(gene="BRCA1", c_notation="c.5217T>A"),
    )


def test_clingen_result_remains_visible_when_clinvar_has_no_record():
    result = _external_model(evidence(
        clinvar={"status": "not_found"},
        clingen={
            "status": "ok",
            "classification": "Likely pathogenic",
            "evidence_codes": [{"code": "PS3", "status": "Met"}],
        },
    ))

    assert result.clinvar_status == "not_found"
    assert result.enigma_ep_class == "Likely pathogenic"
    assert result.enigma_ep_source == "ClinGen ERepo"
    assert result.erepo_evidence_codes == ["PS3"]


def test_external_section_has_explicit_status_when_both_services_fail():
    result = _external_model(evidence(
        clinvar={"status": "api_timeout", "error": "ClinVar timeout"},
        clingen={"status": "api_error", "error": "ERepo unavailable"},
    ))

    assert result is not None
    assert result.clinvar_status == "api_timeout"
    assert result.clingen_status == "api_error"
    assert "unavailable" in result.clinvar_message
    assert result.clinvar_error == "ClinVar timeout"
    assert result.clingen_error == "ERepo unavailable"


def test_ambiguous_clinvar_candidates_are_reported_without_selecting_one():
    result = _external_model(evidence(
        clinvar={
            "status": "ambiguous",
            "candidate_ids": ["123", "456"],
        },
        clingen={"status": "not_found"},
    ))

    assert result.clinvar_status == "ambiguous"
    assert "123, 456" in result.clinvar_message
    assert not result.clinvar_classification
