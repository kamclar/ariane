from types import SimpleNamespace
from unittest.mock import patch

from backend.services.classification_presentation import _criterion_models, _external_model


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
            "guideline_versions": ["1.2.0"],
            "cspec_ids": ["GN092"],
            "assertion_id": "CG:test",
        },
    ))

    assert result.clinvar_status == "not_found"
    assert result.enigma_ep_class == "Likely pathogenic"
    assert result.enigma_ep_source == "Live ClinGen ERepo comparison"
    assert result.erepo_evidence_codes == ["PS3"]
    assert result.erepo_guideline_versions == ["1.2.0"]
    assert result.erepo_cspec_ids == ["GN092"]
    assert result.erepo_assertion_id == "CG:test"


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


def test_historical_erepo_assertion_is_labelled_as_context_only():
    with patch(
        "backend.services.classification_presentation.lookup_erepo_vcep_assertion",
        return_value={
            "status": "historical_vcep_assertion",
            "record": {
                "classification": "Pathogenic",
                "assertion_method_version": "1.0.0",
                "uuid": "11111111-2222-3333-4444-555555555555",
            },
        },
    ):
        result = _external_model(evidence(
            clinvar={"status": "not_found"},
            clingen={"status": "not_found"},
        ))

    assert result.enigma_ep_class == "Pathogenic"
    assert "context only" in result.enigma_ep_source
    assert "not the active v1.2" in result.historical_expert_panel_warning


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


def test_table9_audit_is_preserved_in_public_criterion_model():
    audit = {
        "gene": "BRCA1",
        "c_notation": "c.5074G>A",
        "standardised_text": "Complete reviewed Table 9 text.",
        "assay_results": ["Protein and mRNA assay result."],
    }
    criteria = _criterion_models(
        {
            "PS3": {
                "code": "PS3",
                "strength": "Strong",
                "points": 4,
                "applies": True,
                "table9_audit": audit,
            }
        },
        applies=True,
    )

    assert len(criteria) == 1
    assert criteria[0].table9_audit == audit
