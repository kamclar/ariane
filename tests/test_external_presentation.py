from types import SimpleNamespace
from unittest.mock import patch

from backend.contracts import ClassificationResult
from backend.services.classification_presentation import (
    ClassificationPresentationService,
    _criterion_models,
    _evidence_display_flags,
    _external_model,
)


def evidence(*, clinvar, clingen, predicted_class=3, c_notation="c.5217T>A"):
    return SimpleNamespace(
        clinvar=clinvar,
        clingen=clingen,
        result={"predicted_class": predicted_class},
        variant=SimpleNamespace(gene="BRCA1", c_notation=c_notation),
    )


def test_cached_result_gets_current_external_data_and_replaces_stale_warnings():
    cached = ClassificationResult(
        variant="BRCA1 c.5217T>A p.(Val1739=)",
        gene="BRCA1",
        c_notation="c.5217T>A",
        p_notation="p.(Val1739=)",
        predicted_class=3,
        warnings=[
            "Keep this classification warning.",
            "ClinVar comparison is temporarily unavailable.",
        ],
    )

    refreshed = ClassificationPresentationService().refresh_external(
        cached,
        clinvar={"status": "not_found"},
        clingen={"status": "not_found"},
    )

    assert refreshed.external is not None
    assert refreshed.external.clinvar_status == "not_found"
    assert refreshed.warnings == ["Keep this classification warning."]


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
    assert "assessed variant" in result.clinvar_message
    assert "123" not in result.clinvar_message
    assert result.clinvar_candidate_ids == ["123", "456"]
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


def test_rna_evidence_is_not_presented_as_table9_functional_evidence():
    criteria = _criterion_models(
        {
            "PVS1_RNA": {
                "strength": "Strong",
                "points": 4,
                "applies": True,
                "source": "ENIGMA Supplementary Tables 2 and 3",
            }
        },
        applies=True,
    )

    has_table9, has_rna = _evidence_display_flags(criteria)

    assert has_table9 is False
    assert has_rna is True


def test_table9_ps3_is_presented_as_functional_not_rna_evidence():
    criteria = _criterion_models(
        {
            "PS3": {
                "strength": "Strong",
                "points": 4,
                "applies": True,
                "table9_audit": {"source_row": 17},
            }
        },
        applies=True,
    )

    has_table9, has_rna = _evidence_display_flags(criteria)

    assert has_table9 is True
    assert has_rna is False


def test_historical_multifactorial_enigma_difference_is_prominent_and_explicit():
    with patch(
        "backend.services.classification_presentation.lookup_erepo_vcep_assertion",
        return_value={"status": "not_found", "record": None},
    ):
        result = _external_model(
            evidence(
                clinvar={
                    "status": "ok",
                    "aggregate": {
                        "classification": "Pathogenic",
                        "review_status": "reviewed by expert panel",
                        "n_submitters": 3,
                    },
                    "enigma_submission": {
                        "class": "Pathogenic",
                        "date_eval": "2019-06-18",
                        "comment": (
                            "IARC class based on posterior probability from "
                            "multifactorial likelihood analysis"
                        ),
                    },
                    "submissions": [
                        {
                            "scv": "SCV001161546",
                            "org": "ENIGMA",
                            "class": "Pathogenic",
                            "date_eval": "2019-06-18",
                            "is_enigma_ep": True,
                            "review": "reviewed by expert panel",
                        }
                    ],
                },
                clingen={"status": "not_found"},
                predicted_class=4,
                c_notation="c.4185G>A",
            )
        )

    assert result.expert_panel_difference_message == (
        "Historical ENIGMA expert-panel classification: Pathogenic; current "
        "ARIANE v1.2 automated result: Likely Pathogenic."
    )
    assert "multifactorial posterior probability" in (
        result.historical_expert_panel_warning
    )
    assert "2019-06-18" in result.historical_expert_panel_warning
    assert result.clinvar_submitters[0].curated_status == (
        "ENIGMA expert-panel ClinVar assertion"
    )
