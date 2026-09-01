"""Parity and fail-closed tests for the review-adjusted classification DAG."""

import pytest
from pydantic import ValidationError

from backend.classification_dag.manual import execute_manual_evidence
from backend.classification_dag.engine import DagNodeExecutionError
from backend.classification_dag.runtime import compare_classification_results
from backend.models import ManualCriterionInput, ManualEvidenceRequest
from backend.modules.manual_evidence import evaluate_manual_evidence


PVS1_RNA = [{
    "code": "PVS1_RNA",
    "enabled": True,
    "evidence": {
        "assay_scope": "mrna_only",
        "rna_conclusion": "damaging",
        "functional_transcript_remaining": "absent_or_minimal",
        "curated_strength": "Very Strong",
        "transcript_accession": "NM_007294.4",
        "tissue_or_cell_type": "blood",
        "nmd_assessed": "yes",
    },
}]


def _manual_request_item(**overrides):
    item = {
        "code": "PP1",
        "enabled": True,
        "evidence": {"likelihood_ratio": 2.08},
        "notes": "Quantitative segregation evidence reviewed.",
        "references": ["PMID:1"],
    }
    item.update(overrides)
    return item


def test_manual_request_backend_requires_an_enabled_criterion():
    with pytest.raises(ValidationError, match="Select at least one"):
        ManualEvidenceRequest(
            base_criteria=[],
            manual_criteria=[_manual_request_item(enabled=False)],
            assessor="reviewer",
            assessed_at="2026-08-24",
        )


@pytest.mark.parametrize(
    ("override", "message"),
    (
        ({"notes": ""}, "requires evidence notes"),
        ({"references": []}, "requires at least one evidence reference"),
        ({"references": ["  "]}, "requires at least one evidence reference"),
    ),
)
def test_manual_request_backend_requires_auditable_enabled_records(override, message):
    with pytest.raises(ValidationError, match=message):
        ManualEvidenceRequest(
            base_criteria=[],
            manual_criteria=[_manual_request_item(**override)],
            assessor="reviewer",
            assessed_at="2026-08-24",
        )


def test_manual_request_backend_accepts_complete_auditable_record():
    request = ManualEvidenceRequest(
        base_criteria=[],
        manual_criteria=[_manual_request_item()],
        assessor="reviewer",
        assessed_at="2026-08-24",
    )
    assert request.manual_criteria[0].enabled is True


@pytest.mark.parametrize(
    ("base", "manual", "variant_context"),
    (
        (
            [{"name": "PS3", "applies": True, "strength": "Strong", "points": 4, "reason": "functional"}],
            [
                {
                    "code": "PM3",
                    "enabled": True,
                    "evidence": {
                        "evidence_points": 1,
                        "cooccurring_variant_classification_basis": "vcep_specifications",
                        "vua_benign_population_review": "does_not_meet",
                    },
                },
                {"code": "PP1", "enabled": True, "evidence": {"likelihood_ratio": 2.08}},
            ],
            None,
        ),
        (
            [{"name": "PP3", "applies": True, "strength": "Supporting", "points": 1}],
            PVS1_RNA,
            None,
        ),
        (
            [
                {"name": "BP4", "applies": True, "strength": "Supporting", "points": -1},
                {"name": "BP7", "applies": True, "strength": "Supporting", "points": -1},
            ],
            [{
                "code": "BP7_RNA",
                "enabled": True,
                "evidence": {
                    "assay_scope": "mrna_only",
                    "rna_conclusion": "no_damaging_effect",
                    "transcript_accession": "NM_007294.4",
                    "tissue_or_cell_type": "blood",
                    "nmd_assessed": "not_applicable",
                },
                }],
            {
                "gene": "BRCA1",
                "c_notation": "c.4185G>A",
                "p_notation": "p.(Gln1395=)",
            },
        ),
        (
            [],
            [{
                "code": "PP4",
                "enabled": True,
                "evidence": {
                    "combined_clinical_lr": 350,
                    "source_review_status": "enigma_recognised",
                    "source_pmid": "31853058",
                    "clinical_data_summary": "Variant-specific combined clinical evidence; overlap reviewed.",
                },
                "references": ["PMID:31853058"],
            }],
            None,
        ),
    ),
)
def test_manual_evidence_dag_has_exact_oracle_parity(base, manual, variant_context):
    expected = evaluate_manual_evidence(base, manual, variant_context)
    execution = execute_manual_evidence(base, manual, variant_context)
    assert compare_classification_results(expected, execution.result) == ()
    assert [entry.node_id for entry in execution.trace] == [
        "contract.manual_evidence_inputs",
        "rule.manual_evidence",
    ]


def test_manual_evidence_dag_fails_closed_on_missing_pp4_provenance():
    with pytest.raises(DagNodeExecutionError, match="recorded source"):
        execute_manual_evidence(
            [],
            [{
                "code": "PP4",
                "enabled": True,
                "evidence": {"combined_clinical_lr": 350},
            }],
        )


@pytest.mark.parametrize(
    ("code", "evidence", "override_strength"),
    (
        ("PS4", {"odds_ratio": 1.1, "p_value": 0.9, "lower_ci": 0.5}, "Strong"),
        ("PP1", {"likelihood_ratio": 1.0}, "Very Strong"),
        ("PM3", {"evidence_points": 0}, "Strong"),
        ("BS2", {"evidence_points": 0}, "Strong"),
        ("BS4", {"likelihood_ratio": 1.0}, "Very Strong"),
    ),
)
def test_manual_evidence_dag_rejects_strength_override(
    code, evidence, override_strength
):
    with pytest.raises(DagNodeExecutionError, match="overrides are not permitted"):
        execute_manual_evidence(
            [],
            [{
                "code": code,
                "enabled": True,
                "evidence": evidence,
                "override_strength": override_strength,
            }],
        )


def test_manual_evidence_api_contract_rejects_legacy_strength_override():
    with pytest.raises(ValidationError, match="overrides are not permitted"):
        ManualCriterionInput.model_validate({
            "code": "PP1",
            "enabled": True,
            "evidence": {"likelihood_ratio": 1.0},
            "override_strength": "Very Strong",
        })


@pytest.mark.parametrize(
    ("code", "evidence", "expected_strength", "expected_points"),
    (
        (
            "PS4",
            {
                "odds_ratio": 4.0,
                "p_value": 0.05,
                "lower_ci": 2.01,
                "case_control_country_matched": True,
                "case_control_ethnicity_matched": True,
            },
            "Strong",
            4,
        ),
        ("PP1", {"likelihood_ratio": 18.7}, "Strong", 4),
        (
            "PM3",
            {
                "evidence_points": 4.0,
                "cooccurring_variant_classification_basis": "vcep_specifications",
                "vua_benign_population_review": "does_not_meet",
            },
            "Strong",
            4,
        ),
        (
            "BS2",
            {
                "evidence_points": 4.0,
                "cooccurring_variant_classification_basis": "vcep_specifications",
            },
            "Strong",
            -4,
        ),
        ("BS4", {"likelihood_ratio": 0.05}, "Strong", -4),
    ),
)
def test_manual_evidence_dag_derives_strength_from_enigma_thresholds(
    code, evidence, expected_strength, expected_points
):
    execution = execute_manual_evidence(
        [],
        [{"code": code, "enabled": True, "evidence": evidence}],
    )

    criterion = execution.result["manual_criteria"][0]
    assert criterion["applies"] is True
    assert criterion["suggested_strength"] == expected_strength
    assert criterion["selected_strength"] == expected_strength
    assert criterion["points"] == expected_points
    assert criterion["overridden"] is False


def test_single_aggregate_bs4_strong_does_not_alone_reach_likely_benign():
    execution = execute_manual_evidence(
        [],
        [{
            "code": "BS4",
            "enabled": True,
            "evidence": {"likelihood_ratio": 0.05},
        }],
    )

    criterion = execution.result["manual_criteria"][0]
    assert criterion["selected_strength"] == "Strong"
    assert criterion["single_strong_likely_benign_eligible"] is False
    assert execution.result["predicted_class"] == 3


def test_multiple_independent_bs4_lrs_allow_single_strong_likely_benign():
    execution = execute_manual_evidence(
        [],
        [{
            "code": "BS4",
            "enabled": True,
            "evidence": {
                "likelihood_ratio": 0.05,
                "likelihood_ratio_components": [
                    {
                        "likelihood_ratio": 0.2,
                        "source": "Family A segregation analysis",
                        "independence_group": "family-a",
                    },
                    {
                        "likelihood_ratio": 0.25,
                        "source": "Family B segregation analysis",
                        "independence_group": "family-b",
                    },
                ],
            },
        }],
    )

    criterion = execution.result["manual_criteria"][0]
    assert criterion["selected_strength"] == "Strong"
    assert criterion["single_strong_likely_benign_eligible"] is True
    assert criterion["independent_evidence_contribution_count"] == 2
    assert execution.result["predicted_class"] == 2


def test_bs4_component_product_must_match_reported_combined_lr():
    with pytest.raises(DagNodeExecutionError) as exc_info:
        execute_manual_evidence(
            [],
            [{
                "code": "BS4",
                "enabled": True,
                "evidence": {
                    "likelihood_ratio": 0.04,
                    "likelihood_ratio_components": [
                        {
                            "likelihood_ratio": 0.2,
                            "source": "Family A",
                            "independence_group": "family-a",
                        },
                        {
                            "likelihood_ratio": 0.25,
                            "source": "Family B",
                            "independence_group": "family-b",
                        },
                    ],
                },
            }],
        )
    assert "does not equal the product" in str(exc_info.value.__cause__)


@pytest.mark.parametrize(
    ("code", "evidence"),
    (
        (
            "PS4",
            {
                "odds_ratio": 5.0,
                "p_value": 0.01,
                "lower_ci": 2.1,
                "case_control_country_matched": True,
            },
        ),
        ("PM3", {"evidence_points": 4.0}),
        (
            "PM3",
            {
                "evidence_points": 4.0,
                "cooccurring_variant_classification_basis": "vcep_specifications",
                "vua_benign_population_review": "meets",
            },
        ),
        ("BS2", {"evidence_points": 4.0}),
    ),
)
def test_manual_evidence_dag_does_not_apply_when_enigma_stipulation_is_missing(
    code, evidence
):
    execution = execute_manual_evidence(
        [],
        [{"code": code, "enabled": True, "evidence": evidence}],
    )

    criterion = execution.result["manual_criteria"][0]
    assert criterion["applies"] is False
    assert criterion["selected_strength"] is None
    assert criterion["points"] == 0


def test_pp1_very_strong_requires_predicted_or_experimental_effect():
    without_effect = execute_manual_evidence(
        [],
        [{"code": "PP1", "enabled": True, "evidence": {"likelihood_ratio": 350}}],
    )
    with_effect = execute_manual_evidence(
        [],
        [{
            "code": "PP1",
            "enabled": True,
            "evidence": {
                "likelihood_ratio": 350,
                "very_strong_effect_basis": "experimental_splicing",
            },
        }],
    )

    assert without_effect.result["manual_criteria"][0]["selected_strength"] == "Strong"
    assert with_effect.result["manual_criteria"][0]["selected_strength"] == "Very Strong"


def test_manual_evidence_dag_preserves_standalone_ba1_without_manual_criteria():
    execution = execute_manual_evidence(
        [{
            "name": "BA1",
            "applies": True,
            "strength": "Stand Alone",
            "points": -99,
            "reason": "gnomAD non-founder FAF95 exceeds the ENIGMA BA1 threshold",
        }],
        [],
    )

    assert execution.result["predicted_class"] == 1
    assert execution.result["predicted_label"] == "Benign"
    assert execution.result["total_points"] == -99
    assert "stand-alone" in execution.result["classification_note"]
    assert "Mixed" not in execution.result["classification_note"]


def test_manual_evidence_dag_preserves_standalone_ba1_with_pathogenic_evidence():
    execution = execute_manual_evidence(
        [{
            "name": "BA1",
            "applies": True,
            "strength": "Stand Alone",
            "points": -99,
            "reason": "gnomAD non-founder FAF95 exceeds the ENIGMA BA1 threshold",
        }],
        [{
            "code": "PP1",
            "enabled": True,
            "evidence": {"likelihood_ratio": 2.08},
            "notes": "Synthetic regression evidence for the BA1 precedence test.",
            "references": ["ENIGMA BRCA1/2 VCEP v1.2"],
        }],
    )

    assert execution.result["predicted_class"] == 1
    assert execution.result["predicted_label"] == "Benign"
    assert execution.result["total_points"] == -98
    assert execution.result["manual_criteria"][0]["code"] == "PP1"
    assert execution.result["manual_criteria"][0]["applies"] is True
    assert "stand-alone" in execution.result["classification_note"]
    assert "Mixed" not in execution.result["classification_note"]
