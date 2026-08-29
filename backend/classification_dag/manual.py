"""Fail-closed DAG for review-adjusted classification with manual evidence."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence
import uuid

from backend.classification_dag.engine import DagDefinition, DagExecutor
from backend.classification_dag.policy import classify_by_enigma_combination
from backend.classification_dag.types import DagExecutionContext, DagTraceEntry, NodeResult
from backend.gene_policy import resolve_policy_identity
from backend.modules.bp7_rna import evaluate_bp7_rna_variant_context
from backend.modules.criterion_order import criterion_sort_key
from backend.modules.evidence_interactions import apply_manual_rna_interactions
from backend.modules.manual_evidence import (
    MANUAL_CRITERIA,
    STRENGTH_POINTS,
    STRUCTURED_CURATED_CODES,
    _pp4_source_is_recorded,
    _pp4_value_and_scale,
    evaluate_bs4_likelihood_ratio,
    suggest_strength,
)


@dataclass(frozen=True)
class ManualEvidenceInputs:
    base_criteria: tuple[Mapping[str, Any], ...]
    manual_criteria: tuple[Mapping[str, Any], ...]
    variant_context: Mapping[str, Any] | None

    @classmethod
    def create(
        cls,
        base_criteria: Sequence[Mapping[str, Any]],
        manual_criteria: Sequence[Mapping[str, Any]],
        variant_context: Mapping[str, Any] | None = None,
    ) -> "ManualEvidenceInputs":
        return cls(
            tuple(dict(item) for item in base_criteria),
            tuple(dict(item) for item in manual_criteria),
            dict(variant_context) if variant_context is not None else None,
        )


@dataclass(frozen=True)
class EvaluatedManualEvidence:
    base_combined: Mapping[str, Mapping[str, Any]]
    decisions: tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class InteractedManualEvidence:
    combined: Mapping[str, Mapping[str, Any]]
    decisions: tuple[Mapping[str, Any], ...]
    interactions: tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class ManualEvidenceExecution:
    result: Mapping[str, Any]
    graph_id: str
    graph_version: str
    trace: tuple[DagTraceEntry, ...]


@dataclass(frozen=True)
class ManualInputContractNode:
    id: str = "contract.manual_evidence_inputs"
    version: str = "1"
    requires: frozenset[str] = frozenset({"manual_inputs"})
    provides: frozenset[str] = frozenset({"validated_manual_inputs"})

    def evaluate(self, context, inputs) -> NodeResult:
        value = inputs["manual_inputs"]
        if not isinstance(value, ManualEvidenceInputs):
            raise TypeError("manual_inputs must be ManualEvidenceInputs")
        for item in value.base_criteria:
            if not isinstance(item.get("name"), str):
                raise ValueError("Every base criterion requires a name")
        for item in value.manual_criteria:
            if not isinstance(item.get("code"), str):
                raise ValueError("Every manual criterion requires a code")
        return NodeResult.succeeded({"validated_manual_inputs": value})


@dataclass(frozen=True)
class ManualCriterionEvaluationNode:
    id: str = "rule.manual_evidence"
    version: str = "1"
    requires: frozenset[str] = frozenset({"validated_manual_inputs"})
    provides: frozenset[str] = frozenset({"evaluated_manual_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        source = inputs["validated_manual_inputs"]
        combined = {
            item["name"]: {
                "applies": item.get("applies", True),
                "strength": item.get("strength"),
                "points": item.get("points", 0),
                "reason": item.get("reason", ""),
                "single_strong_likely_benign_eligible": item.get(
                    "single_strong_likely_benign_eligible", False
                ),
                "single_strong_likely_benign_basis": item.get(
                    "single_strong_likely_benign_basis", ""
                ),
                "independent_evidence_contribution_count": item.get(
                    "independent_evidence_contribution_count", 0
                ),
            }
            for item in source.base_criteria
            if item.get("applies", True)
        }
        enabled = [item for item in source.manual_criteria if item.get("enabled")]
        clinical_lr_used = any(code in combined for code in {"PP4", "BP5"}) or any(
            item.get("code") == "PP4" for item in enabled
        )
        if clinical_lr_used:
            for item in enabled:
                if item.get("code") not in {"PP1", "PS4"}:
                    continue
                evidence = item.get("evidence", {})
                rationale = str(evidence.get("independence_rationale") or "").strip()
                if evidence.get("independent_from_pp4_bp5") is not True or not rationale:
                    raise ValueError(
                        f"{item['code']} cannot be combined with PP4/BP5 until the reviewer "
                        "confirms independent observations and records an independence rationale"
                    )

        decisions: list[Mapping[str, Any]] = []
        for item in source.manual_criteria:
            code = item["code"]
            definition = MANUAL_CRITERIA[code]
            if item.get("override_strength") not in {None, ""}:
                raise ValueError(
                    "Manual strength overrides are not permitted. "
                    "ARIANE derives criterion strength from the configured "
                    "VCEP evidence thresholds."
                )
            suggested = suggest_strength(
                code,
                item.get("evidence", {}),
                variant_context=source.variant_context,
                base_criteria=source.base_criteria,
            )
            evidence = item.get("evidence", {})
            pp4_value, pp4_scale = _pp4_value_and_scale(evidence)
            pp4_status = str(evidence.get("source_review_status") or "unreviewed").strip().lower()
            pp4_complete = (
                pp4_value is not None
                and pp4_value >= 0
                and pp4_scale in {"lr", "log10_lr", "acmg_points"}
                and pp4_status in {"appendix_b", "other_reviewed", "unreviewed"}
                and _pp4_source_is_recorded(evidence)
                and bool((evidence.get("clinical_data_summary") or "").strip())
            )
            if code == "PP4" and item.get("enabled") and not pp4_complete:
                raise ValueError(
                    "PP4 requires a clinical LR value and scale, recorded source, and clinical data summary"
                )
            if code == "BP7_RNA" and item.get("enabled") and not suggested:
                context_result = evaluate_bp7_rna_variant_context(
                    source.variant_context, source.base_criteria
                )
                if not context_result["eligible"]:
                    raise ValueError(context_result["reason"])
            if code in STRUCTURED_CURATED_CODES - {"PP4"} and item.get("enabled") and not suggested:
                raise ValueError(f"{code} requires a complete structured curated evidence record")
            selected = suggested
            applies = bool(item.get("enabled") and selected)
            points = STRENGTH_POINTS.get(selected, 0)
            if definition["direction"] == "benign":
                points *= -1
            reason = f"User-provided evidence; ARIANE suggestion: {suggested or 'threshold not met'}"
            single_strong_eligible = False
            single_strong_basis = ""
            contribution_count = 0
            if code == "BS4" and applies:
                _, has_multiple_lrs, contribution_count = evaluate_bs4_likelihood_ratio(
                    evidence
                )
                single_strong_eligible = bool(
                    selected == "Strong" and has_multiple_lrs
                )
                if single_strong_eligible:
                    single_strong_basis = (
                        "Multiple independently identified segregation likelihood ratios "
                        "contribute to BS4 Strong"
                    )
                    reason += "; multiple independent LR components satisfy the ENIGMA Table 3 single-Strong condition"
                elif selected == "Strong":
                    reason += "; BS4 Strong is valid, but the ENIGMA Table 3 single-Strong condition is not documented"
            decisions.append({
                "code": code,
                "applies": applies,
                "suggested_strength": suggested,
                "selected_strength": selected,
                "points": points if applies else 0,
                "reason": reason,
                "threshold_note": definition["threshold"],
                "overridden": False,
                "notes": item.get("notes", ""),
                "references": item.get("references", []),
                "single_strong_likely_benign_eligible": single_strong_eligible,
                "single_strong_likely_benign_basis": single_strong_basis,
                "independent_evidence_contribution_count": contribution_count,
            })
            if applies:
                if code == "PVS1_INIT":
                    combined.pop("PP3", None)
                output_code = "PS1" if code == "PS1_PROTEIN" else code
                if output_code == "PS1" and output_code in combined:
                    raise ValueError(
                        "Protein PS1 is already present in the automated result and cannot be counted twice"
                    )
                combined[output_code] = {
                    "applies": True,
                    "strength": selected,
                    "points": points,
                    "reason": reason,
                    "single_strong_likely_benign_eligible": single_strong_eligible,
                    "single_strong_likely_benign_basis": single_strong_basis,
                    "independent_evidence_contribution_count": contribution_count,
                }
        value = EvaluatedManualEvidence(combined, tuple(decisions))
        return NodeResult.succeeded(
            {"evaluated_manual_evidence": value},
            provenance={"applied": [item["code"] for item in decisions if item["applies"]]},
        )


@dataclass(frozen=True)
class ManualEvidenceInteractionNode:
    id: str = "policy.manual_evidence_interactions"
    version: str = "1"
    requires: frozenset[str] = frozenset({"evaluated_manual_evidence"})
    provides: frozenset[str] = frozenset({"interacted_manual_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        source = inputs["evaluated_manual_evidence"]
        combined = {code: dict(item) for code, item in source.base_combined.items()}
        applied = {item["code"] for item in source.decisions if item["applies"]}
        interactions = apply_manual_rna_interactions(combined, applied)
        value = InteractedManualEvidence(combined, source.decisions, tuple(interactions))
        return NodeResult.succeeded(
            {"interacted_manual_evidence": value},
            provenance={"interaction_count": len(interactions)},
        )


@dataclass(frozen=True)
class ManualClassificationPolicyNode:
    id: str = "policy.manual_enigma_combination"
    version: str = "1"
    requires: frozenset[str] = frozenset({"interacted_manual_evidence"})
    provides: frozenset[str] = frozenset({"manual_classification_result"})

    def evaluate(self, context, inputs) -> NodeResult:
        source = inputs["interacted_manual_evidence"]
        total = sum(item.get("points", 0) for item in source.combined.values())
        cls, label, note = classify_by_enigma_combination(
            source.combined, total, gene=context.metadata.get("gene")
        )
        result = {
            "predicted_class": cls,
            "predicted_label": label,
            "total_points": total,
            "classification_note": note,
            "manual_criteria": sorted(source.decisions, key=lambda item: criterion_sort_key(item["code"])),
            "evidence_interactions": [dict(item) for item in source.interactions],
        }
        return NodeResult.succeeded(
            {"manual_classification_result": result},
            provenance={"method": "ENIGMA-Table-3-or-mixed-points"},
        )


def build_manual_evidence_graph() -> DagDefinition:
    return DagDefinition(
        id="ariane.vcep.manual-evidence",
        version="2.0.0-gene-policy-dag",
        seed_keys={"manual_inputs"},
        nodes=(
            ManualInputContractNode(),
            ManualCriterionEvaluationNode(),
            ManualEvidenceInteractionNode(),
            ManualClassificationPolicyNode(),
        ),
    )


def execute_manual_evidence(
    base_criteria: Sequence[Mapping[str, Any]],
    manual_criteria: Sequence[Mapping[str, Any]],
    variant_context: Mapping[str, Any] | None = None,
) -> ManualEvidenceExecution:
    gene = str((variant_context or {}).get("gene") or "").strip().upper() or None
    selected_policy_id, selected_policy_version = resolve_policy_identity(gene)
    graph = build_manual_evidence_graph()
    run = DagExecutor(graph).run(
        {
            "manual_inputs": ManualEvidenceInputs.create(
                base_criteria, manual_criteria, variant_context
            )
        },
        context=DagExecutionContext(
            run_id=uuid.uuid4().hex,
            variant_key="manual-review",
            policy_id=selected_policy_id,
            metadata={
                "workflow": "review-adjusted-classification",
                "policy_version": selected_policy_version,
                "gene": gene,
            },
        ),
    )
    result = run.values.get("manual_classification_result")
    if not isinstance(result, dict):
        raise RuntimeError("Manual evidence DAG did not produce a result")
    return ManualEvidenceExecution(result, run.graph_id, run.graph_version, run.trace)
