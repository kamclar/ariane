"""Shared immutable values and helpers for native classification nodes."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Mapping

from backend.classification_dag.domain import (
    CriterionDecision,
    CriterionDecisionStatus,
    EvidenceBundle,
)
from backend.gene_policy import policy_version


FIRST_PASS_WARNING = (
    "FIRST PASS - automatable ENIGMA VCEP v1.2 rules only. "
    "The following criteria are NOT automated and require expert review: "
    "PS4 (case-control data), PM3 (Fanconi anemia / trans variants), "
    "PP1 (co-segregation), BS2 (healthy carriers), BS4 (segregation absence). "
    "This automated result must not replace a full expert variant classification."
)
FIRST_PASS_NOTE = (
    "First pass - automatable ENIGMA VCEP v1.2 rules only. "
    "Non-automated criteria (PS4, PM3, PP1, BS2, BS4) may affect final classification."
)


@dataclass(frozen=True)
class SpliceContext:
    service_score: float | None
    effective_score: float | None
    table9_score: float | None
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class RetainedEvidence:
    criteria: tuple[CriterionDecision, ...]
    excluded_criteria: tuple[CriterionDecision, ...]
    not_applicable_criteria: tuple[CriterionDecision, ...]
    warnings: tuple[str, ...]
    evidence_interactions: tuple[Mapping[str, Any], ...]
    has_functional_evidence: bool
    metadata: Mapping[str, Any] = field(default_factory=dict)


def decision(
    code: str,
    value: Mapping[str, Any],
    *,
    gene: str,
    family_id: str,
    status: CriterionDecisionStatus = CriterionDecisionStatus.APPLIED,
    evidence_item_ids: tuple[str, ...] = (),
) -> CriterionDecision:
    payload = dict(value)
    if status == CriterionDecisionStatus.EXCLUDED:
        payload["points"] = 0
    base = CriterionDecision.from_public_mapping(code, payload, status=status)
    return CriterionDecision(
        code=base.code,
        status=base.status,
        strength=base.strength,
        points=base.points,
        reason=base.reason,
        source=base.source,
        rule_id=family_id,
        rule_version=policy_version(gene),
        evidence_item_ids=evidence_item_ids,
        decision_path=base.decision_path,
        raw_payload=base.raw_payload,
    )


def criteria_dict(
    decisions: tuple[CriterionDecision, ...],
) -> dict[str, dict[str, Any]]:
    return {item.code: item.as_public_dict() for item in decisions}


def bundle_value(bundle: EvidenceBundle, evidence_id: str) -> Any:
    item = bundle.get(evidence_id)
    if item is None:
        raise ValueError(f"EvidenceBundle is missing required item {evidence_id!r}")
    return item.value


def excluded_by_policy(decision_item: CriterionDecision, gene: str) -> CriterionDecision:
    reason = (
        f"{decision_item.code} is not applicable under the configured "
        f"VCEP policy for {gene}"
    )
    payload = dict(decision_item.raw_payload)
    payload.update({"applies": False, "points": 0, "reason": reason})
    return replace(
        decision_item,
        status=CriterionDecisionStatus.EXCLUDED,
        points=0,
        reason=reason,
        raw_payload=payload,
    )
