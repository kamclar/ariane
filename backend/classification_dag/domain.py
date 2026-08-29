"""Stable internal records for evidence, criterion decisions, and assertions."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Optional


class EvidenceStatus(str, Enum):
    """State of one source-derived fact before criterion evaluation."""

    AVAILABLE = "available"
    NOT_PROVIDED = "not_provided"
    NOT_APPLICABLE = "not_applicable"
    UNAVAILABLE = "unavailable"
    AMBIGUOUS = "ambiguous"
    REVIEW_REQUIRED = "review_required"


@dataclass(frozen=True)
class NormalizedVariant:
    """Canonical variant identity passed from normalization to classification."""

    gene: str
    reference_transcript: str
    c_notation: str
    p_notation: str
    variant_type: str
    submitted_notation: str = ""
    normalization_source: str = ""
    consequence_status: str = ""
    normalization_provenance: Mapping[str, str] = field(default_factory=dict)
    protein_consequence_explanation: str = ""
    assembly: str = ""
    genomic_notation: str = ""

    @property
    def variant_key(self) -> str:
        return f"{self.gene}:{self.c_notation}"


@dataclass(frozen=True)
class ClassificationInputs:
    """Validated facts consumed by the synchronous classification rule graph."""

    gene: str
    variant_type: str
    p_notation: str
    c_notation: str
    spliceai_score: Optional[float] = None
    bayesdel_score: Optional[float] = None
    gnomad_data: Optional[Mapping[str, Any]] = None
    frequency_policy: Optional[Mapping[str, Any]] = None
    table9_result: Optional[Mapping[str, Any]] = None
    pp4_bp5_result: Optional[Mapping[str, Any]] = None
    ps1_result: Optional[Mapping[str, Any]] = None
    exon_cnv_result: Optional[Mapping[str, Any]] = None
    residue_info: Optional[Mapping[str, Any]] = None
    dup_type: str = "Unknown"
    reference_transcript: str = ""
    submitted_notation: str = ""
    normalization_source: str = ""
    consequence_status: str = ""
    normalization_provenance: Optional[Mapping[str, str]] = None
    protein_consequence_explanation: str = ""
    assembly: str = ""
    genomic_notation: str = ""

    @property
    def variant_key(self) -> str:
        return f"{self.gene}:{self.c_notation}"

    def normalized_variant(self) -> NormalizedVariant:
        return NormalizedVariant(
            gene=self.gene,
            reference_transcript=self.reference_transcript,
            c_notation=self.c_notation,
            p_notation=self.p_notation,
            variant_type=self.variant_type,
            submitted_notation=self.submitted_notation,
            normalization_source=self.normalization_source,
            consequence_status=self.consequence_status,
            normalization_provenance=dict(self.normalization_provenance or {}),
            protein_consequence_explanation=self.protein_consequence_explanation,
            assembly=self.assembly,
            genomic_notation=self.genomic_notation,
        )

    def evidence_bundle(self) -> "EvidenceBundle":
        raw_items = (
            ("spliceai", "splice_prediction", self.spliceai_score),
            ("bayesdel", "protein_prediction", self.bayesdel_score),
            ("gnomad", "population_frequency", self.gnomad_data),
            ("clinical_lr", "clinical_likelihood_ratio", self.pp4_bp5_result),
            ("protein_ps1", "protein_reference", self.ps1_result),
            ("exon_cnv", "copy_number", self.exon_cnv_result),
            ("residue", "protein_residue", self.residue_info),
        )
        return EvidenceBundle(tuple(
            EvidenceItem(
                id=evidence_id,
                kind=kind,
                status=(
                    EvidenceStatus.AVAILABLE
                    if value is not None
                    else EvidenceStatus.NOT_PROVIDED
                ),
                value=value,
                reason=(
                    "Evidence value was supplied to the classification DAG."
                    if value is not None
                    else "No evidence value was supplied to the classification DAG."
                ),
                provenance={"adapter": "classification-inputs-v1"},
            )
            for evidence_id, kind, value in raw_items
        ))


class CriterionDecisionStatus(str, Enum):
    APPLIED = "applied"
    EXCLUDED = "excluded"
    NOT_MET = "not_met"
    NOT_APPLICABLE = "not_applicable"
    UNAVAILABLE = "unavailable"
    AMBIGUOUS = "ambiguous"
    REVIEW_REQUIRED = "review_required"
    ERROR = "error"


@dataclass(frozen=True)
class EvidenceItem:
    """One source-derived fact, before an ACMG/AMP criterion is decided."""

    id: str
    kind: str
    status: EvidenceStatus
    value: Any = None
    source_id: str = ""
    source_version: str = ""
    source_checksum: str = ""
    reason: str = ""
    provenance: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class EvidenceBundle:
    """Immutable collection of uniquely identified evidence inputs."""

    items: tuple[EvidenceItem, ...] = ()

    def __post_init__(self) -> None:
        identifiers = [item.id for item in self.items]
        duplicates = sorted({item_id for item_id in identifiers if identifiers.count(item_id) > 1})
        if duplicates:
            raise ValueError(
                "EvidenceBundle contains duplicate evidence ids: "
                + ", ".join(duplicates)
            )

    def get(self, evidence_id: str) -> EvidenceItem | None:
        return next((item for item in self.items if item.id == evidence_id), None)


@dataclass(frozen=True)
class CriterionDecision:
    """One explicit decision made by a versioned rule from evidence items."""

    code: str
    status: CriterionDecisionStatus
    strength: str = ""
    points: int = 0
    reason: str = ""
    source: str = ""
    rule_id: str = ""
    rule_version: str = ""
    evidence_item_ids: tuple[str, ...] = ()
    decision_path: Mapping[str, Any] = field(default_factory=dict)
    raw_payload: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_public_mapping(
        cls,
        code: str,
        value: Mapping[str, Any],
        *,
        status: CriterionDecisionStatus,
    ) -> "CriterionDecision":
        points = value.get("points", 0)
        if not isinstance(points, int):
            raise TypeError(f"Criterion {code} points must be an integer")
        if status == CriterionDecisionStatus.EXCLUDED and points != 0:
            raise ValueError(f"Excluded criterion {code} must have zero points")
        return cls(
            code=code,
            status=status,
            strength=str(value.get("strength") or ""),
            points=points,
            reason=str(value.get("reason") or ""),
            source=str(value.get("source") or ""),
            rule_id="public-result-contract",
            rule_version="1",
            decision_path=(
                dict(value["decision_path"])
                if isinstance(value.get("decision_path"), Mapping)
                else {}
            ),
            raw_payload=dict(value),
        )

    def as_public_dict(self) -> dict[str, Any]:
        """Return the established API representation without losing rule detail."""
        if self.raw_payload:
            return dict(self.raw_payload)
        value: dict[str, Any] = {
            "applies": self.status == CriterionDecisionStatus.APPLIED,
            "strength": self.strength,
            "points": self.points,
            "reason": self.reason,
        }
        if self.source:
            value["source"] = self.source
        if self.decision_path:
            value["decision_path"] = dict(self.decision_path)
        return value


@dataclass(frozen=True)
class CriterionFamilyResult:
    """Immutable output of one clinically meaningful DAG rule family."""

    family_id: str
    criteria: tuple[CriterionDecision, ...] = ()
    excluded_criteria: tuple[CriterionDecision, ...] = ()
    not_applicable_criteria: tuple[CriterionDecision, ...] = ()
    warnings: tuple[str, ...] = ()
    evidence_interactions: tuple[Mapping[str, Any], ...] = ()
    has_functional_evidence: bool = False
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        applied_codes = [decision.code for decision in self.criteria]
        excluded_codes = [decision.code for decision in self.excluded_criteria]
        not_applicable_codes = [
            decision.code for decision in self.not_applicable_criteria
        ]
        duplicate_applied = sorted({
            code for code in applied_codes if applied_codes.count(code) > 1
        })
        duplicate_excluded = sorted({
            code for code in excluded_codes if excluded_codes.count(code) > 1
        })
        duplicate_not_applicable = sorted({
            code
            for code in not_applicable_codes
            if not_applicable_codes.count(code) > 1
        })
        overlap = sorted(
            (set(applied_codes) & set(excluded_codes))
            | (set(applied_codes) & set(not_applicable_codes))
            | (set(excluded_codes) & set(not_applicable_codes))
        )
        invalid_statuses = [
            f"{decision.code}:{decision.status.value}"
            for decisions, expected in (
                (self.criteria, CriterionDecisionStatus.APPLIED),
                (self.excluded_criteria, CriterionDecisionStatus.EXCLUDED),
                (
                    self.not_applicable_criteria,
                    CriterionDecisionStatus.NOT_APPLICABLE,
                ),
            )
            for decision in decisions
            if decision.status != expected
        ]
        if (
            duplicate_applied
            or duplicate_excluded
            or duplicate_not_applicable
            or overlap
            or invalid_statuses
        ):
            raise ValueError(
                f"Invalid criterion family {self.family_id}: "
                f"duplicate_applied={duplicate_applied}, "
                f"duplicate_excluded={duplicate_excluded}, "
                f"duplicate_not_applicable={duplicate_not_applicable}, "
                f"overlap={overlap}, invalid_statuses={invalid_statuses}"
            )

    @property
    def has_ba1(self) -> bool:
        return any(decision.code == "BA1" for decision in self.criteria)


@dataclass(frozen=True)
class VariantAssertion:
    """Typed classification assertion produced after criterion combination."""

    variant: str
    gene: str
    c_notation: str
    p_notation: str
    criteria: tuple[CriterionDecision, ...]
    excluded_criteria: tuple[CriterionDecision, ...]
    not_applicable_criteria: tuple[CriterionDecision, ...]
    predicted_class: int
    predicted_label: str
    total_points: int
    evidence_direction: str
    mixed_evidence: bool
    classification_note: str
    policy_id: str
    policy_version: str

    @classmethod
    def from_public_result(
        cls,
        result: Mapping[str, Any],
        *,
        policy_id: str,
        policy_version: str,
    ) -> "VariantAssertion":
        required_fields = {
            "variant",
            "gene",
            "c_notation",
            "p_notation",
            "criteria",
            "excluded_criteria",
            "total_points",
            "warnings",
            "predicted_class",
            "predicted_label",
        }
        missing = sorted(required_fields.difference(result))
        if missing:
            raise ValueError(
                "Classification result is missing fields: " + ", ".join(missing)
            )
        if not isinstance(result["criteria"], Mapping):
            raise TypeError("criteria must be a mapping")
        if not isinstance(result["excluded_criteria"], Mapping):
            raise TypeError("excluded_criteria must be a mapping")
        if not isinstance(result.get("not_applicable_criteria", {}), Mapping):
            raise TypeError("not_applicable_criteria must be a mapping")
        if not isinstance(result["warnings"], list):
            raise TypeError("warnings must be a list")
        if not isinstance(result["total_points"], int):
            raise TypeError("total_points must be an integer")
        if result["predicted_class"] not in {1, 2, 3, 4, 5}:
            raise ValueError("predicted_class must be between 1 and 5")

        criteria = tuple(
            CriterionDecision.from_public_mapping(
                code,
                value,
                status=CriterionDecisionStatus.APPLIED,
            )
            for code, value in result["criteria"].items()
        )
        excluded = tuple(
            CriterionDecision.from_public_mapping(
                code,
                {**value, "points": 0},
                status=CriterionDecisionStatus.EXCLUDED,
            )
            for code, value in result["excluded_criteria"].items()
        )
        not_applicable = tuple(
            CriterionDecision.from_public_mapping(
                code,
                {**value, "points": 0},
                status=CriterionDecisionStatus.NOT_APPLICABLE,
            )
            for code, value in result.get("not_applicable_criteria", {}).items()
        )
        calculated_points = sum(criterion.points for criterion in criteria)
        if calculated_points != result["total_points"]:
            raise ValueError(
                "total_points does not equal the sum of applied criteria: "
                f"{result['total_points']} != {calculated_points}"
            )

        return cls(
            variant=str(result["variant"]),
            gene=str(result["gene"]),
            c_notation=str(result["c_notation"]),
            p_notation=str(result["p_notation"]),
            criteria=criteria,
            excluded_criteria=excluded,
            not_applicable_criteria=not_applicable,
            predicted_class=int(result["predicted_class"]),
            predicted_label=str(result["predicted_label"]),
            total_points=result["total_points"],
            evidence_direction=str(result.get("evidence_direction") or "none"),
            mixed_evidence=bool(result.get("mixed_evidence", False)),
            classification_note=str(result.get("classification_note") or ""),
            policy_id=policy_id,
            policy_version=policy_version,
        )
