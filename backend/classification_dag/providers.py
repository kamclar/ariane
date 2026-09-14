"""Versioned evidence providers used by the classification DAG.

Provider nodes acquire and describe facts. They do not assign ACMG/AMP codes.
Every external or local source result is converted into an ``EvidenceItem`` so
missing data remains distinguishable from a rule that was evaluated and did
not meet its threshold.
"""

from __future__ import annotations

import asyncio
from dataclasses import dataclass
import logging
from typing import Any, Callable, Mapping

from backend.domain.classification import (
    ClassificationInputs,
    EvidenceBundle,
    EvidenceItem,
    EvidenceStatus,
    NormalizedVariant,
)
from backend.classification_dag.types import NodeResult
from backend.infrastructure.lookup_execution import lookup_or_unavailable
from backend.policy.spliceai import (
    spliceai_failure_is_retryable,
    spliceai_required_for_classification,
)
from backend.domain.indel import is_indel_allele
from backend.reference_data.table9 import (
    TABLE9_JSON_PATH,
    load_table9_data,
    table9_checksum,
    table9_lookup_ps3_bs3,
)


@dataclass(frozen=True)
class ClassificationRequest:
    """Normalized request entering the evidence and classification DAG."""

    variant: NormalizedVariant
    dup_type: str = "Unknown"


@dataclass(frozen=True)
class CoordinateContext:
    grch37: Mapping[str, Any] | None = None
    grch38: Mapping[str, Any] | None = None
    source: str = ""
    status: str = "not_applicable"
    diagnostics: tuple[str, ...] = ()


@dataclass(frozen=True)
class ProviderDependencies:
    """Injectable source adapters used by production and deterministic tests."""

    resolve_variant: Callable[..., Any]
    get_grch37: Callable[..., Any]
    get_grch38: Callable[..., Any]
    spliceai_lookup: Callable[..., Any]
    spliceai_status: Callable[[str, str], Mapping[str, Any]]
    bayesdel_lookup: Callable[..., Any]
    bayesdel_status: Callable[[str, str], Mapping[str, Any]]
    gnomad_lookup: Callable[..., Any]
    clinical_lr_lookup: Callable[..., Any]
    exon_cnv_lookup: Callable[..., Any]
    residue_lookup: Callable[..., Any]
    ps1_candidate_lookup: Callable[..., Any]
    splice_source_lookup: Callable[..., Any]
    select_ps1_spliceai: Callable[..., Any]
    ps1_lookup: Callable[..., Any]
    erepo_pvs1_rna_lookup: Callable[..., Any]
    st2_pvs1_rna_lookup: Callable[..., Any]

def _request(inputs) -> ClassificationRequest:
    request = inputs["classification_request"]
    if not isinstance(request, ClassificationRequest):
        raise TypeError("classification_request must be a ClassificationRequest")
    return request


def _evidence_status(value: Any, *, applicable: bool = True) -> EvidenceStatus:
    if not applicable:
        return EvidenceStatus.NOT_APPLICABLE
    return EvidenceStatus.AVAILABLE if value is not None else EvidenceStatus.UNAVAILABLE


@dataclass(frozen=True)
class ClassificationRequestContractNode:
    id: str = "contract.classification_request"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_request"})
    provides: frozenset[str] = frozenset({"normalized_variant"})

    def evaluate(self, context, inputs) -> NodeResult:
        request = _request(inputs)
        if request.variant.variant_key != context.variant_key:
            raise ValueError(
                "Normalized variant identity does not match the DAG execution context"
            )
        return NodeResult.succeeded(
            {"normalized_variant": request.variant},
            provenance={"contract": "ariane.classification-request.v1"},
        )


@dataclass(frozen=True)
class CoordinateEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.coordinates"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {"classification_request", "normalized_variant"}
    )
    provides: frozenset[str] = frozenset({"coordinate_context"})

    async def evaluate(self, context, inputs) -> NodeResult:
        request = _request(inputs)
        variant = request.variant
        is_exon_cnv = variant.variant_type.lower() in {
            "exon_deletion", "exon_duplication",
        }
        if is_exon_cnv:
            value = CoordinateContext()
            return NodeResult.succeeded(
                {"coordinate_context": value},
                provenance={"status": value.status, "reason": "uncertain breakpoints"},
            )

        diagnostics: list[str] = []
        try:
            resolved_variant = await asyncio.to_thread(
                self.dependencies.resolve_variant,
                variant.gene,
                variant.c_notation,
            )
        except Exception as exc:
            message = f"Coordinate lookup failed: {type(exc).__name__}: {exc}"
            logging.getLogger(__name__).exception(message)
            diagnostics.append(message)
            resolved_variant = None
        resolved: dict[str, Any] = {}
        source = ""
        status = "unavailable"
        if resolved_variant:
            resolved[variant.variant_key] = resolved_variant
            source = str(getattr(resolved_variant, "source", "") or "")
            status = str(getattr(resolved_variant, "status", "") or "unavailable")
            if status != "ok":
                diagnostics.extend(
                    f"Coordinate resolver: {warning}"
                    for warning in getattr(resolved_variant, "warnings", ())
                )
        grch37 = self.dependencies.get_grch37(
            resolved, variant.gene, variant.c_notation
        )
        grch38 = self.dependencies.get_grch38(
            resolved, variant.gene, variant.c_notation
        )
        value = CoordinateContext(
            grch37=grch37,
            grch38=grch38,
            source=source,
            status=status if (grch37 or grch38) else "unavailable",
            diagnostics=tuple(diagnostics),
        )
        return NodeResult.succeeded(
            {"coordinate_context": value},
            provenance={
                "status": value.status,
                "source": source,
                "grch37_available": grch37 is not None,
                "grch38_available": grch38 is not None,
            },
            warnings=tuple(diagnostics),
        )


@dataclass(frozen=True)
class Ps1CandidateEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.ps1.candidates"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_request"})
    provides: frozenset[str] = frozenset({"ps1_reference_candidates"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        candidates = tuple(self.dependencies.ps1_candidate_lookup(
            variant.gene,
            variant.c_notation,
            variant.p_notation,
            variant.variant_type,
        ))
        return NodeResult.succeeded(
            {"ps1_reference_candidates": candidates},
            provenance={"candidate_count": len(candidates)},
        )


@dataclass(frozen=True)
class ErepoPvs1RnaEvidenceNode:
    """Provide exact, locally approved ERepo RNA assertions to PVS1 rules."""

    dependencies: ProviderDependencies
    id: str = "provider.enigma.erepo_pvs1_rna"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_request"})
    provides: frozenset[str] = frozenset({"erepo_pvs1_rna_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        value = self.dependencies.erepo_pvs1_rna_lookup(
            variant.gene,
            variant.c_notation,
        )
        status = str(value.get("status") or "unavailable")
        evidence = EvidenceItem(
            id="erepo_pvs1_rna",
            kind="curated_rna_assertion",
            status=(
                EvidenceStatus.AVAILABLE
                if status == "eligible"
                else EvidenceStatus.NOT_APPLICABLE
                if status == "not_in_registry"
                else EvidenceStatus.AMBIGUOUS
                if status == "ambiguous"
                else EvidenceStatus.UNAVAILABLE
            ),
            value=value,
            source_id=str(value.get("registry_id") or "enigma-erepo-pvs1-rna"),
            source_version=str(value.get("registry_version") or ""),
            source_checksum=str(value.get("registry_sha256") or ""),
            reason=str(value.get("reason") or ""),
            provenance={
                "status": status,
                "coverage_status": value.get("coverage_status"),
                "source_url": value.get("source_url"),
                "assertion_uuid": (value.get("record") or {}).get("assertion_uuid"),
            },
        )
        return NodeResult.succeeded(
            {"erepo_pvs1_rna_evidence": evidence},
            provenance={
                "evidence_status": evidence.status.value,
                "source_id": evidence.source_id,
                "source_version": evidence.source_version,
                "source_checksum": evidence.source_checksum,
            },
        )


@dataclass(frozen=True)
class St2Pvs1RnaEvidenceNode:
    """Provide exact ENIGMA ST2/ST3 RNA facts and Table 4 context."""

    dependencies: ProviderDependencies
    id: str = "provider.enigma.st2_pvs1_rna"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_request"})
    provides: frozenset[str] = frozenset({"st2_pvs1_rna_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        value = self.dependencies.st2_pvs1_rna_lookup(
            variant.gene,
            variant.c_notation,
        )
        status = str(value.get("status") or "unavailable")
        evidence = EvidenceItem(
            id="st2_pvs1_rna",
            kind="curated_rna_dataset",
            status=(
                EvidenceStatus.AVAILABLE
                if status == "eligible"
                else EvidenceStatus.REVIEW_REQUIRED
                if status == "review_required"
                else EvidenceStatus.NOT_APPLICABLE
                if status in {"not_in_source", "not_applicable"}
                else EvidenceStatus.UNAVAILABLE
            ),
            value=value,
            source_id=str(value.get("source_id") or "enigma-st2-st3-v1.2"),
            source_version=str(value.get("source_version") or ""),
            source_checksum=str(value.get("source_checksum") or ""),
            reason=str(value.get("reason") or ""),
            provenance={
                "status": status,
                "source_url": value.get("source_url"),
                "source_row": (value.get("record") or {}).get("source_row"),
                "table4_exon": value.get("table4_exon"),
            },
        )
        return NodeResult.succeeded(
            {"st2_pvs1_rna_evidence": evidence},
            provenance={
                "evidence_status": evidence.status.value,
                "source_id": evidence.source_id,
                "source_version": evidence.source_version,
                "source_checksum": evidence.source_checksum,
            },
        )


@dataclass(frozen=True)
class SpliceAiEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.spliceai"
    version: str = "2"
    requires: frozenset[str] = frozenset(
        {
            "classification_request",
            "normalized_variant",
            "ps1_reference_candidates",
            "erepo_pvs1_rna_evidence",
            "st2_pvs1_rna_evidence",
        }
    )
    provides: frozenset[str] = frozenset(
        {"spliceai_evidence", "ps1_reference_spliceai_scores"}
    )

    async def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        candidates = tuple(inputs["ps1_reference_candidates"])
        curated_rna = inputs["erepo_pvs1_rna_evidence"]
        st2_rna = inputs["st2_pvs1_rna_evidence"]
        erepo_pvs1_rna = (
            curated_rna.status == EvidenceStatus.AVAILABLE
            and isinstance(curated_rna.value, Mapping)
            and curated_rna.value.get("status") == "eligible"
        )
        st2_pvs1_rna = (
            st2_rna.status == EvidenceStatus.AVAILABLE
            and isinstance(st2_rna.value, Mapping)
            and st2_rna.value.get("status") == "eligible"
        )
        curated_pvs1_rna = erepo_pvs1_rna or st2_pvs1_rna
        required = (
            spliceai_required_for_classification(variant.variant_type)
            and not curated_pvs1_rna
        )
        if not required:
            reason = (
                "SpliceAI is not required because exact curated PVS1 RNA "
                "evidence replaces the Figure 1A prediction branch."
                if curated_pvs1_rna
                else
                "SpliceAI is not required by the automatic ENIGMA path for "
                f"variant type {variant.variant_type}."
            )
            evidence = EvidenceItem(
                id="spliceai",
                kind="splice_prediction",
                status=EvidenceStatus.NOT_APPLICABLE,
                reason=reason,
                provenance={
                    "status": "not_applicable",
                    "required_for_classification": False,
                    "replaced_by_curated_pvs1_rna": curated_pvs1_rna,
                    "pvs1_rna_source": (
                        "erepo" if erepo_pvs1_rna else "st2"
                        if st2_pvs1_rna else None
                    ),
                    "retryable": False,
                    "reason": reason,
                },
            )
            return NodeResult.succeeded(
                {
                    "spliceai_evidence": evidence,
                    "ps1_reference_spliceai_scores": {},
                },
                provenance={
                    "evidence_status": evidence.status.value,
                    "required_for_classification": False,
                    "replaced_by_curated_pvs1_rna": curated_pvs1_rna,
                    "pvs1_rna_source": (
                        "erepo" if erepo_pvs1_rna else "st2"
                        if st2_pvs1_rna else None
                    ),
                },
            )

        async def lookup(c_notation: str):
            local_diagnostics: list[str] = []
            score = await lookup_or_unavailable(
                self.dependencies.spliceai_lookup,
                None,
                "SpliceAI",
                local_diagnostics,
                variant.gene,
                c_notation,
            )
            item_status = dict(self.dependencies.spliceai_status(
                variant.gene, c_notation
            ))
            if score is None and not item_status.get("status"):
                item_status = {
                    "status": "api_error",
                    "reason": (
                        local_diagnostics[0]
                        if local_diagnostics
                        else "SpliceAI score is unavailable"
                    ),
                }
            explicit_retryable = item_status.get("retryable")
            item_status["retryable"] = (
                explicit_retryable
                if isinstance(explicit_retryable, bool)
                else spliceai_failure_is_retryable(
                    str(item_status.get("status") or ""),
                    str(item_status.get("reason") or ""),
                )
            )
            item_status["score"] = score
            return score, item_status, tuple(local_diagnostics)

        lookups = await asyncio.gather(*(
            lookup(c_notation) for c_notation in (variant.c_notation, *candidates)
        ))
        score, status, assessed_diagnostics = lookups[0]
        reference_scores = {
            c_notation: result[0]
            for c_notation, result in zip(candidates, lookups[1:])
        }
        reference_statuses = {
            c_notation: result[1]
            for c_notation, result in zip(candidates, lookups[1:])
        }
        diagnostics = list(assessed_diagnostics)
        for result in lookups[1:]:
            diagnostics.extend(result[2])
        status["required_for_classification"] = True
        status["reference_lookup_required"] = bool(candidates)
        status["reference_lookup_complete"] = all(
            reference_scores[c_notation] is not None
            and str(reference_statuses[c_notation].get("status") or "") == "ok"
            for c_notation in candidates
        )
        status["reference_variant_statuses"] = reference_statuses
        evidence = EvidenceItem(
            id="spliceai",
            kind="splice_prediction",
            status=_evidence_status(score),
            value=score,
            source_id=str(status.get("source") or "configured-spliceai"),
            source_version=str(status.get("scoring_profile_id") or ""),
            source_checksum=str(status.get("scoring_profile_sha256") or ""),
            reason=str(status.get("reason") or (
                "SpliceAI score is available" if score is not None
                else "SpliceAI score is unavailable"
            )),
            provenance=status,
        )
        if status.get("status") not in {None, "ok"}:
            diagnostics.append(
                f"SpliceAI unavailable: status={status.get('status')}; "
                f"{status.get('reason', 'no reason reported')}"
            )
        return NodeResult.succeeded(
            {
                "spliceai_evidence": evidence,
                "ps1_reference_spliceai_scores": reference_scores,
            },
            provenance={
                "evidence_status": evidence.status.value,
                "source_id": evidence.source_id,
                "reference_lookup_count": len(candidates),
                "required_for_classification": True,
            },
            warnings=tuple(dict.fromkeys(diagnostics)),
        )


@dataclass(frozen=True)
class BayesDelEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.bayesdel"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {"classification_request", "normalized_variant"}
    )
    provides: frozenset[str] = frozenset(
        {"bayesdel_evidence", "alphamissense_annotation"}
    )

    async def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        variant_type = variant.variant_type.lower()
        protein_prediction_types = {
            "missense",
            "inframe_deletion",
            "inframe_insertion",
            "inframe_delins",
        }
        if variant_type not in protein_prediction_types:
            evidence = EvidenceItem(
                id="bayesdel",
                kind="protein_prediction",
                status=EvidenceStatus.NOT_APPLICABLE,
                reason=(
                    "BayesDel is only evaluated for missense and confirmed "
                    "in-frame variants in the ENIGMA Figure 1A protein branch."
                ),
            )
            return NodeResult.succeeded(
                {"bayesdel_evidence": evidence, "alphamissense_annotation": None},
                provenance={"evidence_status": evidence.status.value},
            )

        diagnostics: list[str] = []
        score, alphamissense = await lookup_or_unavailable(
            self.dependencies.bayesdel_lookup,
            (None, None),
            "MyVariant/BayesDel",
            diagnostics,
            variant.gene,
            variant.c_notation,
        )
        status = dict(self.dependencies.bayesdel_status(
            variant.gene, variant.c_notation
        ))
        evidence = EvidenceItem(
            id="bayesdel",
            kind="protein_prediction",
            status=_evidence_status(score),
            value=score,
            source_id="MyVariant/BayesDel_noAF",
            reason=str(status.get("reason") or (
                "BayesDel_noAF score is available" if score is not None
                else "BayesDel_noAF score is unavailable"
            )),
            provenance=status,
        )
        if status.get("status") in {"api_error", "no_grch37_coords"}:
            diagnostics.append(
                f"MyVariant/BayesDel unavailable: status={status.get('status')}; "
                f"{status.get('reason', 'no reason reported')}"
            )
        return NodeResult.succeeded(
            {
                "bayesdel_evidence": evidence,
                "alphamissense_annotation": alphamissense,
            },
            provenance={"evidence_status": evidence.status.value},
            warnings=tuple(dict.fromkeys(diagnostics)),
        )


@dataclass(frozen=True)
class GnomadEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.gnomad"
    version: str = "1"
    requires: frozenset[str] = frozenset(
        {"classification_request", "coordinate_context"}
    )
    provides: frozenset[str] = frozenset({"gnomad_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        coordinates = inputs["coordinate_context"]
        coordinates_available = (
            coordinates.grch37 is not None or coordinates.grch38 is not None
        )
        value = self.dependencies.gnomad_lookup(
            gene=variant.gene,
            c_notation=variant.c_notation,
            grch37=coordinates.grch37,
            grch38=coordinates.grch38,
        )
        evidence = EvidenceItem(
            id="gnomad",
            kind="population_frequency",
            status=(
                EvidenceStatus.AVAILABLE
                if coordinates_available
                else EvidenceStatus.NOT_APPLICABLE
            ),
            value=value,
            source_id="gnomAD non-cancer panel datasets",
            reason=(
                str(value.get("status") or "")
                if coordinates_available
                else "Genomic coordinates are unavailable; policy context remains available"
            ),
            provenance={
                "status": value.get("status") if isinstance(value, Mapping) else "not_queried",
                "policy_id": value.get("policy_id") if isinstance(value, Mapping) else None,
            },
        )
        return NodeResult.succeeded(
            {"gnomad_evidence": evidence},
            provenance={"evidence_status": evidence.status.value},
        )


@dataclass(frozen=True)
class ClinicalLrEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.clinical_lr"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_request"})
    provides: frozenset[str] = frozenset({"clinical_lr_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        value = self.dependencies.clinical_lr_lookup(
            variant.gene, variant.c_notation
        )
        evidence = EvidenceItem(
            id="clinical_lr",
            kind="clinical_likelihood_ratio",
            status=EvidenceStatus.AVAILABLE,
            value=value,
            source_id=str(value.get("source_id") or "ENIGMA clinical LR dataset"),
            source_version=str(value.get("source_version") or ""),
            source_checksum=str(value.get("source_checksum") or ""),
            reason=str(value.get("reason") or "No applicable combined clinical LR"),
        )
        return NodeResult.succeeded(
            {"clinical_lr_evidence": evidence},
            provenance={"applies": bool(value.get("applies"))},
        )


@dataclass(frozen=True)
class ExonCnvEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.exon_cnv"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_request"})
    provides: frozenset[str] = frozenset({"exon_cnv_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        applicable = is_indel_allele(variant.c_notation)
        value = (
            self.dependencies.exon_cnv_lookup(variant.gene, variant.c_notation)
            if applicable else None
        )
        evidence = EvidenceItem(
            id="exon_cnv",
            kind="structural_population",
            status=_evidence_status(value, applicable=applicable),
            value=value,
            source_id="ENIGMA Appendix G structural population evidence",
            reason=(
                str(value.get("reason") or "") if isinstance(value, Mapping)
                else "Not an insertion, deletion, duplication, or delins allele"
            ),
        )
        return NodeResult.succeeded(
            {"exon_cnv_evidence": evidence},
            provenance={"evidence_status": evidence.status.value},
        )


@dataclass(frozen=True)
class ResidueEvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.residue"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_request"})
    provides: frozenset[str] = frozenset({"residue_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        value = self.dependencies.residue_lookup(variant.gene, variant.p_notation)
        evidence = EvidenceItem(
            id="residue",
            kind="protein_residue",
            status=EvidenceStatus.AVAILABLE,
            value=value,
            source_id="curated protein residue data",
        )
        return NodeResult.succeeded({"residue_evidence": evidence})


@dataclass(frozen=True)
class ProteinPs1EvidenceNode:
    dependencies: ProviderDependencies
    id: str = "provider.protein_ps1"
    version: str = "1"
    requires: frozenset[str] = frozenset({
        "classification_request",
        "table9_evidence",
        "spliceai_evidence",
        "ps1_reference_spliceai_scores",
    })
    provides: frozenset[str] = frozenset({"protein_ps1_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = _request(inputs).variant
        table9 = inputs["table9_evidence"].value or {}
        spliceai = inputs["spliceai_evidence"]
        splice_sources = self.dependencies.splice_source_lookup(
            variant.gene, variant.c_notation, table9
        )
        score, source = self.dependencies.select_ps1_spliceai(
            spliceai.value, table9
        )
        value = self.dependencies.ps1_lookup(
            variant.gene,
            variant.c_notation,
            variant.p_notation,
            variant_type=variant.variant_type,
            spliceai_score=score,
            vua_spliceai_source=source,
            vua_splice_evidence_status=splice_sources["status"],
            vua_splice_sources_checked=splice_sources["sources_checked"],
            reference_spliceai_scores=dict(inputs["ps1_reference_spliceai_scores"]),
        )
        evidence = EvidenceItem(
            id="protein_ps1",
            kind="protein_reference",
            status=(
                EvidenceStatus.REVIEW_REQUIRED
                if value.get("review_required")
                else EvidenceStatus.AVAILABLE
            ),
            value=value,
            source_id="ENIGMA protein PS1 reference registry",
            reason=str(value.get("reason") or ""),
            provenance={
                "application_status": value.get("application_status"),
                "reference_status": value.get("reference_status"),
            },
        )
        return NodeResult.succeeded(
            {"protein_ps1_evidence": evidence},
            provenance={
                "evidence_status": evidence.status.value,
                "application_status": value.get("application_status"),
            },
        )


@dataclass(frozen=True)
class EvidenceBundleAssemblyNode:
    id: str = "contract.evidence_bundle"
    version: str = "1"
    requires: frozenset[str] = frozenset({
        "classification_request",
        "coordinate_context",
        "table9_evidence",
        "spliceai_evidence",
        "bayesdel_evidence",
        "gnomad_evidence",
        "clinical_lr_evidence",
        "erepo_pvs1_rna_evidence",
        "st2_pvs1_rna_evidence",
        "protein_ps1_evidence",
        "exon_cnv_evidence",
        "residue_evidence",
        "alphamissense_annotation",
    })
    provides: frozenset[str] = frozenset({
        "classification_inputs", "evidence_bundle", "provider_artifacts",
    })

    def evaluate(self, context, inputs) -> NodeResult:
        request = _request(inputs)
        variant = request.variant
        items = tuple(inputs[key] for key in (
            "spliceai_evidence",
            "bayesdel_evidence",
            "gnomad_evidence",
            "clinical_lr_evidence",
            "erepo_pvs1_rna_evidence",
            "st2_pvs1_rna_evidence",
            "protein_ps1_evidence",
            "exon_cnv_evidence",
            "residue_evidence",
            "table9_evidence",
        ))
        bundle = EvidenceBundle(items)

        def value(evidence_id: str) -> Any:
            return bundle.get(evidence_id).value

        gnomad_value = value("gnomad")
        classification_inputs = ClassificationInputs(
            gene=variant.gene,
            variant_type=variant.variant_type,
            p_notation=variant.p_notation,
            c_notation=variant.c_notation,
            spliceai_score=value("spliceai"),
            bayesdel_score=value("bayesdel"),
            gnomad_data=gnomad_value,
            frequency_policy=(
                (gnomad_value or {}).get("classification_policy")
                if isinstance(gnomad_value, Mapping)
                else None
            ),
            table9_result=value("enigma_table9"),
            erepo_pvs1_rna_result=value("erepo_pvs1_rna"),
            st2_pvs1_rna_result=value("st2_pvs1_rna"),
            pp4_bp5_result=value("clinical_lr"),
            ps1_result=value("protein_ps1"),
            exon_cnv_result=value("exon_cnv"),
            residue_info=value("residue"),
            dup_type=request.dup_type,
            reference_transcript=variant.reference_transcript,
            submitted_notation=variant.submitted_notation,
            normalization_source=variant.normalization_source,
            consequence_status=variant.consequence_status,
            normalization_provenance=variant.normalization_provenance,
            protein_consequence_explanation=variant.protein_consequence_explanation,
            assembly=variant.assembly,
            genomic_notation=variant.genomic_notation,
        )
        coordinates = inputs["coordinate_context"]
        provider_artifacts = {
            "grch37": coordinates.grch37,
            "grch38": coordinates.grch38,
            "coordinate_status": coordinates.status,
            "coordinate_source": coordinates.source,
            "alphamissense": inputs["alphamissense_annotation"],
            "spliceai_score": value("spliceai"),
            "bayesdel_score": value("bayesdel"),
            "spliceai_status": dict(inputs["spliceai_evidence"].provenance),
            "bayesdel_status": dict(inputs["bayesdel_evidence"].provenance),
            "gnomad_data": gnomad_value,
            "clinical_lr_result": value("clinical_lr"),
            "erepo_pvs1_rna_result": value("erepo_pvs1_rna"),
            "st2_pvs1_rna_result": value("st2_pvs1_rna"),
            "table9_result": value("enigma_table9"),
            "exon_cnv_result": value("exon_cnv"),
            "evidence_audit": {
                item.id: {
                    "status": item.status.value,
                    "source_id": item.source_id,
                    "source_version": item.source_version,
                    "source_checksum": item.source_checksum,
                    "reason": item.reason,
                }
                for item in bundle.items
            },
        }
        return NodeResult.succeeded(
            {
                "classification_inputs": classification_inputs,
                "evidence_bundle": bundle,
                "provider_artifacts": provider_artifacts,
            },
            provenance={
                "evidence_statuses": {
                    item.id: item.status.value for item in bundle.items
                }
            },
        )


@dataclass(frozen=True)
class Table9EvidenceNode:
    """Read functional evidence from the validated ENIGMA Table 9 dataset."""

    id: str = "provider.enigma.table9"
    version: str = "1"
    requires: frozenset[str] = frozenset({"normalized_variant"})
    provides: frozenset[str] = frozenset({"table9_evidence"})

    def evaluate(self, context, inputs) -> NodeResult:
        variant = inputs["normalized_variant"]
        if not isinstance(variant, NormalizedVariant):
            raise TypeError("normalized_variant must be a NormalizedVariant")

        value = table9_lookup_ps3_bs3(variant.gene, variant.c_notation)
        reviewed = bool(value.get("reviewed"))
        table9_data = load_table9_data()
        dataset_sha256 = table9_checksum()
        evidence = EvidenceItem(
            id="enigma_table9",
            kind="functional_assay",
            status=(
                EvidenceStatus.AVAILABLE
                if reviewed
                else EvidenceStatus.NOT_APPLICABLE
            ),
            value=value,
            source_id="enigma-v1.2-table9",
            source_version=str(table9_data.get("version") or ""),
            source_checksum=dataset_sha256,
            reason=str(value.get("reason") or ""),
            provenance={
                "runtime_dataset": TABLE9_JSON_PATH.name,
                "runtime_rows": table9_data.get("row_count"),
            },
        )
        return NodeResult.succeeded(
            {"table9_evidence": evidence},
            provenance={
                "source_id": evidence.source_id,
                "source_version": evidence.source_version,
                "source_checksum": evidence.source_checksum,
                "reviewed": reviewed,
            },
        )
