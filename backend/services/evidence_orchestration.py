"""Normalize one variant, acquire evidence, and execute the classification DAG."""

from __future__ import annotations

import asyncio
from dataclasses import dataclass
import json
import logging
from typing import Any, Callable, Mapping

from backend.classification_dag import (
    ClassificationRequest,
    ClassifierEngineMode,
    DagNodeExecutionError,
    NormalizedVariant,
    ProviderDependencies,
    execute_classification_request,
)
from backend.classification_dag.runtime import ClassificationExecution
from backend.data_health import get_user_warnings
from backend.lookup_execution import lookup_or_unavailable
from backend.modules.variant_input import NormalizedVariantInput, normalize_variant_input
from backend.modules.variant_type import infer_variant_type


LOGGER = logging.getLogger("ariane.evidence_orchestration")


class VariantPreparationError(ValueError):
    """The submitted variant cannot enter evidence acquisition."""


class EvidenceExecutionError(RuntimeError):
    """The evidence or classification DAG failed closed."""

    def __init__(self, node_id: str, trace: tuple[Any, ...]):
        self.node_id = node_id
        self.trace = trace
        super().__init__(
            f"Classification could not complete at internal step {node_id}. "
            "No classification was returned; the failure was recorded for review."
        )


@dataclass(frozen=True)
class ClassificationCommand:
    gene: str
    c_notation: str
    p_notation: str = ""
    dup_type: str = "Unknown"


@dataclass(frozen=True)
class ExternalEvidenceDependencies:
    clinvar_lookup: Callable[..., Any]
    clingen_lookup: Callable[..., Any]

    @classmethod
    def production(cls) -> "ExternalEvidenceDependencies":
        from backend.lookups.clingen import clingen_erepo_lookup
        from backend.lookups.clinvar import clinvar_lookup

        return cls(
            clinvar_lookup=clinvar_lookup,
            clingen_lookup=clingen_erepo_lookup,
        )


@dataclass(frozen=True)
class OrchestratedEvidence:
    command: ClassificationCommand
    normalized_input: NormalizedVariantInput
    variant: NormalizedVariant
    variant_type: str
    is_exon_cnv: bool
    execution: ClassificationExecution
    result: dict[str, Any]
    artifacts: Mapping[str, Any]
    clinvar: Mapping[str, Any]
    clingen: Mapping[str, Any]
    external_diagnostics: tuple[str, ...]


class EvidenceOrchestrationService:
    """Coordinate normalization, evidence providers, DAG execution and diagnostics."""

    def __init__(
        self,
        *,
        engine_mode: ClassifierEngineMode,
        provider_dependencies: ProviderDependencies | None = None,
        external_dependencies: ExternalEvidenceDependencies | None = None,
    ) -> None:
        self.engine_mode = engine_mode
        self.provider_dependencies = provider_dependencies
        self.external_dependencies = external_dependencies

    def _external_dependencies(self) -> ExternalEvidenceDependencies:
        return self.external_dependencies or ExternalEvidenceDependencies.production()

    @staticmethod
    def prepare(
        command: ClassificationCommand,
    ) -> tuple[NormalizedVariantInput, NormalizedVariant]:
        try:
            normalized = normalize_variant_input(
                command.gene,
                command.c_notation,
                p_notation=command.p_notation or None,
            )
        except (ValueError, RuntimeError) as exc:
            raise VariantPreparationError(str(exc)) from exc
        variant_type = infer_variant_type(
            normalized.c_notation,
            normalized.p_notation,
        )
        variant = NormalizedVariant(
            gene=normalized.gene,
            reference_transcript=normalized.reference_transcript,
            c_notation=normalized.c_notation,
            p_notation=normalized.p_notation,
            variant_type=variant_type,
            submitted_notation=normalized.submitted_notation,
            normalization_source=normalized.normalization_source,
            consequence_status=normalized.consequence_status,
            normalization_provenance=normalized.normalization_provenance or {},
            protein_consequence_explanation=normalized.protein_consequence_explanation,
            assembly=normalized.assembly,
            genomic_notation=normalized.genomic_notation,
        )
        return normalized, variant

    async def orchestrate(self, command: ClassificationCommand) -> OrchestratedEvidence:
        normalized, variant = self.prepare(command)
        is_exon_cnv = variant.variant_type.lower() in {
            "exon_deletion",
            "exon_duplication",
        }
        request = ClassificationRequest(variant=variant, dup_type=command.dup_type)
        external_dependencies = self._external_dependencies()
        external_diagnostics: list[str] = []
        classification_task = execute_classification_request(
            request,
            dependencies=self.provider_dependencies,
            mode=self.engine_mode,
        )
        external_tasks = (
            lookup_or_unavailable(
                external_dependencies.clinvar_lookup,
                {"status": "api_timeout", "error": "ClinVar lookup timed out"},
                "ClinVar",
                external_diagnostics,
                variant.gene,
                variant.c_notation,
            ),
            lookup_or_unavailable(
                external_dependencies.clingen_lookup,
                {
                    "status": "api_timeout",
                    "error": "ClinGen ERepo lookup timed out",
                },
                "ClinGen ERepo",
                external_diagnostics,
                variant.gene,
                variant.c_notation,
            ),
        )
        try:
            execution, clinvar, clingen = await asyncio.gather(
                classification_task,
                *external_tasks,
            )
        except DagNodeExecutionError as exc:
            LOGGER.exception(
                "Classification DAG failed for %s at %s: %s; trace=%s",
                variant.variant_key,
                exc.node_id,
                exc,
                json.dumps(
                    [entry.as_dict() for entry in exc.trace],
                    ensure_ascii=True,
                    separators=(",", ":"),
                ),
            )
            raise EvidenceExecutionError(exc.node_id, exc.trace) from exc

        artifacts = dict(execution.provider_artifacts)
        result = dict(execution.result)
        result["warnings"] = list(result.get("warnings", []))
        self._record_diagnostics(variant, execution, external_diagnostics)
        self._append_availability_warnings(
            result,
            artifacts,
            is_exon_cnv=is_exon_cnv,
            clinvar=clinvar,
            clingen=clingen,
        )
        return OrchestratedEvidence(
            command=command,
            normalized_input=normalized,
            variant=variant,
            variant_type=variant.variant_type,
            is_exon_cnv=is_exon_cnv,
            execution=execution,
            result=result,
            artifacts=artifacts,
            clinvar=dict(clinvar),
            clingen=dict(clingen),
            external_diagnostics=tuple(external_diagnostics),
        )

    @staticmethod
    def _record_diagnostics(
        variant: NormalizedVariant,
        execution: ClassificationExecution,
        external_diagnostics: list[str],
    ) -> None:
        if execution.trace:
            LOGGER.info(
                "Classification DAG execution: %s",
                json.dumps(
                    {
                        "variant": variant.variant_key,
                        **execution.audit_record(),
                    },
                    ensure_ascii=True,
                    separators=(",", ":"),
                ),
            )
        for diagnostic in (*execution.provider_warnings, *external_diagnostics):
            LOGGER.warning("External lookup diagnostic: %s", diagnostic)

    @staticmethod
    def _append_availability_warnings(
        result: dict[str, Any],
        artifacts: Mapping[str, Any],
        *,
        is_exon_cnv: bool,
        clinvar: Mapping[str, Any],
        clingen: Mapping[str, Any],
    ) -> None:
        warnings = result["warnings"]
        if is_exon_cnv:
            warnings.append(
                "Small-variant coordinate-dependent evidence was not evaluated because "
                "this exon-level copy-number variant has uncertain genomic breakpoints. "
                "Population evidence was evaluated through the general ENIGMA Appendix G "
                "exon-CNV decision path."
            )
            exon_cnv = artifacts.get("exon_cnv_result") or {}
            if exon_cnv.get("reason"):
                warnings.append(exon_cnv["reason"])
            table9 = artifacts.get("table9_result") or {}
            if not table9.get("applies"):
                warnings.append(
                    "PS3/BS3 was not applied because this variant has no applicable "
                    "variant-specific recommendation in ENIGMA Table 9."
                )
        elif not artifacts.get("grch37") and not artifacts.get("grch38"):
            warnings.append(
                "Coordinate-dependent evidence was not evaluated because genomic "
                "coordinates could not be resolved."
            )
        for warning in get_user_warnings():
            if warning not in warnings:
                warnings.append(warning)
        if clinvar.get("status") == "ambiguous":
            warnings.append(
                "ClinVar lookup was ambiguous; no external ClinVar record was selected. "
                f"Candidate IDs: {', '.join(clinvar.get('candidate_ids', [])) or 'not reported'}."
            )
        elif clinvar.get("status") not in {"ok", "not_found"}:
            warnings.append("ClinVar comparison is temporarily unavailable.")
        if clingen.get("status") not in {"ok", "not_found"}:
            warnings.append("ClinGen ERepo comparison is temporarily unavailable.")
