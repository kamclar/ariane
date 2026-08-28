"""Production runtime for the native ARIANE classification DAG."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import os
from typing import Any, Dict, Mapping, Optional
import uuid

from backend.gene_policy import implementation_profile, policy_version, runtime_policy_id

from backend.classification_dag.engine import DagDefinition, DagExecutor
from backend.classification_dag.domain import (
    EvidenceBundle,
    EvidenceItem,
    EvidenceStatus,
    NormalizedVariant,
    VariantAssertion,
)
from backend.classification_dag.types import (
    DagExecutionContext,
    DagTraceEntry,
    NodeResult,
)
from backend.classification_dag.providers import (
    BayesDelEvidenceNode,
    ClassificationRequest,
    ClassificationRequestContractNode,
    ClinicalLrEvidenceNode,
    CoordinateEvidenceNode,
    EvidenceBundleAssemblyNode,
    ExonCnvEvidenceNode,
    GnomadEvidenceNode,
    ProteinPs1EvidenceNode,
    ProviderDependencies,
    Ps1CandidateEvidenceNode,
    ResidueEvidenceNode,
    SpliceAiEvidenceNode,
    Table9EvidenceNode,
)
from backend.classification_dag.nodes import (
    BioinformaticCriteriaNode,
    ClassificationPolicyNode,
    ClinicalLrCriteriaNode,
    EvidenceInteractionNode,
    ExonCnvCriteriaNode,
    FrequencyCriteriaNode,
    FunctionalCriteriaNode,
    ProteinPs1CriteriaNode,
    Pvs1CriteriaNode,
    ReviewTriageNode,
    SpliceContextNode,
)


class ClassifierEngineMode(str, Enum):
    DAG = "dag"


SUPPORTED_CLASSIFICATION_PROFILES = frozenset({"enigma_brca_vcep_1_2"})


def _require_supported_profile(gene: str) -> str:
    profile = implementation_profile(gene)
    if profile not in SUPPORTED_CLASSIFICATION_PROFILES:
        raise RuntimeError(
            f"No classification DAG is implemented for policy profile {profile!r} "
            f"configured for {gene}"
        )
    return profile


@dataclass(frozen=True)
class ClassificationInputs:
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

    def evidence_bundle(self) -> EvidenceBundle:
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

@dataclass(frozen=True)
class ClassificationExecution:
    result: Dict[str, Any]
    engine_mode: ClassifierEngineMode
    graph_id: str = ""
    graph_version: str = ""
    trace: tuple[DagTraceEntry, ...] = ()
    provider_artifacts: Mapping[str, Any] = field(default_factory=dict)

    def audit_record(self) -> dict[str, Any]:
        return {
            "engine_mode": self.engine_mode.value,
            "graph_id": self.graph_id,
            "graph_version": self.graph_version,
            "provider_evidence": dict(
                self.provider_artifacts.get("evidence_audit", {})
            ),
            "trace": [entry.as_dict() for entry in self.trace],
        }

    @property
    def provider_warnings(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(
            warning
            for entry in self.trace
            if entry.node_id.startswith("provider.")
            for warning in entry.warnings
        ))


@dataclass(frozen=True)
class _ClassificationInputContractNode:
    id: str = "contract.classification_inputs"
    version: str = "1"
    requires: frozenset[str] = frozenset({"classification_inputs"})
    provides: frozenset[str] = frozenset({"normalized_variant", "evidence_bundle"})

    def evaluate(self, context, inputs) -> NodeResult:
        classification_inputs = inputs["classification_inputs"]
        if not isinstance(classification_inputs, ClassificationInputs):
            raise TypeError("classification_inputs must be ClassificationInputs")
        normalized_variant = classification_inputs.normalized_variant()
        if normalized_variant.variant_key != context.variant_key:
            raise ValueError(
                "Normalized variant identity does not match the DAG execution context"
            )
        return NodeResult.succeeded(
            {
                "normalized_variant": normalized_variant,
                "evidence_bundle": classification_inputs.evidence_bundle(),
            },
            provenance={"contract": "ariane.classification-inputs.v1"},
        )


@dataclass(frozen=True)
class _StructuredAssertionNode:
    id: str = "contract.variant_assertion"
    version: str = "1"
    requires: frozenset[str] = frozenset({"unvalidated_classification_result"})
    provides: frozenset[str] = frozenset({"variant_assertion"})

    def evaluate(self, context, inputs) -> NodeResult:
        result = inputs["unvalidated_classification_result"]
        assertion = VariantAssertion.from_public_result(
            result,
            policy_id=context.policy_id,
            policy_version=str(context.metadata["policy_version"]),
        )
        return NodeResult.succeeded(
            {"variant_assertion": assertion},
            provenance={"contract": "ariane.variant-assertion.v1"},
        )


@dataclass(frozen=True)
class _PublicProjectionNode:
    id: str = "projection.public_result"
    version: str = "1"
    requires: frozenset[str] = frozenset({
        "unvalidated_classification_result",
        "variant_assertion",
    })
    provides: frozenset[str] = frozenset({"classification_result"})

    def evaluate(self, context, inputs) -> NodeResult:
        # The typed assertion proves the clinical contract. The established API
        # dictionary is preserved as a presentation boundary.
        assertion = inputs["variant_assertion"]
        if not isinstance(assertion, VariantAssertion):
            raise TypeError("variant_assertion must be a VariantAssertion")
        return NodeResult.succeeded(
            {"classification_result": inputs["unvalidated_classification_result"]},
            provenance={"projection": "public-result-v1"},
        )


def build_native_graph() -> DagDefinition:
    """Build the synchronous classification graph from explicit evidence inputs."""
    return DagDefinition(
        id="ariane.vcep.classification",
        version="4.0.0-gene-policy",
        seed_keys={"classification_inputs"},
        nodes=(
            _ClassificationInputContractNode(),
            Table9EvidenceNode(),
            SpliceContextNode(),
            FrequencyCriteriaNode(),
            ExonCnvCriteriaNode(),
            FunctionalCriteriaNode(),
            Pvs1CriteriaNode(),
            ClinicalLrCriteriaNode(),
            ProteinPs1CriteriaNode(),
            BioinformaticCriteriaNode(),
            EvidenceInteractionNode(),
            ClassificationPolicyNode(),
            ReviewTriageNode(),
            _StructuredAssertionNode(),
            _PublicProjectionNode(),
        ),
    )


def build_provider_graph(
    dependencies: ProviderDependencies | None = None,
) -> DagDefinition:
    """Build the production graph including all classification evidence providers."""
    providers = dependencies or ProviderDependencies.production()
    return DagDefinition(
        id="ariane.vcep.classification",
        version="4.0.0-gene-policy-provider-dag",
        seed_keys={"classification_request"},
        nodes=(
            ClassificationRequestContractNode(),
            CoordinateEvidenceNode(providers),
            Table9EvidenceNode(),
            Ps1CandidateEvidenceNode(providers),
            SpliceAiEvidenceNode(providers),
            BayesDelEvidenceNode(providers),
            GnomadEvidenceNode(providers),
            ClinicalLrEvidenceNode(providers),
            ExonCnvEvidenceNode(providers),
            ResidueEvidenceNode(providers),
            ProteinPs1EvidenceNode(providers),
            EvidenceBundleAssemblyNode(),
            SpliceContextNode(),
            FrequencyCriteriaNode(),
            ExonCnvCriteriaNode(),
            FunctionalCriteriaNode(),
            Pvs1CriteriaNode(),
            ClinicalLrCriteriaNode(),
            ProteinPs1CriteriaNode(),
            BioinformaticCriteriaNode(),
            EvidenceInteractionNode(),
            ClassificationPolicyNode(),
            ReviewTriageNode(),
            _StructuredAssertionNode(),
            _PublicProjectionNode(),
        ),
    )


def get_configured_engine_mode(value: str | None = None) -> ClassifierEngineMode:
    raw = (value if value is not None else os.getenv(
        "ARIANE_CLASSIFIER_ENGINE", "dag"
    )).strip().lower()
    try:
        return ClassifierEngineMode(raw)
    except ValueError as exc:
        allowed = ", ".join(mode.value for mode in ClassifierEngineMode)
        raise ValueError(
            f"Invalid ARIANE_CLASSIFIER_ENGINE={raw!r}; expected one of {allowed}"
        ) from exc


def execute_classification(
    inputs: ClassificationInputs,
    *,
    mode: ClassifierEngineMode | str | None = None,
) -> ClassificationExecution:
    selected_mode = (
        mode if isinstance(mode, ClassifierEngineMode)
        else get_configured_engine_mode(mode)
    )
    _require_supported_profile(inputs.gene)
    graph = build_native_graph()
    dag_run = DagExecutor(graph).run(
        {"classification_inputs": inputs},
        context=DagExecutionContext(
            run_id=uuid.uuid4().hex,
            variant_key=inputs.variant_key,
            policy_id=runtime_policy_id(inputs.gene),
            metadata={
                "engine_mode": selected_mode.value,
                "policy_version": policy_version(inputs.gene),
            },
        ),
    )
    dag_result = dag_run.values.get("classification_result")
    if not isinstance(dag_result, dict):
        raise RuntimeError("Classification DAG did not produce classification_result")

    return ClassificationExecution(
        result=dag_result,
        engine_mode=selected_mode,
        graph_id=dag_run.graph_id,
        graph_version=dag_run.graph_version,
        trace=dag_run.trace,
        provider_artifacts={},
    )


async def execute_classification_request(
    request: ClassificationRequest,
    *,
    dependencies: ProviderDependencies | None = None,
    mode: ClassifierEngineMode | str | None = None,
) -> ClassificationExecution:
    """Acquire evidence and classify one normalized variant in one DAG run."""
    selected_mode = (
        mode if isinstance(mode, ClassifierEngineMode)
        else get_configured_engine_mode(mode)
    )
    _require_supported_profile(request.variant.gene)
    graph = build_provider_graph(dependencies)
    dag_run = await DagExecutor(graph).run_async(
        {"classification_request": request},
        context=DagExecutionContext(
            run_id=uuid.uuid4().hex,
            variant_key=request.variant.variant_key,
            policy_id=runtime_policy_id(request.variant.gene),
            metadata={
                "engine_mode": selected_mode.value,
                "policy_version": policy_version(request.variant.gene),
            },
        ),
    )
    dag_result = dag_run.values.get("classification_result")
    if not isinstance(dag_result, dict):
        raise RuntimeError("Classification provider DAG did not produce classification_result")
    artifacts = dag_run.values.get("provider_artifacts")
    if not isinstance(artifacts, Mapping):
        raise RuntimeError("Classification provider DAG did not produce provider_artifacts")
    return ClassificationExecution(
        result=dag_result,
        engine_mode=selected_mode,
        graph_id=dag_run.graph_id,
        graph_version=dag_run.graph_version,
        trace=dag_run.trace,
        provider_artifacts=dict(artifacts),
    )


def compare_classification_results(
    expected: Any,
    actual: Any,
    *,
    path: str = "$",
) -> tuple[str, ...]:
    differences: list[str] = []
    _compare_values(expected, actual, path, differences)
    return tuple(differences)


def _compare_values(
    expected: Any,
    actual: Any,
    path: str,
    differences: list[str],
) -> None:
    if len(differences) >= 50:
        return
    if type(expected) is not type(actual):
        differences.append(
            f"{path}: type {type(expected).__name__} != {type(actual).__name__}"
        )
        return
    if isinstance(expected, Mapping):
        expected_keys = set(expected)
        actual_keys = set(actual)
        for key in sorted(expected_keys - actual_keys, key=str):
            differences.append(f"{path}.{key}: missing from DAG result")
        for key in sorted(actual_keys - expected_keys, key=str):
            differences.append(f"{path}.{key}: only in DAG result")
        for key in sorted(expected_keys & actual_keys, key=str):
            _compare_values(
                expected[key], actual[key], f"{path}.{key}", differences
            )
        return
    if isinstance(expected, (list, tuple)):
        if len(expected) != len(actual):
            differences.append(
                f"{path}: length {len(expected)} != {len(actual)}"
            )
            return
        for index, (expected_item, actual_item) in enumerate(zip(expected, actual)):
            _compare_values(
                expected_item,
                actual_item,
                f"{path}[{index}]",
                differences,
            )
        return
    if expected != actual:
        differences.append(f"{path}: {expected!r} != {actual!r}")
