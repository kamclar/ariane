"""Explicit application-data initialization and readiness reporting."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import logging
from pathlib import Path
from typing import Any, Callable

from backend.classification_dag.provider_wiring import production_provider_dependencies
from backend.reference_data.paths import (
    ENIGMA_EREPO_VCEP_METADATA_PATH,
    ENIGMA_EREPO_VCEP_REGISTRY_PATH,
    ENIGMA_REFERENCE_TABLES_PATH,
    ENIGMA_RULE_CATALOG_PATH,
    EXON_CNV_EVIDENCE_MANIFEST_PATH,
    EXON_CNV_EVIDENCE_PATH,
    PS1_PROTEIN_REGISTRY_PATH,
    ST2_SPLICE_EVIDENCE_PATH,
    ST7_PATH,
    TABLE4_PATH,
    TABLE9_PATH,
)
from backend.policy.gene import GENE_POLICY_MANIFEST_PATH, GENE_POLICY_METADATA_PATH
from backend.reference_data.validation import validate_required_datasets
from backend.infrastructure.health import DataHealthRegistry
from backend.lookups import bayesdel, coordinates, spliceai
from backend.lookups.founder_variants import load_founder_variant_snapshot
from backend.policy.gene import validate_policy_source_bindings
from backend.population_frequency import PopulationFrequencyService
from backend.reference_data.enigma_rules import validate_rule_catalog
from backend.reference_data.erepo_pvs1_rna import load_erepo_pvs1_rna_registry
from backend.reference_data.pp4_bp5 import load_pp4_bp5_snapshot
from backend.reference_data.residues import initialize_residue_data
from backend.reference_data.table4 import load_table4_data
from backend.reference_data.table9 import load_table9_data
from backend.services import EvidenceOrchestrationService
from backend.policy.spliceai_profile import validate_spliceai_profile
from backend.variant_processing.hgvs_engine import validate_hgvs_engine


LOGGER = logging.getLogger("ariane.startup")


@dataclass(frozen=True)
class StartupCheck:
    component: str
    required: bool
    status: str
    detail: str = ""


class StartupStatus:
    """Collect startup checks without turning data failures into import errors."""

    def __init__(self) -> None:
        self._checks: dict[str, StartupCheck] = {}

    def run(
        self,
        component: str,
        operation: Callable[[], Any],
        *,
        required: bool,
    ) -> Any | None:
        try:
            value = operation()
        except Exception as exc:
            detail = f"{type(exc).__name__}: {exc}"
            self._checks[component] = StartupCheck(
                component=component,
                required=required,
                status="failed",
                detail=detail,
            )
            LOGGER.exception("Startup check failed: %s", component)
            return None
        self._checks[component] = StartupCheck(
            component=component,
            required=required,
            status="ok",
        )
        return value

    @property
    def ready(self) -> bool:
        return not any(
            check.required and check.status != "ok"
            for check in self._checks.values()
        )

    @property
    def degraded(self) -> bool:
        return any(
            not check.required and check.status != "ok"
            for check in self._checks.values()
        )

    def report(self) -> dict[str, Any]:
        failures = [
            asdict(check)
            for check in self._checks.values()
            if check.status != "ok"
        ]
        return {
            "ready": self.ready,
            "status": "ready" if self.ready else "not_ready",
            "checks": [asdict(check) for check in self._checks.values()],
            "failures": failures,
        }


@dataclass(frozen=True)
class ApplicationRuntime:
    status: StartupStatus
    health: DataHealthRegistry
    population_frequency_service: PopulationFrequencyService | None
    classification_orchestration: EvidenceOrchestrationService | None


def _required_dataset_paths() -> dict[str, Path]:
    return {
        "table4": TABLE4_PATH,
        "table9": TABLE9_PATH,
        "enigma_rule_catalog": ENIGMA_RULE_CATALOG_PATH,
        "enigma_reference_tables": ENIGMA_REFERENCE_TABLES_PATH,
        "st7": ST7_PATH,
        "ps1_protein_registry": PS1_PROTEIN_REGISTRY_PATH,
        "enigma_erepo_vcep_registry": ENIGMA_EREPO_VCEP_REGISTRY_PATH,
        "enigma_erepo_vcep_metadata": ENIGMA_EREPO_VCEP_METADATA_PATH,
        "st2_splice_evidence": ST2_SPLICE_EVIDENCE_PATH,
        "exon_cnv_evidence": EXON_CNV_EVIDENCE_PATH,
        "exon_cnv_evidence_manifest": EXON_CNV_EVIDENCE_MANIFEST_PATH,
        "gene_policy_manifest": GENE_POLICY_MANIFEST_PATH,
        "gene_policy_metadata": GENE_POLICY_METADATA_PATH,
    }


def initialize_application_data() -> ApplicationRuntime:
    """Validate immutable inputs and warm runtime sources in one visible place."""
    status = StartupStatus()
    health = DataHealthRegistry()
    status.run(
        "required classification datasets",
        lambda: validate_required_datasets(_required_dataset_paths()),
        required=True,
    )
    status.run("policy source bindings", validate_policy_source_bindings, required=True)
    status.run("ENIGMA Table 4", load_table4_data, required=True)
    status.run("ENIGMA Table 9", load_table9_data, required=True)
    status.run(
        "ERepo PVS1 RNA registry", load_erepo_pvs1_rna_registry, required=True
    )
    status.run("SpliceAI scoring profile", validate_spliceai_profile, required=True)
    status.run("PP4/BP5 clinical LR snapshot", load_pp4_bp5_snapshot, required=True)
    status.run("ENIGMA rule catalog", validate_rule_catalog, required=True)
    status.run("HGVS reference engine", validate_hgvs_engine, required=True)

    # These sources have criterion-level unavailable states. Their failures are
    # visible, but do not make unrelated criteria unsafe.
    status.run(
        "local coordinate sources",
        lambda: coordinates.load_local_coordinate_sources(health),
        required=False,
    )
    status.run(
        "pathogenic founder snapshot",
        lambda: load_founder_variant_snapshot(health=health),
        required=False,
    )
    status.run("BayesDel runtime cache", lambda: bayesdel._load_cache(health), required=False)
    status.run("SpliceAI runtime cache", lambda: spliceai._load_api_cache(health), required=False)
    status.run(
        "clinically important residues",
        lambda: initialize_residue_data(health),
        required=False,
    )

    population_service = status.run(
        "gnomAD population snapshots",
        lambda: PopulationFrequencyService.load_default(health=health),
        required=True,
    )
    orchestration = None
    if isinstance(population_service, PopulationFrequencyService):
        orchestration = EvidenceOrchestrationService(
            provider_dependencies=production_provider_dependencies(
                population_frequency_lookup=population_service.get_frequencies,
                health=health,
            ),
            health=health,
        )
    return ApplicationRuntime(
        status=status,
        health=health,
        population_frequency_service=population_service,
        classification_orchestration=orchestration,
    )
