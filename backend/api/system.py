"""Health, reference-resource, and maintenance HTTP routes."""

from __future__ import annotations

from collections.abc import Callable
import os
import secrets
from typing import Any, Protocol

from fastapi import APIRouter, Header, HTTPException, Query, Request
from fastapi.responses import JSONResponse

from backend.api.routing import PairedApiRoutes
from backend.domain.classification import CLASSIFICATION_ENGINE_ID
from backend.infrastructure.cache_registry import clear_runtime_caches
from backend.lookups.spliceai import spliceai_runtime_health
from backend.policy.gene import active_genes, get_gene_policy, not_used_criteria
from backend.population_frequency import PopulationFrequencyService
from backend.reference_data.enigma_rules import (
    get_decision_tree,
    public_catalog,
    search_reference_table,
    search_table9,
)
from backend.reference_data.paths import (
    ENIGMA_EREPO_VCEP_REGISTRY_PATH,
    ENIGMA_RULE_CATALOG_PATH,
    EXON_CNV_EVIDENCE_PATH,
    PS1_PROTEIN_REGISTRY_PATH,
    ST2_SPLICE_EVIDENCE_PATH,
    ST7_PATH,
    TABLE4_PATH,
    TABLE9_PATH,
)
from backend.reference_data.ps1_splice_evidence import list_splice_ps1_candidate_discovery
from backend.review.definitions import manual_criteria_for_gene, resource_links_for_gene
from backend.variant_processing.hgvs_provider import load_panel_provider
from backend.version import ARIANE_VERSION


AuditWriter = Callable[..., None]


class RuntimeHealth(Protocol):
    def issues(self) -> list[dict[str, str]]: ...


class RuntimeStartupStatus(Protocol):
    degraded: bool

    def report(self) -> dict[str, Any]: ...


class ApplicationRuntimeView(Protocol):
    status: RuntimeStartupStatus
    health: RuntimeHealth


def create_system_router(
    *,
    runtime: ApplicationRuntimeView,
    population_frequency_service: PopulationFrequencyService | None,
    audit: AuditWriter,
) -> APIRouter:
    router = APIRouter()
    paired = PairedApiRoutes(router)

    @router.get("/api/health")
    async def health():
        issues = runtime.health.issues()
        startup = runtime.status.report()
        panel_provenance: dict[str, str] = {}
        try:
            panel_provenance = load_panel_provider().provenance
        except Exception:
            pass
        spliceai = spliceai_runtime_health()
        status = (
            "not_ready"
            if not startup["ready"]
            else "degraded"
            if issues
            or runtime.status.degraded
            or (spliceai["local"] and spliceai["status"] != "ok")
            else "ok"
        )
        payload = {
            "status": status,
            "ready": startup["ready"],
            "version": ARIANE_VERSION,
            "classification_engine": CLASSIFICATION_ENGINE_ID,
            "startup": startup,
            "data": {
                "table4": TABLE4_PATH.exists(),
                "table9": TABLE9_PATH.exists(),
                "enigma_rule_catalog": ENIGMA_RULE_CATALOG_PATH.exists(),
                "st7": ST7_PATH.exists(),
                "ps1_protein_registry": PS1_PROTEIN_REGISTRY_PATH.exists(),
                "enigma_erepo_vcep_registry": ENIGMA_EREPO_VCEP_REGISTRY_PATH.exists(),
                "st2_splice_evidence": ST2_SPLICE_EVIDENCE_PATH.exists(),
                "exon_cnv_evidence": EXON_CNV_EVIDENCE_PATH.exists(),
                "reference_bundle": panel_provenance.get("reference_bundle", ""),
                "normalization_engine": panel_provenance.get("normalization_engine", ""),
                "spliceai": spliceai,
            },
            "data_issues": issues,
        }
        return JSONResponse(status_code=200 if startup["ready"] else 503, content=payload)

    @paired.get("/resources")
    async def resources(gene: str | None = None):
        return {
            "version": ARIANE_VERSION,
            "build_revision": os.getenv("ARIANE_BUILD_REVISION", "").strip(),
            "issue_tracker_url": os.getenv("ARIANE_ISSUE_TRACKER_URL", "").strip(),
            "manual_criteria": manual_criteria_for_gene(gene) if gene else {},
            "genes": [
                {
                    "symbol": symbol,
                    "reference_transcript": get_gene_policy(symbol)["gene_config"]["reference_transcript"],
                    "reference_protein": get_gene_policy(symbol)["gene_config"]["reference_protein"],
                    "policy_id": get_gene_policy(symbol)["policy"]["runtime_policy_id"],
                    "policy_name": get_gene_policy(symbol)["policy"]["name"],
                    "policy_version": get_gene_policy(symbol)["policy"]["version"],
                    "policy_source_url": get_gene_policy(symbol)["policy"]["source_url"],
                    "not_used_criteria": not_used_criteria(symbol),
                }
                for symbol in active_genes()
            ],
            "links": resource_links_for_gene(gene),
            "splice_ps1_candidates": list_splice_ps1_candidate_discovery(),
        }

    @paired.get("/rules")
    async def enigma_rules_catalog():
        return public_catalog()

    @paired.get("/rules/trees/{tree_id}")
    async def enigma_decision_tree(tree_id: str):
        tree = get_decision_tree(tree_id)
        if tree is None:
            raise HTTPException(status_code=404, detail="Decision tree not found")
        return tree

    @paired.get("/rules/tables/table9")
    async def enigma_table9_records(
        gene: str | None = Query(default=None, max_length=40),
        query: str = Query(default="", max_length=200),
        code: str | None = Query(default=None, max_length=20),
        page: int = Query(default=1, ge=1),
        page_size: int = Query(default=25, ge=1, le=100),
    ):
        if gene is not None:
            gene = gene.strip().upper()
            if gene not in set(active_genes()):
                raise HTTPException(
                    status_code=422,
                    detail=f"Gene must be one of: {', '.join(active_genes())}",
                )
        return search_table9(gene=gene, query=query, code=code, page=page, page_size=page_size)

    @paired.get("/rules/tables/{table_id}")
    async def enigma_reference_table_records(
        table_id: str,
        section: str | None = Query(default=None, max_length=80),
        query: str = Query(default="", max_length=200),
        page: int = Query(default=1, ge=1),
        page_size: int = Query(default=25, ge=1, le=100),
    ):
        payload = search_reference_table(
            table_id,
            section_id=section,
            query=query,
            page=page,
            page_size=page_size,
        )
        if payload is None:
            raise HTTPException(status_code=404, detail="ENIGMA table or section not found")
        return payload

    @router.post("/api/clear-cache")
    async def clear_cache(
        request: Request,
        x_ariane_admin_token: str | None = Header(default=None),
    ):
        admin_token = os.getenv("ARIANE_ADMIN_TOKEN", "")
        if not admin_token:
            raise HTTPException(status_code=503, detail="Administrative API is disabled")
        if not x_ariane_admin_token or not secrets.compare_digest(
            x_ariane_admin_token,
            admin_token,
        ):
            raise HTTPException(status_code=403, detail="Invalid administrative token")
        clear_runtime_caches()
        if population_frequency_service is None:
            raise HTTPException(
                status_code=503,
                detail="Population frequency service is not initialized",
            )
        population_frequency_service.reload()
        audit(request, "cache_cleared")
        return {"status": "ok", "message": "All caches cleared"}

    return router
