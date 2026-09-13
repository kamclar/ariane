"""ARIANE FastAPI composition root."""

from __future__ import annotations

import logging
from pathlib import Path

from fastapi import FastAPI, HTTPException

from backend.api.admin import router as admin_router
from backend.api.audit import RequestAuditLogger
from backend.api.classification import ClassificationApi
from backend.api.errors import install_exception_handlers
from backend.api.frontend import install_frontend
from backend.api.manual import create_manual_router
from backend.api.middleware import install_http_middleware
from backend.api.review import router as review_router
from backend.api.system import create_system_router
from backend.bootstrap import initialize_application_data
from backend.classification_runtime import (
    ClassificationCacheRepository,
    ClassificationUsageRepository,
    public_api_daily_classification_limit,
)
from backend.presentation.frontend import load_frontend_assets
from backend.services import (
    EvidenceOrchestrationService,
    ManualEvidenceService,
    VariantClassificationService,
)
from backend.version import ARIANE_VERSION


APPLICATION_RUNTIME = initialize_application_data()
STARTUP_STATUS = APPLICATION_RUNTIME.status
POPULATION_FREQUENCY_SERVICE = APPLICATION_RUNTIME.population_frequency_service
CLASSIFICATION_ORCHESTRATION = APPLICATION_RUNTIME.classification_orchestration
CLASSIFICATION_CACHE: ClassificationCacheRepository | None = None
CLASSIFICATION_USAGE: ClassificationUsageRepository | None = None
PUBLIC_API_QUOTA: ClassificationUsageRepository | None = None
PUBLIC_API_DAILY_CLASSIFICATION_LIMIT = public_api_daily_classification_limit()
BATCH_CONCURRENCY = 3


def _optional_classification_cache() -> ClassificationCacheRepository | None:
    try:
        return ClassificationCacheRepository()
    except Exception:
        logging.getLogger("ariane.classification_cache").exception(
            "Classification cache could not be initialized"
        )
        return None


def _optional_usage_repository() -> ClassificationUsageRepository | None:
    try:
        return ClassificationUsageRepository()
    except Exception:
        logging.getLogger("ariane.classification_usage").exception(
            "Classification usage storage could not be initialized"
        )
        return None


CLASSIFICATION_CACHE = _optional_classification_cache()
CLASSIFICATION_USAGE = _optional_usage_repository()
PUBLIC_API_QUOTA = CLASSIFICATION_USAGE
CLASSIFICATION_SERVICE = VariantClassificationService(
    batch_concurrency=BATCH_CONCURRENCY,
)
AUDIT = RequestAuditLogger.configured()
_audit = AUDIT


def _classification_orchestration() -> EvidenceOrchestrationService:
    if not STARTUP_STATUS.ready or CLASSIFICATION_ORCHESTRATION is None:
        raise HTTPException(
            status_code=503,
            detail=(
                "ARIANE classification is unavailable because required startup data "
                "did not pass validation. See /api/health for the failing component."
            ),
            headers={
                "X-ARIANE-Error-Code": "application_not_ready",
                "X-ARIANE-Retryable": "false",
            },
        )
    return CLASSIFICATION_ORCHESTRATION


CLASSIFICATION_API = ClassificationApi(
    service=CLASSIFICATION_SERVICE,
    orchestration=_classification_orchestration,
    cache_repository=lambda: CLASSIFICATION_CACHE,
    usage_repository=lambda: CLASSIFICATION_USAGE,
    quota_repository=lambda: PUBLIC_API_QUOTA,
    daily_classification_limit=lambda: PUBLIC_API_DAILY_CLASSIFICATION_LIMIT,
    audit=AUDIT,
)

app = FastAPI(
    title="ARIANE",
    description="Automated ACMG Rule-based Interpretation and Annotation ENgine",
    version=ARIANE_VERSION,
)
install_http_middleware(app, audit=AUDIT)
install_exception_handlers(app, audit=AUDIT)

app.include_router(admin_router)
app.include_router(review_router)
app.include_router(CLASSIFICATION_API.router())
app.include_router(
    create_manual_router(
        service=ManualEvidenceService(),
        audit=AUDIT,
    )
)
app.include_router(
    create_system_router(
        runtime=APPLICATION_RUNTIME,
        population_frequency_service=POPULATION_FREQUENCY_SERVICE,
        audit=AUDIT,
    )
)

FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"
FRONTEND_ASSETS = load_frontend_assets(FRONTEND_DIR, ARIANE_VERSION)
INDEX_HTML_TEMPLATE = FRONTEND_ASSETS.html
FRONTEND_ASSET_VERSION = FRONTEND_ASSETS.version
install_frontend(app, frontend_dir=FRONTEND_DIR, assets=FRONTEND_ASSETS)
