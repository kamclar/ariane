"""Explicit registry of mutable runtime caches owned by the application."""

from backend.lookups.bayesdel import BAYESDEL_CACHE
from backend.lookups.clingen import EREPO_CACHE
from backend.lookups.clinvar import CLINVAR_CACHE
from backend.lookups.spliceai import SPLICEAI_CACHE, SPLICEAI_STATUS_CACHE


def clear_runtime_caches() -> None:
    """Clear API-derived in-memory caches without touching immutable datasets."""
    SPLICEAI_CACHE.clear()
    SPLICEAI_STATUS_CACHE.clear()
    BAYESDEL_CACHE.clear()
    CLINVAR_CACHE.clear()
    EREPO_CACHE.clear()
