"""Registry of mutable runtime caches owned by higher-level components."""

from __future__ import annotations

from collections.abc import Callable


CacheClearer = Callable[[], None]
_CACHE_CLEARERS: dict[str, CacheClearer] = {}


def register_runtime_cache(name: str, clearer: CacheClearer) -> None:
    """Register or replace one named cache clearer without importing its owner."""
    normalized = name.strip()
    if not normalized:
        raise ValueError("Runtime cache registration requires a non-empty name")
    _CACHE_CLEARERS[normalized] = clearer


def registered_runtime_caches() -> tuple[str, ...]:
    """Expose deterministic registration state for startup diagnostics and tests."""
    return tuple(sorted(_CACHE_CLEARERS))


def clear_runtime_caches() -> None:
    """Clear every cache registered by an imported application component."""
    for name in registered_runtime_caches():
        _CACHE_CLEARERS[name]()
