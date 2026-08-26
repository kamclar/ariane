"""Paths for mutable runtime caches.

Runtime data must never be mixed with immutable, versioned reference datasets.
"""

from pathlib import Path
import os


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def choose_runtime_cache_dir() -> Path:
    """Return the configured persistent cache directory or the local dev path."""
    configured = os.environ.get("ARIANE_RUNTIME_CACHE_DIR")
    if configured:
        return Path(configured)

    railway_volume = os.environ.get("RAILWAY_VOLUME_MOUNT_PATH")
    if railway_volume:
        return Path(railway_volume) / "ariane-runtime-cache"

    return PROJECT_ROOT / ".runtime-cache"


def runtime_cache_path(filename: str) -> Path:
    """Return a path below the single mutable runtime-cache directory."""
    directory = choose_runtime_cache_dir()
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename
