"""Paths for mutable application records that are not runtime caches."""

from pathlib import Path
import os


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def choose_runtime_data_dir() -> Path:
    """Return the configured persistent record directory or the local dev path."""
    configured = os.environ.get("ARIANE_RUNTIME_DATA_DIR")
    if configured:
        return Path(configured)

    railway_volume = os.environ.get("RAILWAY_VOLUME_MOUNT_PATH")
    if railway_volume:
        return Path(railway_volume) / "ariane-runtime-data"

    return PROJECT_ROOT / ".runtime-data"


def runtime_data_path(filename: str) -> Path:
    directory = choose_runtime_data_dir()
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename
