"""Typed runtime state for population-frequency evidence."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class GnomadRepository:
    """One validated, read-only pair of gnomAD frequency and coverage snapshots."""

    variants: Mapping[str, Sequence[Mapping[str, Any]]]
    metadata: Mapping[str, Any]
    coverage_by_position: Mapping[str, Mapping[str, Any]]
    frequency_path: Path | None
    frequency_status: str
    coverage_path: Path | None
    coverage_status: str

    @property
    def classification_ready(self) -> bool:
        return (
            self.frequency_status == "approved_snapshot"
            and self.coverage_status == "approved_snapshot"
        )


EMPTY_REPOSITORY = GnomadRepository(
    variants={},
    metadata={},
    coverage_by_position={},
    frequency_path=None,
    frequency_status="not_loaded",
    coverage_path=None,
    coverage_status="not_loaded",
)
