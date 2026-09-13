"""Application-owned population-frequency evidence service."""

from __future__ import annotations

from threading import RLock
from typing import Any, Callable, Mapping

from backend.lookups.founder_variants import lookup_pathogenic_founder_variant
from backend.population_frequency.lookup import get_gnomad_frequencies
from backend.population_frequency.models import GnomadRepository
from backend.population_frequency.snapshot_repository import load_gnomad_repository
from backend.infrastructure.health import DataHealthRegistry


FounderLookup = Callable[[str, str], Mapping[str, Any]]


class PopulationFrequencyService:
    """Own one immutable snapshot binding and replace it atomically on reload."""

    def __init__(
        self,
        repository: GnomadRepository,
        *,
        founder_lookup: FounderLookup = lookup_pathogenic_founder_variant,
        health: DataHealthRegistry | None = None,
    ) -> None:
        self._repository = repository
        self._founder_lookup = founder_lookup
        self._health = health
        self._lock = RLock()

    @classmethod
    def load_default(
        cls,
        *,
        health: DataHealthRegistry | None = None,
    ) -> "PopulationFrequencyService":
        return cls(load_gnomad_repository(health=health), health=health)

    @property
    def repository(self) -> GnomadRepository:
        with self._lock:
            return self._repository

    def reload(self) -> GnomadRepository:
        repository = load_gnomad_repository(health=self._health)
        with self._lock:
            self._repository = repository
        return repository

    def get_frequencies(
        self,
        *,
        gene: str | None,
        c_notation: str = "",
        grch37: Any | None = None,
        grch38: Any | None = None,
    ) -> dict[str, Any]:
        repository = self.repository
        result = get_gnomad_frequencies(
            repository,
            gene=gene,
            grch37=grch37,
            grch38=grch38,
        )
        founder_exception = dict(
            self._founder_lookup(gene or "", c_notation or "")
        )
        result["founder_exception"] = founder_exception
        audit = result.get("population_frequency_audit")
        if isinstance(audit, dict):
            audit["founder_exception"] = founder_exception
        return result
