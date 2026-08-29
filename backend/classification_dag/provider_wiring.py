"""Production composition root for classification evidence providers.

This module is the only place where the abstract provider dependency contract
is connected to concrete lookup, cache, and local dataset implementations.
Keeping the wiring explicit makes the code dependency graph statically visible
without coupling individual DAG nodes to production infrastructure.
"""

from __future__ import annotations

from typing import Any, Callable

from backend.classification_dag.providers import ProviderDependencies
from backend.lookups import bayesdel, coordinates, spliceai
from backend.modules import (
    exon_cnv_evidence,
    pp4_bp5,
    ps1,
    ps1_splice_evidence,
    residues,
)
from backend.population_frequency import PopulationFrequencyService


def production_provider_dependencies(
    *,
    population_frequency_lookup: Callable[..., Any] | None = None,
) -> ProviderDependencies:
    """Bind the provider interface to the production evidence implementations."""
    if population_frequency_lookup is None:
        population_frequency_lookup = (
            PopulationFrequencyService.load_default().get_frequencies
        )
    return ProviderDependencies(
        resolve_variant=lambda *args, **kwargs: coordinates.resolve_variant(
            *args, **kwargs
        ),
        get_grch37=lambda *args, **kwargs: coordinates.get_grch37(*args, **kwargs),
        get_grch38=lambda *args, **kwargs: coordinates.get_grch38(*args, **kwargs),
        spliceai_lookup=lambda *args, **kwargs: spliceai.get_spliceai_score(
            *args, **kwargs
        ),
        spliceai_status=lambda gene, c: dict(
            spliceai.SPLICEAI_STATUS_CACHE.get(f"{gene}:{c}", {})
        ),
        bayesdel_lookup=lambda *args, **kwargs: bayesdel.get_bayesdel_and_alphamissense(
            *args, **kwargs
        ),
        bayesdel_status=lambda gene, c: dict(
            bayesdel.BAYESDEL_STATUS_CACHE.get(f"{gene}:{c}", {})
        ),
        gnomad_lookup=population_frequency_lookup,
        clinical_lr_lookup=lambda *args, **kwargs: pp4_bp5.evaluate_pp4_bp5(
            *args, **kwargs
        ),
        exon_cnv_lookup=lambda *args, **kwargs: exon_cnv_evidence.lookup_exon_cnv_evidence(
            *args, **kwargs
        ),
        residue_lookup=lambda *args, **kwargs: residues.check_important_residue(
            *args, **kwargs
        ),
        ps1_candidate_lookup=lambda *args, **kwargs: ps1.discover_ps1_reference_variants(
            *args, **kwargs
        ),
        splice_source_lookup=lambda *args, **kwargs: ps1_splice_evidence.evaluate_defined_splice_sources(
            *args, **kwargs
        ),
        select_ps1_spliceai=lambda *args, **kwargs: ps1.select_vua_spliceai_for_ps1(
            *args, **kwargs
        ),
        ps1_lookup=lambda *args, **kwargs: ps1.evaluate_ps1(*args, **kwargs),
    )
