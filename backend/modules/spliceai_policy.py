"""ENIGMA BRCA v1.2 policy for SpliceAI provenance comparisons.

Specifications Table 9 records the SpliceAI context used when ENIGMA reviewed
published functional evidence for PS3 and BS3.  It is not a replacement source
for the configured SpliceAI result used by Figure 1A and protein PS1 rules.
"""

from __future__ import annotations

from typing import Any, Mapping

from backend.gene_policy import spliceai_thresholds


def enigma_spliceai_band(gene: str, score: float | None) -> str:
    """Return the configured VCEP prediction band for a SpliceAI score."""
    if score is None:
        return "unavailable"
    thresholds = spliceai_thresholds(gene)
    if score <= thresholds["bp4"]:
        return "no_impact"
    if score < thresholds["pp3"]:
        return "not_informative"
    return "impact"


def compare_table9_spliceai(
    gene: str,
    configured_score: float | None,
    table9_result: Mapping[str, Any] | None,
) -> tuple[float | None, tuple[str, ...]]:
    """Compare Table 9 context without overriding the configured score.

    The returned Table 9 value is audit information only.  If the configured
    score is unavailable, callers must keep SpliceAI-dependent rules unavailable.
    """
    table9_result = table9_result or {}
    if not table9_result.get("reviewed"):
        return None, ()

    warnings: list[str] = []
    raw_score = table9_result.get("spliceai_prediction")
    table9_score = float(raw_score) if isinstance(raw_score, (int, float)) else None

    if table9_score is not None and configured_score is None:
        warnings.append(
            "The configured SpliceAI result is unavailable. ENIGMA Table 9 records "
            f"SpliceAI {table9_score:.3f} as context for its PS3/BS3 review, but this "
            "value does not replace the missing prediction. SpliceAI-dependent "
            "automated criteria were not evaluated."
        )
    elif (
        table9_score is not None
        and configured_score is not None
        and abs(configured_score - table9_score) > 1e-9
    ):
        threshold_crossed = (
            enigma_spliceai_band(gene, configured_score)
            != enigma_spliceai_band(gene, table9_score)
        )
        warning = (
            f"SpliceAI differs from ENIGMA Table 9: configured source="
            f"{configured_score:.3f}, Table 9={table9_score:.3f}. Automated "
            "bioinformatic criteria use the configured ENIGMA-compatible SpliceAI "
            "result. The Table 9 value is retained only as audit context for the "
            "PS3/BS3 review."
        )
        if threshold_crossed:
            warning += (
                " The values fall in different ENIGMA prediction bands; expert "
                "review of SpliceAI provenance is required."
            )
        warnings.append(warning)

    splice_result = table9_result.get("splice_result_published")
    if splice_result:
        warnings.append(f"ENIGMA Table 9 published splice result: {splice_result}.")

    return table9_score, tuple(warnings)
