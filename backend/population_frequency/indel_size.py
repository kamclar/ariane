"""Population-frequency adapter for shared HGVS indel-size assessment."""

from __future__ import annotations

from typing import Any

from backend.domain.indel import (
    assess_indel_size_for_transcript,
    indel_operation as indel_operation,
    is_indel_allele as is_indel_allele,
)
from backend.policy.catalog import TRANSCRIPTS


def assess_indel_size(gene: str, c_notation: str) -> dict[str, Any]:
    """Measure an exactly described c. HGVS indel without sequence lookups.

    For delins alleles, the deleted and inserted lengths are retained
    separately. If they lie on opposite sides of the Appendix G boundary, this
    parser does not invent a single event-size convention. Symbolic, intronic,
    or uncertain descriptions remain unknown and require structural review.
    """
    transcript = TRANSCRIPTS.get((gene or "").upper(), "")
    return assess_indel_size_for_transcript(transcript, c_notation)
