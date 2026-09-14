"""Public automatic-classification scope for normalized variant types.

Recognizing an HGVS consequence is not the same as having a complete ENIGMA
Module 1 path for it.  This policy is checked before a result may be published
or stored in the completed-result cache.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


SUPPORTED_AUTOMATIC_VARIANT_TYPES = (
    "nonsense",
    "frameshift",
    "splice_site",
    "initiation_codon",
    "exon_deletion",
    "exon_duplication",
    "missense",
    "synonymous",
    "intronic",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_delins",
)

# These are valid, recognizable consequence families, but the current ARIANE
# graph has no complete ENIGMA v1.2 automatic path for them.
OUT_OF_SCOPE_VARIANT_TYPES = (
    "5utr",
    "3utr",
    "stop_lost",
)

# These tokens mean that normalization did not establish enough consequence
# information to choose a PTC, Figure 1A, or structural path.
UNRESOLVED_VARIANT_TYPES = (
    "unknown",
    "delins",
    "deletion",
    "insertion",
    "duplication",
)

_SUPPORTED_TYPES = frozenset((*SUPPORTED_AUTOMATIC_VARIANT_TYPES, "silent"))
_OUT_OF_SCOPE_TYPES = frozenset(OUT_OF_SCOPE_VARIANT_TYPES)
_UNRESOLVED_TYPES = frozenset(UNRESOLVED_VARIANT_TYPES)


@dataclass(frozen=True)
class ClassificationScopeDecision:
    variant_type: str
    status: Literal["supported", "out_of_scope", "unresolved"]
    error_code: str = ""
    reason: str = ""

    @property
    def supported(self) -> bool:
        return self.status == "supported"


def automatic_classification_scope(variant_type: str) -> ClassificationScopeDecision:
    """Return whether a normalized type has a complete public Module 1 path."""
    normalized = (variant_type or "").strip().lower()
    if normalized in _SUPPORTED_TYPES:
        return ClassificationScopeDecision(normalized, "supported")

    if normalized in _UNRESOLVED_TYPES:
        return ClassificationScopeDecision(
            normalized or "unknown",
            "unresolved",
            "protein_consequence_unresolved",
            (
                "The normalized protein consequence does not identify a complete "
                "ENIGMA PTC, Figure 1A, or structural-variant path. No "
                "classification was returned."
            ),
        )

    display_type = {
        "5utr": "5' UTR",
        "3utr": "3' UTR",
        "stop_lost": "stop-loss",
    }.get(normalized, normalized or "unknown")
    return ClassificationScopeDecision(
        normalized or "unknown",
        "out_of_scope",
        "unsupported_variant_type",
        (
            f"{display_type} variants do not have a complete automatic path in "
            "the current ARIANE implementation of ENIGMA BRCA1/2 VCEP v1.2. "
            "No classification was returned."
        ),
    )


def automatic_classification_is_supported(variant_type: str) -> bool:
    """Return whether a public automatic result may be published and cached."""
    return automatic_classification_scope(variant_type).supported
