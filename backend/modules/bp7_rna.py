"""ENIGMA BP7 Strong (RNA) variant-context eligibility.

The RNA assay record and the variant-context stipulations are deliberately
evaluated separately.  A reviewer can describe a valid mRNA-only assay while
the assessed variant is still ineligible for BP7 Strong (RNA).
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence
import re

from backend.config import FUNCTIONAL_DOMAINS
from backend.modules.variant_type import infer_variant_type


_INFRAME_TYPES = {
    "inframe_deletion",
    "inframe_insertion",
    "inframe_delins",
    "delins",
}


def _protein_interval(p_notation: str) -> tuple[int, int] | None:
    """Return the complete affected amino-acid interval when it is explicit."""
    positions = [
        int(value)
        for value in re.findall(r"[A-Z][a-z]{2}(\d+)", p_notation or "")
    ]
    if not positions:
        return None
    return min(positions), max(positions)


def _overlapping_domains(
    gene: str,
    protein_interval: tuple[int, int],
) -> tuple[str, ...]:
    start, end = protein_interval
    return tuple(
        name
        for name, (domain_start, domain_end) in FUNCTIONAL_DOMAINS.get(gene, {}).items()
        if start <= domain_end and end >= domain_start
    )


def _has_table9_bs3(base_criteria: Sequence[Mapping[str, Any]]) -> bool:
    """Accept only an applied BS3 carrying ENIGMA Table 9 provenance."""
    for criterion in base_criteria:
        if (
            str(criterion.get("name") or "").upper() != "BS3"
            or criterion.get("applies") is not True
            or not criterion.get("strength")
        ):
            continue
        path = criterion.get("decision_path")
        if not isinstance(path, Mapping):
            continue
        sources = path.get("sources")
        if not isinstance(sources, list):
            continue
        if any(
            isinstance(source, Mapping)
            and source.get("source_id") == "enigma-v1.2-table9"
            for source in sources
        ):
            return True
    return False


def evaluate_bp7_rna_variant_context(
    variant_context: Mapping[str, Any] | None,
    base_criteria: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Apply the BP7 Strong (RNA) location and BS3 stipulations fail-closed."""
    if not isinstance(variant_context, Mapping):
        return {
            "eligible": False,
            "status": "unavailable",
            "reason": "Variant context is required for BP7 Strong (RNA).",
        }

    gene = str(variant_context.get("gene") or "").strip().upper()
    c_notation = str(variant_context.get("c_notation") or "").strip()
    p_notation = str(variant_context.get("p_notation") or "").strip()
    if gene not in FUNCTIONAL_DOMAINS or not c_notation or not p_notation:
        return {
            "eligible": False,
            "status": "unavailable",
            "reason": "Gene, c. notation and p. notation are required for BP7 Strong (RNA).",
        }

    variant_type = infer_variant_type(c_notation, p_notation).lower()
    details: dict[str, Any] = {
        "eligible": False,
        "status": "not_eligible",
        "gene": gene,
        "variant_type": variant_type,
        "p_notation": p_notation,
        "requires_bs3": False,
        "bs3_met": False,
        "functional_domains": [],
    }

    if variant_type in {"intronic", "synonymous", "silent"}:
        details.update(
            eligible=True,
            status="eligible",
            reason=(
                f"ENIGMA permits BP7 Strong (RNA) for {variant_type} variants "
                "when the mRNA-only assay requirements are met."
            ),
        )
        return details

    if variant_type not in {"missense", *_INFRAME_TYPES}:
        details["reason"] = (
            f"ENIGMA BP7 Strong (RNA) is not applicable to {variant_type} variants."
        )
        return details

    protein_interval = _protein_interval(p_notation)
    if protein_interval is None:
        details.update(
            status="unavailable",
            reason=(
                "The complete affected amino-acid position could not be determined; "
                "BP7 Strong (RNA) was not applied."
            ),
        )
        return details

    domains = _overlapping_domains(gene, protein_interval)
    details["protein_interval"] = list(protein_interval)
    details["functional_domains"] = list(domains)
    if not domains:
        details.update(
            eligible=True,
            status="eligible",
            reason=(
                "The missense/in-frame variant is outside the ENIGMA clinically "
                "important functional domains."
            ),
        )
        return details

    if variant_type != "missense":
        details["reason"] = (
            "ENIGMA permits BP7 Strong (RNA) for in-frame variants outside a "
            "clinically important functional domain; this variant overlaps "
            f"{', '.join(domains)}."
        )
        return details

    details["requires_bs3"] = True
    details["bs3_met"] = _has_table9_bs3(base_criteria)
    if details["bs3_met"]:
        details.update(
            eligible=True,
            status="eligible",
            reason=(
                "The missense variant is inside an ENIGMA functional domain and "
                "BS3 is met by calibrated ENIGMA Table 9 evidence."
            ),
        )
    else:
        details["reason"] = (
            "The missense variant is inside an ENIGMA functional domain "
            f"({', '.join(domains)}), but BS3 is not met by calibrated ENIGMA "
            "Table 9 evidence. BP7 Strong (RNA) cannot be applied."
        )
    return details
