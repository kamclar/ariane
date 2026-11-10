"""Source precedence for exact PVS1 RNA evidence."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


_EVIDENCE_ITEM_IDS = {
    "erepo": "erepo_pvs1_rna",
    "table4": "table4_pvs1_rna",
    "st2": "st2_pvs1_rna",
}


@dataclass(frozen=True)
class Pvs1RnaSourceSelection:
    """The first available source in the approved precedence order."""

    source: str | None
    result: Mapping[str, Any] | None

    @property
    def evidence_item_id(self) -> str | None:
        return _EVIDENCE_ITEM_IDS.get(self.source or "")

    @property
    def replaces_prediction(self) -> bool:
        if self.result is None:
            return False
        status = str(self.result.get("status") or "")
        return (
            (self.source == "erepo" and status == "eligible")
            or (
                self.source == "table4"
                and status in {"eligible", "not_applicable_rna"}
            )
            or (self.source == "st2" and status == "eligible")
        )


def select_pvs1_rna_source(
    erepo_result: Mapping[str, Any] | None,
    table4_result: Mapping[str, Any] | None,
    st2_result: Mapping[str, Any] | None,
) -> Pvs1RnaSourceSelection:
    """Select ERepo, Table 4, or ST2 without reinterpreting the evidence."""
    if erepo_result and erepo_result.get("status") != "not_in_registry":
        return Pvs1RnaSourceSelection("erepo", erepo_result)
    if table4_result and table4_result.get("status") not in {
        "not_in_source",
        "not_applicable",
    }:
        return Pvs1RnaSourceSelection("table4", table4_result)
    if st2_result and st2_result.get("record") is not None:
        return Pvs1RnaSourceSelection("st2", st2_result)
    return Pvs1RnaSourceSelection(None, None)
