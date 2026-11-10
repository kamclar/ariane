"""Resolved PVS1 RNA evidence used across the classification DAG."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


_EVIDENCE_ITEM_IDS = {
    "erepo": "erepo_pvs1_rna",
    "table4": "table4_pvs1_rna",
    "st2": "st2_pvs1_rna",
}


@dataclass(frozen=True)
class Pvs1RnaSelection:
    """One selected source and its final PVS1 RNA decision."""

    source: str | None
    decision: Mapping[str, Any]
    source_result: Mapping[str, Any] | None

    @property
    def evidence_item_id(self) -> str | None:
        return _EVIDENCE_ITEM_IDS.get(self.source or "")

    @property
    def replaces_prediction(self) -> bool:
        return bool(self.decision.get("applies")) or (
            self.decision.get("application_status")
            == "not_applicable_from_table4"
        )

    def matches_source_result(
        self,
        result: Mapping[str, Any] | None,
    ) -> bool:
        """Return whether the resolved decision still matches its input."""
        return self.source_result == result
