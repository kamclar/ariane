"""Keep the reviewed PVS1 RNA registry aligned with current ERepo assertions."""

from __future__ import annotations

import json
from pathlib import Path
import re


DATA_DIR = Path(__file__).resolve().parents[1] / "backend" / "data"
RNA_CODE_RE = re.compile(
    r"PVS1(?:_(?:Strong|Moderate|Supporting))? \(RNA\)"
)


def _load(name: str) -> dict:
    return json.loads((DATA_DIR / name).read_text(encoding="utf-8"))


def test_reviewed_rna_registry_covers_all_current_explicit_v12_assertions():
    """A source refresh must not silently omit an explicit final RNA weight.

    This is deliberately a release-time consistency guard, not runtime text
    interpretation.  Records still enter the checksum-bound registry only
    after review; the test makes a newly exposed current v1.2 assertion fail
    CI until that review is completed.
    """
    erepo = _load("enigma_erepo_vcep_registry.json")
    approved = _load("enigma_erepo_pvs1_rna_registry.json")

    current_explicit: dict[str, tuple[str, str, str]] = {}
    for record in erepo["records"]:
        if record.get("source_status") != "current_vcep_assertion":
            continue
        if record.get("assertion_method_version") != "1.2.0":
            continue
        source_codes = set(RNA_CODE_RE.findall(record.get("interpretation_summary", "")))
        if not source_codes:
            continue
        assert len(source_codes) == 1, (
            "Current ERepo assertion has more than one final PVS1 RNA weight: "
            f"{record['uuid']} {sorted(source_codes)}"
        )
        source_code = source_codes.pop()
        criterion_code = source_code.removesuffix(" (RNA)")
        assert criterion_code in record.get("criteria_met", []), (
            "The explicit PVS1 RNA weight disagrees with criteria_met for "
            f"{record['uuid']}"
        )
        current_explicit[record["uuid"]] = (
            record["gene"],
            record["c_notation"],
            source_code,
        )

    approved_by_uuid = {
        record["assertion_uuid"]: (
            record["gene"],
            record["c_notation"],
            record["source_code"],
        )
        for record in approved["records"]
    }

    assert approved_by_uuid == current_explicit

