"""Versioned BRCA1/2 pathogenic-founder exception lookup for BA1/BS1.

ENIGMA VCEP v1.2 prohibits BA1 and BS1 for well-established pathogenic
founder variants, but does not publish a machine-readable exhaustive list.
ARIANE therefore uses a small, provenance-bearing policy snapshot and never
falls back to an inferred or free-text runtime classification.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Dict

from backend.data_health import clear_issue, register_issue
from backend.gene_policy import active_genes


FOUNDER_VARIANT_SNAPSHOT = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "brca_pathogenic_founder_variants.json"
)

FOUNDER_VARIANT_INDEX: Dict[str, Dict[str, Any]] = {}
FOUNDER_VARIANT_METADATA: Dict[str, Any] = {}
FOUNDER_VARIANT_STATUS = "not_loaded"


def _records_sha256(records: list[dict[str, Any]]) -> str:
    payload = json.dumps(
        records,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _normalise_lookup_c(c_notation: str) -> str:
    value = re.sub(r"\s+", "", str(c_notation or ""))
    # Accepted sequence-explicit del/dup forms are equivalent only after the
    # upstream HGVS/reference validator has checked them. This lookup receives
    # that validated notation and merely makes the policy index alias-tolerant.
    value = re.sub(r"(del|dup)[ACGT]+$", r"\1", value, flags=re.IGNORECASE)
    if value.startswith("C."):
        value = "c." + value[2:]
    return value


def load_founder_variant_snapshot(path: Path | None = None) -> None:
    global FOUNDER_VARIANT_INDEX, FOUNDER_VARIANT_METADATA, FOUNDER_VARIANT_STATUS

    selected = Path(path) if path is not None else FOUNDER_VARIANT_SNAPSHOT
    try:
        payload = json.loads(selected.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        FOUNDER_VARIANT_INDEX = {}
        FOUNDER_VARIANT_METADATA = {}
        FOUNDER_VARIANT_STATUS = "unavailable"
        register_issue(
            "BRCA pathogenic founder variant snapshot",
            f"could not load {selected}: {type(exc).__name__}: {exc}",
        )
        return

    records = payload.get("variants")
    metadata = payload.get("metadata") or {}
    error = None
    if not isinstance(records, list) or not records:
        error = "variant records are missing or empty"
    elif metadata.get("enigma_version") != "1.2":
        error = "ENIGMA policy version is not 1.2"
    elif metadata.get("records_sha256") != _records_sha256(records):
        error = "records checksum mismatch"
    elif any(
        record.get("gene") not in set(active_genes())
        or not record.get("canonical_c_notation")
        or not record.get("source_assertions")
        for record in records
    ):
        error = "one or more records lack required identity or provenance"

    if error:
        FOUNDER_VARIANT_INDEX = {}
        FOUNDER_VARIANT_METADATA = metadata
        FOUNDER_VARIANT_STATUS = "invalid"
        register_issue("BRCA pathogenic founder variant snapshot", error)
        return

    index: Dict[str, Dict[str, Any]] = {}
    for record in records:
        gene = record["gene"].upper()
        aliases = {
            record["canonical_c_notation"],
            *(record.get("input_aliases") or []),
        }
        for alias in aliases:
            index[f"{gene}:{_normalise_lookup_c(alias)}"] = record

    FOUNDER_VARIANT_INDEX = index
    FOUNDER_VARIANT_METADATA = metadata
    FOUNDER_VARIANT_STATUS = "ok"
    clear_issue("BRCA pathogenic founder variant snapshot")


def lookup_pathogenic_founder_variant(gene: str, c_notation: str) -> Dict[str, Any]:
    """Return a fail-closed founder-exception decision for BA1/BS1."""
    if FOUNDER_VARIANT_STATUS != "ok":
        return {
            "status": "unavailable",
            "is_pathogenic_founder": None,
            "reason": (
                "the versioned BRCA pathogenic founder variant snapshot is "
                f"{FOUNDER_VARIANT_STATUS}"
            ),
        }
    if not gene or not c_notation:
        return {
            "status": "unavailable",
            "is_pathogenic_founder": None,
            "reason": "gene and canonical c. notation are required for the founder exception check",
        }

    key = f"{gene.upper()}:{_normalise_lookup_c(c_notation)}"
    record = FOUNDER_VARIANT_INDEX.get(key)
    if record is None:
        return {
            "status": "not_found",
            "is_pathogenic_founder": False,
            "reason": "variant is not present in the approved pathogenic-founder snapshot",
            "snapshot_version": FOUNDER_VARIANT_METADATA.get("snapshot_version"),
        }
    return {
        "status": "found",
        "is_pathogenic_founder": True,
        "reason": record.get("founder_context") or "well-established pathogenic founder variant",
        "record": record,
        "snapshot_version": FOUNDER_VARIANT_METADATA.get("snapshot_version"),
    }


load_founder_variant_snapshot()
