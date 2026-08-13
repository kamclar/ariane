"""Build the protein-PS1 reference registry from official ENIGMA v1.2 data.

ST7 supplies P/LP missense reference candidates.  Table 9, complete ST2 and
the Table 9 SpliceAI value determine whether each candidate is currently
eligible, excluded, or requires review for protein-level PS1.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any, Dict


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from backend.modules.ps1_splice_evidence import (  # noqa: E402
    DEFINED_SOURCES,
    evaluate_defined_splice_sources,
)
from backend.modules.table9 import table9_lookup_ps3_bs3  # noqa: E402


DATA = PROJECT_ROOT / "backend" / "data"
ST7_PATH = DATA / "st7_reference_set.json"
TABLE9_PATH = DATA / "enigma_table9.json"
ST2_PATH = DATA / "enigma_st2_splice_evidence.json"
OUTPUT_PATH = DATA / "ps1_protein_reference_registry.json"
EXTENSIONS_PATH = DATA / "ps1_protein_reference_extensions.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record_checksum(record: Dict[str, Any]) -> str:
    basis = {key: value for key, value in record.items() if key != "approval_basis_checksum"}
    encoded = json.dumps(
        basis, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _is_normalized_missense(p_notation: str) -> bool:
    value = p_notation.replace("p.", "").replace("(", "").replace(")", "").strip()
    return bool(re.fullmatch(r"[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", value))


def _score(table9: Dict[str, Any]) -> float | None:
    raw = table9.get("spliceai_prediction")
    if raw in (None, ""):
        return None
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def _status(score: float | None, splice_status: str) -> tuple[str, str, str]:
    if splice_status == "conflicting":
        return (
            "review_required",
            "missense_mechanism_unresolved",
            "Defined ENIGMA sources contain conflicting splice evidence.",
        )
    if splice_status == "abnormal":
        return (
            "excluded",
            "missense_with_splice_effect",
            "Defined ENIGMA sources record a predicted or confirmed splice effect.",
        )
    if score is None:
        return (
            "review_required",
            "missense_mechanism_unresolved",
            "The reference-variant SpliceAI value is unavailable.",
        )
    if score > 0.1:
        return (
            "excluded",
            "missense_with_predicted_splice_effect",
            f"Reference-variant SpliceAI {score:.3f} is greater than 0.1.",
        )
    if splice_status not in {"normal", "none_identified"}:
        return (
            "review_required",
            "missense_mechanism_unresolved",
            "The defined RNA/splice source review is incomplete.",
        )
    return (
        "eligible",
        "missense_no_splice_effect",
        "ST7 P/LP missense reference meets the recorded ENIGMA protein-PS1 splice conditions.",
    )


def build() -> Dict[str, Any]:
    st7 = json.loads(ST7_PATH.read_text(encoding="utf-8"))
    source_checksums = {
        "st7_sha256": _sha256(ST7_PATH),
        "table9_sha256": _sha256(TABLE9_PATH),
        "st2_sha256": _sha256(ST2_PATH),
        "curated_extensions_sha256": _sha256(EXTENSIONS_PATH),
    }
    records = []
    for source in st7["variants"]:
        if source.get("iarc_class") not in {4, 5}:
            continue
        p_notation = str(source.get("p_notation") or "")
        if not _is_normalized_missense(p_notation):
            continue

        gene = source["gene"]
        c_notation = source["c_notation"]
        table9 = table9_lookup_ps3_bs3(gene, c_notation)
        splice = evaluate_defined_splice_sources(gene, c_notation, table9)
        spliceai_score = _score(table9)
        status, protein_branch, status_reason = _status(
            spliceai_score, splice["status"]
        )
        classification = (
            "Pathogenic" if source["iarc_class"] == 5 else "Likely Pathogenic"
        )
        record = {
            "reference_id": f"ENIGMA_ST7_V1_2|{gene}|{c_notation}",
            "gene": gene,
            "transcript": "NM_007294.4" if gene == "BRCA1" else "NM_000059.4",
            "c_notation": c_notation,
            "p_notation": p_notation,
            "classification": classification,
            "iarc_class": source["iarc_class"],
            "classification_verification": "enigma_st7_v1_2_reference_set",
            "classification_source": source.get("source") or "ENIGMA ST7 v1.2",
            "candidate_source": "ENIGMA Supplementary Table 7 v1.2",
            "status": status,
            "status_reason": status_reason,
            "protein_branch": protein_branch,
            "protein_mechanism_evidence": {
                "basis": (
                    "enigma_table9_ps3_functional_evidence"
                    if str(table9.get("code") or "").upper() == "PS3"
                    else "pathogenic_missense_with_no_predicted_or_confirmed_splice_effect"
                ),
                "table9_code": table9.get("code"),
                "table9_strength": table9.get("strength"),
                "table9_summary": table9.get("text"),
            },
            "reference_splice_evidence": {
                "spliceai_score": spliceai_score,
                "threshold": 0.1,
                "confirmed_status": splice["status"],
                "sources_checked": list(DEFINED_SOURCES),
                "checked_at": date.today().isoformat(),
                "source_details": splice,
                "provenance": {
                    "provider": "official_enigma_v1_2_snapshots",
                    "input_variant": f"{gene}:{c_notation}",
                    **source_checksums,
                },
            },
            "classification_ps1_dependency": {
                "used": "unknown",
                "reference_ids": [],
            },
        }
        record["approval_basis_checksum"] = _record_checksum(record)
        records.append(record)

    extensions = json.loads(EXTENSIONS_PATH.read_text(encoding="utf-8"))
    if extensions.get("schema_version") != 1 or extensions.get("status") != "active":
        raise RuntimeError("Protein PS1 curated extension file has unsupported metadata")
    extension_records = extensions.get("records")
    if not isinstance(extension_records, list):
        raise RuntimeError("Protein PS1 curated extension file has no records list")
    for source_record in extension_records:
        record = dict(source_record)
        record["approval_basis_checksum"] = _record_checksum(record)
        records.append(record)

    record_keys = [(item.get("gene"), item.get("c_notation")) for item in records]
    if len(record_keys) != len(set(record_keys)):
        raise RuntimeError("Protein PS1 sources contain a duplicate gene/c. reference")

    records.sort(key=lambda item: (item["gene"], item["c_notation"]))
    counts = Counter(record["status"] for record in records)
    return {
        "schema_version": 2,
        "registry_version": date.today().isoformat() + ".1",
        "status": "active",
        "description": (
            "ENIGMA ST7 v1.2 P/LP missense references for protein-level PS1, "
            "with explicit eligibility derived from the recorded ENIGMA splice checks."
        ),
        "rule_source": {
            "name": "ClinGen ENIGMA BRCA1/2 VCEP PS1 specification",
            "version": "1.2.0",
            "url": "https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2.0",
        },
        "candidate_source": {
            "name": "ClinGen ENIGMA BRCA1/2 VCEP Supplementary Table 7",
            "version": "1.2.0",
            "usage": "enigma_reference_set_with_per_record_ps1_eligibility",
        },
        "classification_policy": (
            "P/LP classifications in official ENIGMA ST7 v1.2 are accepted as the "
            "reference classification basis. Protein-PS1 eligibility additionally "
            "requires a missense consequence, SpliceAI <= 0.1 and no damaging "
            "splice evidence in the defined sources."
        ),
        "reference_source_policy": {
            "accepted_classification_bases": [
                {
                    "id": "enigma_st7_v1_2_reference_set",
                    "source": "ENIGMA Supplementary Table 7 v1.2",
                    "use": "Official P/LP reference classification; each missense record still requires the protein-PS1 splice checks."
                },
                {
                    "id": "external_vcep_assertion",
                    "source": "Versioned ENIGMA/ClinGen VCEP assertion, normally ClinGen Evidence Repository or a ClinVar expert-panel assertion",
                    "use": "May enter through the curated extension file after identity, transcript, mechanism, splice and provenance validation."
                },
                {
                    "id": "locally_recurated_under_enigma_vcep",
                    "source": "Documented local reclassification under a named ENIGMA VCEP specification version",
                    "use": "May enter through the curated extension file; it must not be labelled as an official expert-panel assertion."
                }
            ],
            "discovery_only_not_sufficient_for_eligibility": [
                "ClinVar record without an ENIGMA/ClinGen expert-panel assertion",
                "CANVarUK",
                "BRCA Exchange",
                "individual publication without a complete VCEP classification record",
                "computational prediction alone"
            ],
            "supporting_and_exclusion_sources": [
                "ENIGMA Specifications Table 9 v1.2",
                "ENIGMA Supplementary Table 2 v1.2",
                "versioned SpliceAI result with reference genome, transcript and model provenance",
                "canonical RefSeq transcript and normalized protein consequence"
            ]
        },
        "defined_splice_sources": list(DEFINED_SOURCES),
        "source_checksums": source_checksums,
        "reference_count": len(records),
        "status_counts": dict(sorted(counts.items())),
        "references": records,
    }


def main() -> None:
    data = build()
    OUTPUT_PATH.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Wrote {data['reference_count']} references to {OUTPUT_PATH}: "
        f"{data['status_counts']}"
    )


if __name__ == "__main__":
    main()
