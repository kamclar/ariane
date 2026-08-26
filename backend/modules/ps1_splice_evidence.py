"""Defined-source check for confirmed RNA/splice evidence used by protein PS1."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


ST2_EVIDENCE_PATH = (
    Path(__file__).resolve().parents[1] / "data" / "enigma_st2_splice_evidence.json"
)
DEFINED_SOURCES = [
    "ENIGMA Specifications Table 9 v1.2",
    "ENIGMA Supplementary Table 2 v1.2",
]

_ST2_BY_VARIANT: Optional[Dict[Tuple[str, str], Dict[str, Any]]] = None
_ST2_PAYLOAD: Optional[Dict[str, Any]] = None


def _load_st2_payload() -> Dict[str, Any]:
    global _ST2_PAYLOAD
    if _ST2_PAYLOAD is not None:
        return _ST2_PAYLOAD
    if not ST2_EVIDENCE_PATH.is_file():
        raise RuntimeError(
            f"Required complete ENIGMA ST2 splice evidence snapshot is missing: {ST2_EVIDENCE_PATH}"
        )
    data = json.loads(ST2_EVIDENCE_PATH.read_text(encoding="utf-8"))
    if (
        data.get("schema_version") != 1
        or data.get("source_columns") != 11
        or data.get("total_variants") != 220
        or len(data.get("variants", [])) != 220
    ):
        raise RuntimeError("ENIGMA ST2 splice evidence snapshot is incomplete")
    _ST2_PAYLOAD = data
    return data


def _load_st2_evidence() -> Dict[Tuple[str, str], Dict[str, Any]]:
    global _ST2_BY_VARIANT
    if _ST2_BY_VARIANT is not None:
        return _ST2_BY_VARIANT
    data = _load_st2_payload()
    _ST2_BY_VARIANT = {
        (record["gene"], record["c_notation"]): record
        for record in data["variants"]
    }
    if len(_ST2_BY_VARIANT) != 220:
        raise RuntimeError("ENIGMA ST2 splice evidence snapshot has duplicate variants")
    return _ST2_BY_VARIANT


def get_st2_splice_record(gene: str, c_notation: str) -> Optional[Dict[str, Any]]:
    """Return the exact official ST2 row for a normalized BRCA variant."""
    return _load_st2_evidence().get((gene, c_notation))


def list_splice_ps1_candidate_discovery(gene: Optional[str] = None) -> Dict[str, Any]:
    """Return factual P/LP splice candidates derived directly from official ST2.

    These records support candidate discovery only. They do not establish that
    a reference is eligible for PS1(splicing), that a VUA has the same event,
    or which Appendix J/Table 17 strength applies.
    """
    payload = _load_st2_payload()
    requested_gene = gene.upper() if gene else None
    candidates = []
    for record in payload["variants"]:
        record_gene = str(record.get("gene") or "")
        result = str(record.get("result") or "").strip()
        classification_numeric = record.get("final_multifactorial_class")
        if requested_gene and record_gene != requested_gene:
            continue
        if classification_numeric not in {4, 5}:
            continue
        if not result or result.lower() == "no aberration":
            continue
        classification = "Pathogenic" if classification_numeric == 5 else "Likely Pathogenic"
        source_row = int(record["source_row"])
        c_notation = str(record.get("c_notation") or "")
        candidates.append({
            "key": f"{record_gene}|{c_notation}|{source_row}",
            "gene": record_gene,
            "reference_variant": c_notation,
            "p_notation": str(record.get("p_notation") or ""),
            "classification": classification,
            "classification_numeric": classification_numeric,
            "classification_basis": "ENIGMA ST2 final multifactorial class",
            "reference_splice_event": result,
            "assay_result_category": str(record.get("splicing_assay_result_category") or ""),
            "assay_context": str(record.get("variant_assay_summary") or ""),
            "included_in_analysis": str(record.get("included_in_analysis") or ""),
            "prior_probability": record.get("prior_probability"),
            "source_row": source_row,
            "source_label": f"ENIGMA Supplementary Table 2 v1.2, source row {source_row}",
            "source_url": str(payload.get("source_url") or ""),
            "source_file_sha256": str(payload.get("source_file_sha256") or ""),
            "eligibility_status": "candidate_discovery_only",
            "eligibility_note": (
                "The ST2 record does not by itself confirm PS1(splicing) eligibility, "
                "same-event matching, prediction-strength comparison, or Appendix J/Table 17 strength."
            ),
        })
    return {
        "status": "candidate_discovery_only",
        "source": "ENIGMA Supplementary Table 2 v1.2",
        "source_url": str(payload.get("source_url") or ""),
        "source_file_sha256": str(payload.get("source_file_sha256") or ""),
        "candidate_count": len(candidates),
        "candidates": candidates,
    }


def evaluate_defined_splice_sources(
    gene: str,
    c_notation: str,
    table9_result: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Return a qualified, auditable known-splice-evidence status.

    ``none_identified`` means only that these two named, versioned ENIGMA
    sources were checked for the exact variant and no abnormal result was
    found. It is not an absolute assertion that no splice evidence exists.
    """
    st2_record = get_st2_splice_record(gene, c_notation)
    table9_result = table9_result or {}
    table9_reviewed = bool(table9_result.get("reviewed"))
    table9_flag = str(
        table9_result.get("predicted_or_observed_splicing") or ""
    ).strip()
    table9_published = str(
        table9_result.get("splice_result_published") or ""
    ).strip()
    st2_result = str((st2_record or {}).get("result") or "").strip()

    # Table 9 has a structured conclusion.  In particular, all rows labelled
    # ``N, no aberration`` also carry a publication result such as
    # ``no aberration (PMID: ...)``.  A non-empty publication field therefore
    # cannot itself be interpreted as abnormal splicing.
    table9_reports_splice = table9_flag.upper().startswith("Y")
    table9_reports_no_splice = table9_flag.upper().startswith("N")
    st2_reports_no_splice = st2_result.lower() == "no aberration"
    st2_reports_splice = bool(st2_record) and bool(st2_result) and not st2_reports_no_splice

    if table9_reports_no_splice and st2_reports_splice:
        status = "conflicting"
        reason = "ENIGMA Table 9 and Supplementary Table 2 report conflicting splice conclusions"
    elif table9_reports_splice or st2_reports_splice:
        status = "abnormal"
        locations = []
        if table9_reports_splice:
            locations.append("ENIGMA Specifications Table 9")
        if st2_reports_splice:
            locations.append("ENIGMA Supplementary Table 2")
        reason = "Predicted or confirmed abnormal splice evidence is recorded in " + " and ".join(locations)
    elif table9_reports_no_splice or st2_reports_no_splice:
        status = "normal"
        reason = "A defined ENIGMA source explicitly reports no splice aberration"
    else:
        status = "none_identified"
        reason = (
            "No confirmed abnormal splice effect was identified for the exact variant "
            "in the complete ENIGMA Table 9 and Supplementary Table 2 snapshots"
        )

    return {
        "status": status,
        "sources_checked": list(DEFINED_SOURCES),
        "reason": reason,
        "table9_reviewed": table9_reviewed,
        "table9_splice_flag": table9_flag,
        "table9_published_splice_result": table9_published,
        "supplementary_table2_match": bool(st2_record),
        "supplementary_table2_result": st2_result,
        "absence_semantics": (
            "No exact record means only that no evidence was identified in these "
            "defined source versions, not that no splice evidence exists."
        ),
    }


def reset_splice_source_cache_for_tests() -> None:
    global _ST2_BY_VARIANT, _ST2_PAYLOAD
    _ST2_BY_VARIANT = None
    _ST2_PAYLOAD = None
