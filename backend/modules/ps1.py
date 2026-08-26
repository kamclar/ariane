"""ENIGMA v1.2 protein-level PS1 evaluation.

ST7 is an official, trusted ENIGMA source for discovering candidate reference
variants.  It is not treated as an automatic PS1 allowlist.  Automatic scoring
uses only the separate, versioned registry whose records document the VCEP
classification and PS1-specific splice review required by ENIGMA.
"""

from __future__ import annotations

import json
import hashlib
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from backend.gene_policy import (
    active_genes,
    reference_transcript,
    spliceai_thresholds,
    vcep_specification,
)


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
ST7_PATH = DATA_DIR / "st7_reference_set.json"
PS1_REGISTRY_PATH = DATA_DIR / "ps1_protein_reference_registry.json"

PS1_POINTS = {"Strong": 4, "Moderate": 2}
APPROVED_CLASSIFICATION_VERIFICATIONS = {
    "external_vcep_assertion",
    "locally_recurated_under_enigma_vcep",
    "enigma_st7_v1_2_reference_set",
}
NO_DAMAGING_SPLICE_STATUSES = {"none_identified", "normal"}
_ST7_LOOKUP: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
_APPROVED_LOOKUP: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
_LOADED = False


def _extract_aa_change(p_notation: str) -> Optional[str]:
    """Return a canonical three-letter missense key, for example Arg170Gln."""
    if not p_notation:
        return None
    clean = (
        p_notation.replace("p.", "")
        .replace("(", "")
        .replace(")", "")
        .strip()
    )
    if any(token in clean.lower() for token in ("fs", "ter", "del", "ins", "dup", "=", "?")):
        return None
    return clean if re.fullmatch(r"[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", clean) else None


def select_vua_spliceai_for_ps1(
    service_or_cache_score: Optional[float],
    table9_result: Optional[Dict[str, Any]],
) -> tuple[Optional[float], str]:
    """Return the configured SpliceAI result used for protein PS1.

    Table 9 is used separately for curated PS3/BS3 and published splice
    evidence. Its recorded prediction does not replace the configured result.
    The second argument remains in the public helper signature for callers that
    already provide the reviewed Table 9 record.
    """
    return service_or_cache_score, "configured SpliceAI source"


def _candidate_from_st7(record: Dict[str, Any]) -> Dict[str, Any]:
    iarc_class = record.get("iarc_class")
    classification = "Pathogenic" if iarc_class == 5 else "Likely Pathogenic"
    return {
        "key": f"ST7|{record['gene']}|{record['c_notation']}",
        "reference_id": "",
        "gene": record["gene"],
        "transcript": reference_transcript(record["gene"]),
        "c_notation": record["c_notation"],
        "p_notation": record.get("p_notation") or "",
        "classification": classification,
        "iarc_class": iarc_class,
        "classification_basis": "enigma_multifactorial_likelihood_reference_set",
        "classification_source": record.get("source") or "",
        "reference_status": "review_required",
        "status_reason": (
            "Official ENIGMA ST7 candidate. A qualifying VCEP classification "
            "and completed PS1 splice review are not recorded in ST7."
        ),
        "source_dataset": "ENIGMA Supplementary Table 7 v1.2",
    }


def _validate_sha256(value: Any, label: str) -> None:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise RuntimeError(f"PS1 reference registry has invalid {label}")


def compute_approval_basis_checksum(record: Dict[str, Any]) -> str:
    """Checksum every approved-record field except the checksum itself."""
    basis = {
        key: value
        for key, value in record.items()
        if key != "approval_basis_checksum"
    }
    canonical = json.dumps(
        basis,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def validate_ps1_reference_registry(data: Dict[str, Any]) -> None:
    """Validate the curated automatic-scoring registry and known dependencies."""
    if data.get("schema_version") != 3 or data.get("status") != "active":
        raise RuntimeError("PS1 protein reference registry has unsupported metadata")
    if not str(data.get("registry_version") or "").strip():
        raise RuntimeError("PS1 protein reference registry has no registry_version")
    defined_sources = data.get("defined_splice_sources")
    if not isinstance(defined_sources, list) or not defined_sources:
        raise RuntimeError("PS1 protein reference registry has no defined splice sources")
    source_policy = data.get("reference_source_policy")
    accepted_bases = (
        source_policy.get("accepted_classification_bases")
        if isinstance(source_policy, dict)
        else None
    )
    if not isinstance(accepted_bases, list) or {
        item.get("id") for item in accepted_bases if isinstance(item, dict)
    } != APPROVED_CLASSIFICATION_VERIFICATIONS:
        raise RuntimeError("PS1 protein reference registry has invalid reference source policy")
    source_checksums = data.get("source_checksums")
    required_checksum_keys = {
        "st7_sha256",
        "table9_sha256",
        "st2_sha256",
        "curated_extensions_sha256",
    }
    if not isinstance(source_checksums, dict) or set(source_checksums) != required_checksum_keys:
        raise RuntimeError("PS1 protein reference registry has invalid source checksums")
    for label, value in source_checksums.items():
        _validate_sha256(value, label)
    records = data.get("references")
    if not isinstance(records, list):
        raise RuntimeError("PS1 protein reference registry has no references list")
    if data.get("reference_count") != len(records):
        raise RuntimeError("PS1 protein reference registry has an invalid reference_count")
    actual_counts: Dict[str, int] = {}
    for record in records:
        record_status = str(record.get("status") or "")
        actual_counts[record_status] = actual_counts.get(record_status, 0) + 1
    if data.get("status_counts") != dict(sorted(actual_counts.items())):
        raise RuntimeError("PS1 protein reference registry has invalid status_counts")

    by_id: Dict[str, Dict[str, Any]] = {}
    variant_keys = set()
    for index, record in enumerate(records):
        prefix = f"PS1 approved reference #{index}"
        reference_id = str(record.get("reference_id") or "")
        gene = record.get("gene")
        c_notation = str(record.get("c_notation") or "")
        p_notation = str(record.get("p_notation") or "")
        if not reference_id or reference_id in by_id:
            raise RuntimeError(f"{prefix} has missing or duplicate reference_id")
        if gene not in set(active_genes()) or not c_notation.startswith("c."):
            raise RuntimeError(f"{prefix} has invalid variant identity")
        if not _extract_aa_change(p_notation):
            raise RuntimeError(f"{prefix} is not a normalized missense variant")
        expected_transcript = reference_transcript(gene)
        if record.get("transcript") != expected_transcript:
            raise RuntimeError(f"{prefix} does not use the ENIGMA reference transcript")
        variant_key = (gene, c_notation)
        if variant_key in variant_keys:
            raise RuntimeError(f"PS1 registry has duplicate variant {gene}:{c_notation}")
        variant_keys.add(variant_key)

        status = record.get("status")
        if status not in {"eligible", "excluded", "review_required"}:
            raise RuntimeError(f"{prefix} has invalid eligibility status")
        if record.get("classification") not in {"Pathogenic", "Likely Pathogenic"}:
            raise RuntimeError(f"{prefix} lacks a P/LP classification")
        verification = record.get("classification_verification")
        if verification not in APPROVED_CLASSIFICATION_VERIFICATIONS:
            raise RuntimeError(f"{prefix} lacks qualifying VCEP classification verification")
        if not str(record.get("classification_source") or "").strip():
            raise RuntimeError(f"{prefix} lacks classification_source")
        if verification == "external_vcep_assertion":
            assertion = record.get("classification_assertion")
            if not isinstance(assertion, dict) or any(
                not str(assertion.get(field) or "").strip()
                for field in ("organization", "assertion_id", "ruleset_version", "accessed_at")
            ):
                raise RuntimeError(f"{prefix} lacks external VCEP assertion provenance")
        elif verification == "locally_recurated_under_enigma_vcep":
            local_review = record.get("local_reclassification")
            if not isinstance(local_review, dict) or any(
                not str(local_review.get(field) or "").strip()
                for field in ("reviewer_id", "reviewed_at", "ruleset_version", "evidence_record_id")
            ):
                raise RuntimeError(f"{prefix} lacks local reclassification provenance")
        if not str(record.get("status_reason") or "").strip():
            raise RuntimeError(f"{prefix} lacks status_reason")
        mechanism = record.get("protein_mechanism_evidence")
        if not isinstance(mechanism, dict) or mechanism.get("basis") not in {
            "enigma_table9_ps3_functional_evidence",
            "pathogenic_missense_with_no_predicted_or_confirmed_splice_effect",
            "curated_protein_mechanism_assessment",
        }:
            raise RuntimeError(f"{prefix} lacks auditable protein mechanism evidence")

        splice = record.get("reference_splice_evidence")
        if not isinstance(splice, dict):
            raise RuntimeError(f"{prefix} lacks reference_splice_evidence")
        threshold = splice.get("threshold")
        expected_threshold = spliceai_thresholds(record["gene"])["bp4"]
        if threshold != expected_threshold or splice.get("prediction_policy") != "runtime_required":
            raise RuntimeError(f"{prefix} has invalid runtime SpliceAI policy")
        if "spliceai_score" in splice:
            raise RuntimeError(f"{prefix} embeds a SpliceAI score in the registry")
        if not splice.get("sources_checked") or not splice.get("checked_at"):
            raise RuntimeError(f"{prefix} lacks an auditable confirmed-splice source review")
        if not set(defined_sources).issubset(set(splice["sources_checked"])):
            raise RuntimeError(f"{prefix} has not checked every registry-defined splice source")
        if not isinstance(splice.get("provenance"), dict) or not splice["provenance"]:
            raise RuntimeError(f"{prefix} lacks SpliceAI provenance")
        confirmed_status = splice.get("confirmed_status")
        if status == "eligible":
            if record.get("protein_branch") != "missense_runtime_spliceai_check_required":
                raise RuntimeError(f"{prefix} lacks protein-missense branch approval")
            if confirmed_status not in NO_DAMAGING_SPLICE_STATUSES:
                raise RuntimeError(f"{prefix} has unresolved confirmed splice evidence")
        elif status == "excluded":
            if confirmed_status != "abnormal":
                raise RuntimeError(f"{prefix} has no recorded reason for exclusion")
        _validate_sha256(record.get("approval_basis_checksum"), f"approval_basis_checksum for {reference_id}")
        if record["approval_basis_checksum"] != compute_approval_basis_checksum(record):
            raise RuntimeError(f"{prefix} approval_basis_checksum does not match its content")

        dependency = record.get("classification_ps1_dependency", {})
        used = dependency.get("used", "unknown")
        dependencies = dependency.get("reference_ids", [])
        if used not in {True, False, "unknown"} or not isinstance(dependencies, list):
            raise RuntimeError(f"{prefix} has invalid PS1 dependency metadata")
        if used is True and not dependencies:
            raise RuntimeError(f"{prefix} used PS1 but does not identify its reference")
        if reference_id in dependencies:
            raise RuntimeError(f"{prefix} directly depends on itself")
        by_id[reference_id] = record

    # Known PS1 dependencies must not contain direct or longer registry cycles.
    visiting = set()
    visited = set()

    def visit(reference_id: str) -> None:
        if reference_id in visiting:
            raise RuntimeError(f"PS1 reference registry contains a circular dependency at {reference_id}")
        if reference_id in visited:
            return
        visiting.add(reference_id)
        record = by_id[reference_id]
        dependency = record.get("classification_ps1_dependency", {})
        if dependency.get("used") is True:
            for child in dependency.get("reference_ids", []):
                if child in by_id:
                    visit(child)
        visiting.remove(reference_id)
        visited.add(reference_id)

    for reference_id in by_id:
        visit(reference_id)


def _load_references() -> None:
    global _LOADED, _ST7_LOOKUP, _APPROVED_LOOKUP
    if _LOADED:
        return
    if not ST7_PATH.is_file():
        raise RuntimeError(f"Required ENIGMA ST7 dataset is missing: {ST7_PATH}")
    if not PS1_REGISTRY_PATH.is_file():
        raise RuntimeError(f"Required PS1 protein reference registry is missing: {PS1_REGISTRY_PATH}")

    st7 = json.loads(ST7_PATH.read_text(encoding="utf-8"))
    registry = json.loads(PS1_REGISTRY_PATH.read_text(encoding="utf-8"))
    validate_ps1_reference_registry(registry)

    expected_st7 = {
        (record.get("gene"), record.get("c_notation")): record
        for record in st7.get("variants", [])
        if record.get("iarc_class") in {4, 5}
        and _extract_aa_change(record.get("p_notation") or "")
    }
    registry_records = registry.get("references", [])
    registry_by_variant = {
        (record.get("gene"), record.get("c_notation")): record
        for record in registry_records
    }
    st7_registry = {
        key: record
        for key, record in registry_by_variant.items()
        if record.get("classification_verification") == "enigma_st7_v1_2_reference_set"
    }
    if set(st7_registry) != set(expected_st7):
        raise RuntimeError(
            "PS1 protein reference registry does not contain the complete ST7 P/LP missense set"
        )
    for key, source in expected_st7.items():
        record = st7_registry[key]
        expected_class = "Pathogenic" if source["iarc_class"] == 5 else "Likely Pathogenic"
        if (
            record.get("p_notation") != source.get("p_notation")
            or record.get("classification") != expected_class
            or record.get("classification_source") != (source.get("source") or "ENIGMA ST7 v1.2")
        ):
            raise RuntimeError(f"PS1 registry record does not match ST7 for {key[0]}:{key[1]}")

    st7_lookup: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    for record in registry_records:
        aa_key = _extract_aa_change(record["p_notation"])
        st7_lookup.setdefault(record["gene"], {}).setdefault(aa_key, []).append(
            _public_registry_candidate(record)
        )

    approved_lookup: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    for record in registry_records:
        if record.get("status") != "eligible":
            continue
        aa_key = _extract_aa_change(record["p_notation"])
        approved_lookup.setdefault(record["gene"], {}).setdefault(aa_key, []).append(record)

    _ST7_LOOKUP = st7_lookup
    _APPROVED_LOOKUP = approved_lookup
    _LOADED = True


def reset_ps1_reference_cache_for_tests() -> None:
    global _LOADED, _ST7_LOOKUP, _APPROVED_LOOKUP
    _LOADED = False
    _ST7_LOOKUP = {}
    _APPROVED_LOOKUP = {}


def discover_ps1_reference_variants(
    gene: str,
    c_notation: str,
    p_notation: str,
    variant_type: str,
) -> List[str]:
    """Return candidate reference c. variants that require runtime SpliceAI.

    Discovery is local and deterministic. It does not apply PS1 and it does not
    call an external service. The caller can therefore schedule reference and
    assessed-variant SpliceAI lookups together.
    """
    if (variant_type or "").lower() != "missense":
        return []
    aa_key = _extract_aa_change(p_notation)
    if not aa_key:
        return []
    _load_references()
    return sorted(
        {
            str(record["c_notation"])
            for record in _APPROVED_LOOKUP.get(gene, {}).get(aa_key, [])
            if record.get("c_notation") != c_notation
        }
    )


def _deduplicate_candidates(records: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    deduplicated: Dict[tuple, Dict[str, Any]] = {}
    for record in records:
        key = (record.get("gene"), record.get("c_notation"))
        current = deduplicated.get(key)
        if current is None or record.get("reference_status") == "approved":
            deduplicated[key] = record
    return list(deduplicated.values())


def _public_registry_candidate(record: Dict[str, Any]) -> Dict[str, Any]:
    registry_status = record.get("status", "review_required")
    return {
        "key": f"REGISTRY|{record['reference_id']}",
        "reference_id": record["reference_id"],
        "gene": record["gene"],
        "transcript": record["transcript"],
        "c_notation": record["c_notation"],
        "p_notation": record["p_notation"],
        "classification": record["classification"],
        "iarc_class": 5 if record["classification"] == "Pathogenic" else 4,
        "classification_basis": record["classification_verification"],
        "classification_source": record["classification_source"],
        "reference_status": "approved" if registry_status == "eligible" else registry_status,
        "status_reason": record.get("status_reason", ""),
        "source_dataset": record.get("candidate_source", "Curated PS1 protein reference registry"),
    }


def _public_approved_candidate(record: Dict[str, Any]) -> Dict[str, Any]:
    return _public_registry_candidate(record)


def evaluate_ps1(
    gene: str,
    c_notation: str,
    p_notation: str,
    variant_type: str,
    spliceai_score: Optional[float] = None,
    vua_spliceai_source: str = "configured SpliceAI source",
    vua_splice_evidence_status: str = "not_assessed",
    vua_splice_sources_checked: Optional[List[str]] = None,
    reference_spliceai_scores: Optional[Dict[str, Optional[float]]] = None,
) -> Dict[str, Any]:
    """Apply PS1 only from an approved reference; otherwise return a review candidate."""
    _load_references()
    result: Dict[str, Any] = {
        "applies": False,
        "code": "PS1",
        "strength": None,
        "points": 0,
        "reason": "",
        "reference_variant": None,
        "reference_status": "not_found",
        "application_status": "not_applied",
        "review_required": False,
        "blocking_reasons": [],
        "candidates": [],
        "vua_splice_evidence_status": vua_splice_evidence_status,
        "vua_splice_sources_checked": vua_splice_sources_checked or [],
        "vua_spliceai_score": spliceai_score,
        "vua_spliceai_source": vua_spliceai_source,
        "reference_spliceai_scores": reference_spliceai_scores or {},
        "source_url": vcep_specification(gene)["url"],
    }
    splice_low = spliceai_thresholds(gene)["bp4"]

    if (variant_type or "").lower() != "missense":
        result["reason"] = "Protein-level PS1 applies only to missense variants"
        return result
    aa_key = _extract_aa_change(p_notation)
    if not aa_key:
        result["reason"] = f"Could not extract a normalized missense change from {p_notation}"
        return result

    st7_matches = [
        item for item in _ST7_LOOKUP.get(gene, {}).get(aa_key, [])
        if item["c_notation"] != c_notation
    ]
    approved_records = [
        item for item in _APPROVED_LOOKUP.get(gene, {}).get(aa_key, [])
        if item["c_notation"] != c_notation
    ]
    candidates = _deduplicate_candidates(
        [*st7_matches, *(_public_approved_candidate(item) for item in approved_records)]
    )
    result["candidates"] = candidates
    if not candidates:
        result["reason"] = f"No P/LP reference with the same missense change {aa_key} was found"
        return result

    review_candidates = [
        item for item in candidates if item.get("reference_status") == "review_required"
    ]
    excluded_candidates = [
        item for item in candidates if item.get("reference_status") == "excluded"
    ]
    result["reference_status"] = (
        "approved"
        if approved_records
        else "review_required"
        if review_candidates
        else "excluded"
    )

    if spliceai_score is None:
        result["blocking_reasons"].append(
            "SpliceAI is unavailable for the variant under assessment"
        )
    elif spliceai_score > splice_low:
        result["reason"] = (
            f"Protein-level PS1 not applicable: assessed-variant SpliceAI "
            f"{spliceai_score:.3f} > {splice_low}"
        )
        result["reference_status"] = "ineligible_for_this_application"
        result["application_status"] = "not_applicable"
        return result

    if vua_splice_evidence_status == "abnormal":
        result["reason"] = (
            "Protein-level PS1 not applicable: a confirmed splice effect is recorded "
            "for the variant under assessment"
        )
        result["reference_status"] = "ineligible_for_this_application"
        result["application_status"] = "not_applicable"
        return result
    if vua_splice_evidence_status == "conflicting":
        result["blocking_reasons"].append(
            "Confirmed splice evidence for the variant under assessment is conflicting"
        )
    elif vua_splice_evidence_status not in NO_DAMAGING_SPLICE_STATUSES:
        result["blocking_reasons"].append(
            "The defined confirmed RNA/splice evidence sources have not been completely checked"
        )

    if not approved_records and review_candidates:
        result["blocking_reasons"].append(
            "At least one matching ST7 reference still requires protein-PS1 eligibility review"
        )
    elif not approved_records and excluded_candidates:
        result["application_status"] = "reference_ineligible"
        result["reason"] = "; ".join(
            item.get("status_reason") or "Matching ST7 reference is ineligible"
            for item in excluded_candidates
        )
        return result

    runtime_scores = reference_spliceai_scores or {}
    runtime_approved_records: List[Dict[str, Any]] = []
    unavailable_references: List[str] = []
    high_score_references: List[tuple[str, float]] = []
    for reference in approved_records:
        reference_c = str(reference["c_notation"])
        reference_score = runtime_scores.get(reference_c)
        if reference_score is None:
            unavailable_references.append(reference_c)
        elif reference_score > splice_low:
            high_score_references.append((reference_c, reference_score))
        else:
            runtime_approved_records.append(reference)

    if approved_records and unavailable_references:
        result["blocking_reasons"].append(
            "SpliceAI is unavailable for matching reference variant(s): "
            + ", ".join(unavailable_references)
        )

    if result["blocking_reasons"]:
        result["review_required"] = True
        result["application_status"] = "manual_review_required"
        result["reason"] = "; ".join(result["blocking_reasons"])
        return result

    if approved_records and not runtime_approved_records:
        result["reference_status"] = "ineligible_for_this_application"
        result["application_status"] = "not_applicable"
        result["reason"] = (
            "Protein-level PS1 not applicable: matching reference variant "
            f"SpliceAI exceeds {splice_low} ("
            + ", ".join(
                f"{reference_c}={score:.3f}"
                for reference_c, score in high_score_references
            )
            + ")"
        )
        return result

    best = max(
        runtime_approved_records,
        key=lambda item: 5 if item["classification"] == "Pathogenic" else 4,
    )
    best_reference_score = runtime_scores[str(best["c_notation"])]
    strength = "Strong" if best["classification"] == "Pathogenic" else "Moderate"
    result.update(
        {
            "applies": True,
            "strength": strength,
            "points": PS1_POINTS[strength],
            "reference_variant": _public_approved_candidate(best),
            "reference_status": "approved",
            "application_status": "auto_applied",
            "reason": (
                f"Same normalized missense change {aa_key} as approved "
                f"{best['classification']} reference {best['c_notation']} "
                f"{best['p_notation']} ({best['classification_source']}); both variants "
                "meet the recorded ENIGMA protein-level PS1 splice conditions "
                f"(assessed-variant SpliceAI {spliceai_score:.3f}, reference-variant "
                f"SpliceAI {best_reference_score:.3f}, source {vua_spliceai_source})"
            ),
        }
    )
    return result
