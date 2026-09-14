"""Validation helpers for the curated protein-level PS1 registry.

This module owns the registry schema and integrity checks.  It deliberately
does not import the PS1 criterion evaluator, so reference-data validation can
run before rule modules are constructed without creating a package cycle.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Dict, Optional

from backend.policy.gene import active_genes, reference_transcript, spliceai_thresholds


KNOWN_CLASSIFICATION_VERIFICATIONS = {
    "external_vcep_assertion",
    "locally_recurated_under_enigma_vcep",
    "enigma_st7_v1_2_reference_set",
}
AUTOMATIC_CLASSIFICATION_VERIFICATIONS = {
    "external_vcep_assertion",
    "locally_recurated_under_enigma_vcep",
    "enigma_st7_v1_2_reference_set",
}
NO_DAMAGING_SPLICE_STATUSES = {"none_identified", "normal"}


def extract_missense_aa_change(p_notation: str) -> Optional[str]:
    """Return a canonical three-letter missense key, for example Arg170Gln."""
    if not p_notation:
        return None
    clean = p_notation.replace("p.", "").replace("(", "").replace(")", "").strip()
    if any(token in clean.lower() for token in ("fs", "ter", "del", "ins", "dup", "=", "?")):
        return None
    return clean if re.fullmatch(r"[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", clean) else None


def _validate_sha256(value: Any, label: str) -> None:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise RuntimeError(f"PS1 reference registry has invalid {label}")


def compute_approval_basis_checksum(record: Dict[str, Any]) -> str:
    """Checksum every approved-record field except the checksum itself."""
    basis = {key: value for key, value in record.items() if key != "approval_basis_checksum"}
    canonical = json.dumps(
        basis,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def validate_ps1_reference_registry(data: Dict[str, Any]) -> None:
    """Validate the curated automatic-scoring registry and known dependencies."""
    if data.get("schema_version") != 6 or data.get("status") != "active":
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
    } != KNOWN_CLASSIFICATION_VERIFICATIONS:
        raise RuntimeError("PS1 protein reference registry has invalid reference source policy")
    decision = data.get("methodological_decision")
    if not isinstance(decision, dict) or any(
        not str(decision.get(field) or "").strip()
        for field in ("id", "date", "status", "decision", "scope")
    ):
        raise RuntimeError("PS1 protein reference registry lacks its methodological decision")
    source_checksums = data.get("source_checksums")
    required_checksum_keys = {
        "st7_sha256",
        "table9_sha256",
        "st2_sha256",
        "curated_extensions_sha256",
        "erepo_vcep_registry_sha256",
        "erepo_vcep_metadata_sha256",
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
    configured_genes = set(active_genes())
    for index, record in enumerate(records):
        prefix = f"PS1 approved reference #{index}"
        reference_id = str(record.get("reference_id") or "")
        gene = record.get("gene")
        c_notation = str(record.get("c_notation") or "")
        p_notation = str(record.get("p_notation") or "")
        if not reference_id or reference_id in by_id:
            raise RuntimeError(f"{prefix} has missing or duplicate reference_id")
        if gene not in configured_genes or not c_notation.startswith("c."):
            raise RuntimeError(f"{prefix} has invalid variant identity")
        if not extract_missense_aa_change(p_notation):
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
        if verification not in KNOWN_CLASSIFICATION_VERIFICATIONS:
            raise RuntimeError(f"{prefix} has an unknown classification verification")
        if status == "eligible" and verification not in AUTOMATIC_CLASSIFICATION_VERIFICATIONS:
            raise RuntimeError(
                f"{prefix} is eligible without an accepted reference classification basis"
            )
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
        elif status == "excluded" and confirmed_status != "abnormal":
            raise RuntimeError(f"{prefix} has no recorded reason for exclusion")
        _validate_sha256(
            record.get("approval_basis_checksum"),
            f"approval_basis_checksum for {reference_id}",
        )
        if record["approval_basis_checksum"] != compute_approval_basis_checksum(record):
            raise RuntimeError(f"{prefix} approval_basis_checksum does not match its content")

        dependency = record.get("classification_ps1_dependency", {})
        used = dependency.get("used", "unknown")
        dependencies = dependency.get("reference_ids", [])
        if used not in {True, False, "unknown"} or not isinstance(dependencies, list):
            raise RuntimeError(f"{prefix} has invalid PS1 dependency metadata")
        if status == "eligible" and used not in {True, False}:
            raise RuntimeError(f"{prefix} is eligible without a resolved PS1 dependency review")
        if used is True and not dependencies:
            raise RuntimeError(f"{prefix} used PS1 but does not identify its reference")
        if used is not True and dependencies:
            raise RuntimeError(
                f"{prefix} lists PS1 dependencies without recording that PS1 was used"
            )
        if reference_id in dependencies:
            raise RuntimeError(f"{prefix} directly depends on itself")
        by_id[reference_id] = record

    visiting = set()
    visited = set()

    def visit(reference_id: str) -> None:
        if reference_id in visiting:
            raise RuntimeError(
                f"PS1 reference registry contains a circular dependency at {reference_id}"
            )
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
