"""Generic ENIGMA Appendix G evaluation for exon-level BRCA CNVs."""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict

from backend.config import EXON_CNV_EVIDENCE_PATH, EXON_CNV_EVIDENCE_MANIFEST_PATH
from backend.gene_policy import active_genes
from backend.modules.table4 import parse_exon_from_deletion_notation


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def validate_exon_cnv_evidence_snapshot(
    payload: Dict[str, Any], manifest_path: Path = EXON_CNV_EVIDENCE_MANIFEST_PATH
) -> None:
    if payload.get("schema_version") != 2:
        raise RuntimeError("Required exon-CNV population snapshot has an unsupported schema")
    if payload.get("builder") != "scripts/build_exon_cnv_evidence_snapshot.py":
        raise RuntimeError("Required exon-CNV population snapshot has unknown provenance")
    if not manifest_path.is_file():
        raise RuntimeError(f"Required exon-CNV population manifest is missing: {manifest_path}")
    if payload.get("manifest_sha256") != _file_sha256(manifest_path):
        raise RuntimeError("Required exon-CNV population snapshot manifest checksum mismatch")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Required exon-CNV population manifest cannot be loaded: {exc}") from exc

    expected_source = manifest["population_source"]["expected_source_identity"]
    source = payload.get("source_identity") or {}
    if (
        source.get("bytes") != expected_source.get("bytes")
        or source.get("sha256") != expected_source.get("sha256")
        or source.get("url") != manifest["population_source"]["url"]
        or (
            source.get("etag")
            and expected_source.get("etag")
            and source.get("etag") != expected_source.get("etag")
        )
    ):
        raise RuntimeError("Required exon-CNV population source identity is invalid")

    exons = payload.get("exons")
    if not isinstance(exons, dict) or len(exons) != 50:
        raise RuntimeError("Required exon-CNV population snapshot must cover all 50 Table 4 exons")
    if payload.get("exons_sha256") != _canonical_sha256(exons):
        raise RuntimeError("Required exon-CNV population snapshot exon checksum mismatch")
    if payload.get("pm2_policy") != manifest["population_source"]["pm2_policy"]:
        raise RuntimeError("Required exon-CNV population snapshot policy mismatch")

    for key, record in exons.items():
        if key != f"{record.get('gene')}:{record.get('exon')}":
            raise RuntimeError(f"Invalid exon-CNV population key: {key}")
        if record.get("gene") not in set(active_genes()):
            raise RuntimeError(f"Invalid exon-CNV gene: {key}")
        if record.get("coordinate_status") not in {"ok", "unavailable"}:
            raise RuntimeError(f"Invalid exon-CNV coordinate status: {key}")
        if record.get("coordinate_status") == "ok":
            interval = record.get("grch37_coding_interval") or {}
            if not all(interval.get(field) is not None for field in (
                "chrom", "start_1_based_inclusive", "end_1_based_inclusive"
            )):
                raise RuntimeError(f"Incomplete exon-CNV interval: {key}")
        for match in record.get("pass_matching_deletions", []):
            if match not in record.get("all_matching_deletions", []):
                raise RuntimeError(f"PASS exon-CNV match is absent from full match set: {key}")


@lru_cache(maxsize=1)
def load_exon_cnv_evidence_snapshot() -> Dict[str, Any]:
    if not EXON_CNV_EVIDENCE_PATH.is_file():
        raise RuntimeError(
            f"Required exon-CNV population snapshot is missing: {EXON_CNV_EVIDENCE_PATH}"
        )
    try:
        payload = json.loads(EXON_CNV_EVIDENCE_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Required exon-CNV population snapshot cannot be loaded: {exc}"
        ) from exc
    validate_exon_cnv_evidence_snapshot(payload)
    return payload


def lookup_exon_cnv_evidence(gene: str, c_notation: str) -> Dict[str, Any]:
    """Run the general Appendix G path; never look up the exact variant."""
    payload = load_exon_cnv_evidence_snapshot()
    is_deletion = c_notation.endswith("del") and not c_notation.endswith("delins")
    trace = [{
        "step": "variant_type",
        "status": "pass" if is_deletion else "fail",
        "detail": "exon-level deletion" if is_deletion else "not an exon-level deletion",
    }]
    if not is_deletion:
        return {
            "found": False,
            "criteria": [],
            "decision_trace": trace,
            "reason": (
                "PM2 was not applied: the current ENIGMA Appendix G population "
                "path is validated for exon-level deletions, not this CNV type."
            ),
            "snapshot_id": payload["snapshot_id"],
        }
    exon = parse_exon_from_deletion_notation(c_notation, gene)
    if exon is None:
        trace.append({
            "step": "table4_exon_mapping",
            "status": "fail",
            "detail": "notation does not identify one complete Table 4 exon",
        })
        return {
            "found": False,
            "criteria": [],
            "decision_trace": trace,
            "reason": (
                "PM2 was not applied: the exon-CNV notation could not be mapped "
                "unambiguously to one complete ENIGMA Table 4 exon."
            ),
            "snapshot_id": payload["snapshot_id"],
        }
    trace.append({"step": "table4_exon_mapping", "status": "pass", "detail": exon})

    record = payload["exons"].get(f"{gene}:{exon}")
    if record is None or record.get("coordinate_status") != "ok":
        trace.append({
            "step": "grch37_exon_interval",
            "status": "fail",
            "detail": "reproducible exon interval unavailable",
        })
        return {
            "found": True,
            "criteria": [],
            "decision_trace": trace,
            "reason": "PM2 was not applied: the GRCh37 exon interval is unavailable.",
            "exon": exon,
            "snapshot_id": payload["snapshot_id"],
        }
    trace.append({
        "step": "grch37_exon_interval", "status": "pass",
        "detail": record["grch37_coding_interval"],
    })

    policy = payload["pm2_policy"]
    eligible_size = int(record["coding_length_bp"]) >= int(policy["minimum_variant_size_bp"])
    trace.append({
        "step": "appendix_g_size",
        "status": "pass" if eligible_size else "fail",
        "detail": f"minimum deleted coding sequence {record['coding_length_bp']} bp",
    })
    if not eligible_size:
        return {
            "found": True,
            "criteria": [],
            "decision_trace": trace,
            "reason": "PM2 was not applied: the deletion is not >50 bp.",
            "exon": exon,
            "snapshot_id": payload["snapshot_id"],
        }

    matches = record.get("all_matching_deletions", [])
    pass_matches = record.get("pass_matching_deletions", [])
    absence_established = not matches
    trace.append({
        "step": "gnomad_sv_exon_match",
        "status": "pass" if absence_established else "fail",
        "detail": {
            "all_matching_deletions": len(matches),
            "pass_matching_deletions": len(pass_matches),
        },
    })
    criteria = []
    if absence_established:
        criteria.append({
            "code": policy["criterion"],
            "strength": policy["strength"],
            "points": policy["points"],
            "reason": (
                f"ENIGMA Appendix G: this >50 bp {gene} {exon} deletion is absent "
                "from gnomAD-SV v2.1 using exon-overlap matching; PM2 Supporting applies."
            ),
            "source": policy["source_url"],
        })
        reason = "PM2 Supporting was assigned by the ENIGMA Appendix G decision path."
    else:
        reason = (
            "PM2 was not applied: gnomAD-SV contains a deletion spanning the "
            f"complete coding interval of {gene} {exon}."
        )
    return {
        "found": True,
        "criteria": criteria,
        "decision_trace": trace,
        "reason": reason,
        "exon": exon,
        "population_evidence": {
            "dataset": "gnomAD-SV v2.1",
            "absence_established": absence_established,
            "matching_records": matches,
            "matching_pass_records": pass_matches,
        },
        "snapshot_id": payload["snapshot_id"],
        "source_identity": payload["source_identity"],
    }
