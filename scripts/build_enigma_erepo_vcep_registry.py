"""Build the local ENIGMA BRCA1/2 VCEP assertion registry from ClinGen ERepo.

The full ERepo export is used to discover records. Each selected record is then
resolved through the record endpoint so that the assertion-method identifier
and version are stored explicitly. Runtime classification never downloads this
data and never infers VCEP eligibility from ClinVar review stars.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import time
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "backend" / "data"
REGISTRY_PATH = DATA_DIR / "enigma_erepo_vcep_registry.json"
METADATA_PATH = DATA_DIR / "enigma_erepo_vcep_registry.metadata.json"

EXPORT_URL = (
    "https://erepo.clinicalgenome.org/evrepo/api/summary/"
    "classifications/download?type=csv"
)
SERVICE_URL = "https://erepo.clinicalgenome.org/evrepo/api/summary/srvc"
DETAIL_URL = "https://erepo.clinicalgenome.org/evrepo/api/classification/{uuid}"
PANEL_NAME = "ENIGMA BRCA1 and BRCA2 VCEP"
GUIDELINE_AFFILIATION = "50087"
REFERENCE_TRANSCRIPTS = {"BRCA1": "NM_007294.4", "BRCA2": "NM_000059.4"}
CURRENT_METHOD_VERSION = "1.2.0"
_VARIATION_RE = re.compile(
    r"^(NM_\d+\.\d+)(?:\((BRCA[12])\))?:(c\.[^ ]+)(?: \((p\.[^)]+)\))?$"
)


def _download_json(url: str, *, attempts: int = 4) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "ARIANE ERepo snapshot builder/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.load(response)
        except Exception as exc:  # network retry belongs only in the build tool
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Could not download {url}: {last_error}")


def _download_bytes(url: str, *, attempts: int = 4) -> bytes:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "ARIANE ERepo snapshot builder/1.0"},
            )
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read()
        except Exception as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Could not download {url}: {last_error}")


def _normalize_method_version(value: Any) -> str:
    text = str(value or "").strip()
    if re.fullmatch(r"\d+\.\d+", text):
        return f"{text}.0"
    return text


def _criterion_list(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def _detail(uuid: str) -> dict[str, Any]:
    return _download_json(DETAIL_URL.format(uuid=uuid))


def _record(export_row: dict[str, str], detail: dict[str, Any]) -> dict[str, Any]:
    uuid = export_row["Uuid"].strip()
    if str(detail.get("uuid") or "") != uuid:
        raise RuntimeError(f"ERepo detail UUID mismatch for {uuid}")
    match = _VARIATION_RE.fullmatch(export_row["Variation"].strip())
    if not match:
        raise RuntimeError(f"Unsupported preferred variation label: {export_row['Variation']!r}")
    transcript, label_gene, c_notation, p_notation = match.groups()
    gene = export_row["HGNC Gene Symbol"].strip()
    if label_gene and label_gene != gene:
        raise RuntimeError(
            f"ERepo record {uuid} gene mismatch: {label_gene!r} != {gene!r}"
        )
    if transcript != REFERENCE_TRANSCRIPTS[gene]:
        raise RuntimeError(
            f"ERepo record {uuid} uses {transcript}, expected {REFERENCE_TRANSCRIPTS[gene]}"
        )
    method = dict(detail.get("assertionMethod") or {})
    method_version = _normalize_method_version(method.get("version"))
    method_label = str(method.get("label") or "").strip()
    if method_label and (
        gene not in method_label or "ENIGMA BRCA1 and BRCA2" not in method_label
    ):
        raise RuntimeError(f"Unexpected assertion method for {uuid}: {method_label!r}")
    guideline = export_row["Guideline"].strip()
    if GUIDELINE_AFFILIATION not in guideline:
        raise RuntimeError(f"Unexpected guideline affiliation for {uuid}: {guideline!r}")

    classification = export_row["Assertion"].strip()
    outcome = str((detail.get("statementOutcome") or {}).get("label") or "").strip()
    if outcome != classification:
        raise RuntimeError(
            f"ERepo classification mismatch for {uuid}: {classification!r} != {outcome!r}"
        )
    retracted = export_row["Retracted"].strip().lower() == "true"
    return {
        "uuid": uuid,
        "gene": gene,
        "transcript": transcript,
        "c_notation": c_notation,
        "p_notation": f"p.({p_notation[2:]})" if p_notation else "",
        "classification": classification,
        "caid": export_row["Allele Registry Id"].strip(),
        "clinvar_variation_id": export_row["ClinVar Variation Id"].strip(),
        "disease": export_row["Disease"].strip(),
        "mondo_id": export_row["Mondo Id"].strip(),
        "mode_of_inheritance": export_row["Mode of Inheritance"].strip(),
        "criteria_met": _criterion_list(export_row["Applied Evidence Codes (Met)"]),
        "criteria_not_met": _criterion_list(
            export_row["Applied Evidence Codes (Not Met)"]
        ),
        "interpretation_summary": export_row["Summary of interpretation"].strip(),
        "pubmed_articles": _criterion_list(export_row["PubMed Articles"]),
        "expert_panel": export_row["Expert Panel"].strip(),
        "guideline_url": guideline,
        "assertion_method_id": str(method.get("@id") or "").strip(),
        "assertion_method_label": method_label,
        "assertion_method_version": method_version,
        "approval_date": export_row["Approval Date"].strip(),
        "published_date": export_row["Published Date"].strip(),
        "retracted": retracted,
        "erepo_url": export_row["Evidence Repo Link"].strip(),
        "source_status": (
            "current_vcep_assertion"
            if method_version == CURRENT_METHOD_VERSION and not retracted
            else "unversioned_vcep_assertion"
            if not method_version and not retracted
            else "historical_vcep_assertion"
        ),
    }


def _canonical_bytes(data: dict[str, Any]) -> bytes:
    return (json.dumps(data, ensure_ascii=False, indent=2) + "\n").encode("utf-8")


def build() -> tuple[dict[str, Any], dict[str, Any]]:
    exported = _download_bytes(EXPORT_URL)
    service = _download_json(SERVICE_URL)
    text = exported.decode("utf-8-sig")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    selected = [
        row
        for row in rows
        if row.get("Expert Panel", "").strip() == PANEL_NAME
        and row.get("HGNC Gene Symbol", "").strip() in REFERENCE_TRANSCRIPTS
    ]
    if not selected:
        raise RuntimeError("The ERepo export contains no ENIGMA BRCA1/2 VCEP records")

    details: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {
            executor.submit(_detail, row["Uuid"].strip()): row["Uuid"].strip()
            for row in selected
        }
        for future in as_completed(futures):
            uuid = futures[future]
            details[uuid] = future.result()

    records = [_record(row, details[row["Uuid"].strip()]) for row in selected]
    records.sort(key=lambda item: (item["gene"], item["c_notation"], item["uuid"]))
    if len({record["uuid"] for record in records}) != len(records):
        raise RuntimeError("The ERepo selection contains duplicate UUIDs")
    if any(record["retracted"] for record in records):
        raise RuntimeError("The ERepo selection contains a retracted assertion")

    classification_counts = Counter(record["classification"] for record in records)
    version_counts = Counter(record["assertion_method_version"] for record in records)
    status_counts = Counter(record["source_status"] for record in records)
    registry = {
        "schema_version": 1,
        "registry_version": datetime.now(timezone.utc).date().isoformat(),
        "status": "active",
        "description": (
            "Complete local snapshot of published ENIGMA BRCA1/2 VCEP assertions "
            "from ClinGen Evidence Repository. Current v1.2 assertions may support "
            "automatic workflows. Historical assertions are display-only."
        ),
        "expert_panel": PANEL_NAME,
        "guideline_affiliation_id": GUIDELINE_AFFILIATION,
        "current_assertion_method_version": CURRENT_METHOD_VERSION,
        "runtime_policy": {
            "primary_source": "local_checksum_validated_clingen_erepo_snapshot",
            "current_records": "eligible_for_rule_specific_runtime_review",
            "historical_or_unversioned_records": "display_only_manual_review_warning",
            "clinvar_review_stars": "display_only_not_a_vcep_eligibility_source",
            "network_fallback": False,
        },
        "record_count": len(records),
        "classification_counts": dict(sorted(classification_counts.items())),
        "assertion_method_version_counts": dict(sorted(version_counts.items())),
        "source_status_counts": dict(sorted(status_counts.items())),
        "records": records,
    }
    registry_bytes = _canonical_bytes(registry)
    metadata = {
        "schema_version": 1,
        "status": "active",
        "registry_file": REGISTRY_PATH.name,
        "registry_sha256": hashlib.sha256(registry_bytes).hexdigest(),
        "source_export_url": EXPORT_URL,
        "source_export_sha256": hashlib.sha256(exported).hexdigest(),
        "source_service_url": SERVICE_URL,
        "source_api_version": str((service.get("data") or {}).get("version") or ""),
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "expert_panel": PANEL_NAME,
        "guideline_affiliation_id": GUIDELINE_AFFILIATION,
        "record_count": len(records),
        "classification_counts": dict(sorted(classification_counts.items())),
        "assertion_method_version_counts": dict(sorted(version_counts.items())),
        "source_status_counts": dict(sorted(status_counts.items())),
    }
    return registry, metadata


def main() -> None:
    registry, metadata = build()
    REGISTRY_PATH.write_bytes(_canonical_bytes(registry))
    METADATA_PATH.write_bytes(_canonical_bytes(metadata))
    print(
        f"Wrote {registry['record_count']} ENIGMA ERepo assertions; "
        f"versions={registry['assertion_method_version_counts']}; "
        f"statuses={registry['source_status_counts']}"
    )


if __name__ == "__main__":
    main()
