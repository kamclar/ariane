"""Build the BRCA1/2 PP4/BP5 snapshot from the current ENIGMA data track.

The UCSC ENIGMA BRCAmfa track publishes one curator-combined clinical LR per
variant. ARIANE preserves that combined value as one evidence item and applies
the ENIGMA VCEP v1.2 thresholds. It never multiplies the publication-level
components again.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DEFAULT_SOURCE_MANIFEST = ROOT / "data/sources/enigma/clinical_lr_sources.manifest.json"
DEFAULT_OUTPUT = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.index.json"
DEFAULT_METADATA = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json"
INDEL_INDEX = ROOT / "data/precomputed/brca_normalized_indel_snapshot.index.json"
INDEL_METADATA = ROOT / "data/precomputed/brca_normalized_indel_snapshot.metadata.json"

TRANSCRIPTS = {"NM_007294.4": "BRCA1", "NM_000059.4": "BRCA2"}
BED_COLUMNS = (
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "thickStart", "thickEnd", "reserved", "combinedLR", "ACMGcode",
    "familyHistoryCombinedLR", "cooccurrenceCombinedLR",
    "segregationCombinedLR", "pathologyCombinedLR", "caseControlLR",
    "bridgesLR", "carriersLR", "ukbLR", "zantiSuggestedCode",
    "caputoLRs", "parsonsLRs", "liLRs", "eastonLRs", "HGVSp",
    "mouseOver",
)


@dataclass(frozen=True)
class SourceContext:
    manifest: dict
    dataset_id: str
    dataset: dict
    path: Path
    release_date: str
    track_url: str
    description_url: str
    build_url: str

CLINICAL_TYPE_COLUMNS = {
    "family_history": "familyHistoryCombinedLR",
    "cooccurrence": "cooccurrenceCombinedLR",
    "segregation": "segregationCombinedLR",
    "pathology": "pathologyCombinedLR",
    "case_control": "caseControlLR",
}

PUBLICATIONS = (
    {"pmid": "17924331", "citation": "Easton et al. 2007", "field": "eastonLRs"},
    {"pmid": "31131967", "citation": "Parsons et al. 2019", "field": "parsonsLRs"},
    {"pmid": "31853058", "citation": "Li et al. 2020", "field": "liLRs"},
    {"pmid": "34597585", "citation": "Caputo et al. 2021", "field": "caputoLRs"},
    {"pmid": "40413188", "citation": "Zanti et al. 2025", "field": "caseControlLR"},
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json_lf(path: Path, value: dict) -> None:
    """Write canonical UTF-8 JSON with LF endings on every operating system."""
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    path.write_bytes(payload.encode("utf-8"))


def _load_json_object(path: Path, description: str) -> dict:
    if not path.is_file():
        raise RuntimeError(f"{description} is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"{description} must contain a JSON object")
    return value


def load_source_manifest(
    manifest_path: Path,
    source: Path | None = None,
) -> SourceContext:
    manifest = _load_json_object(manifest_path, "Clinical LR source manifest")
    if manifest.get("schema_version") != 3:
        raise RuntimeError("Clinical LR source manifest schema is unsupported")
    if manifest.get("status") != "validated_source_manifest":
        raise RuntimeError("Clinical LR source manifest is not validated")
    if manifest.get("update_policy", {}).get("automatic_release_activation") is not False:
        raise RuntimeError("Clinical LR update policy must prohibit automatic activation")
    datasets = manifest.get("datasets") or {}
    if not isinstance(datasets, dict) or len(datasets) != 1:
        raise RuntimeError("Clinical LR source manifest must pin one active dataset")
    dataset_id, dataset = next(iter(datasets.items()))
    if not isinstance(dataset, dict):
        raise RuntimeError("Clinical LR source dataset contract is invalid")
    source_path = source or manifest_path.parent / str(dataset.get("file") or "")
    if not source_path.is_file():
        raise RuntimeError(f"Required clinical LR source is missing: {source_path}")
    if source_path.name != dataset.get("file"):
        raise RuntimeError("Clinical LR source filename does not match the manifest")
    actual_sha256 = sha256(source_path)
    if dataset.get("sha256") != actual_sha256:
        raise RuntimeError(
            "Clinical LR source checksum mismatch: "
            f"expected {dataset.get('sha256')}, found {actual_sha256}"
        )
    if dataset.get("schema") != "bed9+17" or dataset.get("columns") != list(BED_COLUMNS):
        raise RuntimeError("Clinical LR source schema does not match the validated contract")
    if dataset.get("combination_status") != "publisher_combined":
        raise RuntimeError("Clinical LR source is not declared as publisher-combined")
    release = manifest.get("clinical_lr_data_release") or {}
    required_release_fields = ("release_date", "description_url", "build_record_url")
    if any(not str(release.get(field) or "").strip() for field in required_release_fields):
        raise RuntimeError("Clinical LR data release provenance is incomplete")
    bigbed_path = manifest_path.parent / str(
        dataset.get("derived_from_bigbed_file") or ""
    )
    if not bigbed_path.is_file():
        raise RuntimeError(f"Pinned clinical LR BigBed is missing: {bigbed_path}")
    if sha256(bigbed_path) != dataset.get("derived_from_bigbed_sha256"):
        raise RuntimeError("Pinned clinical LR BigBed checksum mismatch")
    if bigbed_path.stat().st_size != int(dataset.get("remote_content_length") or -1):
        raise RuntimeError("Pinned clinical LR BigBed size does not match the manifest")
    return SourceContext(
        manifest=manifest,
        dataset_id=dataset_id,
        dataset=dataset,
        path=source_path,
        release_date=str(release.get("release_date") or ""),
        track_url=str(dataset.get("url") or ""),
        description_url=str(release.get("description_url") or ""),
        build_url=str(release.get("build_record_url") or ""),
    )


def _optional_float(value: str) -> float | None:
    text = str(value or "").strip()
    if not text or text.upper() in {"NULL", "NA", "N/A"}:
        return None
    parsed = float(text)
    if not math.isfinite(parsed) or parsed < 0:
        raise ValueError(f"Invalid likelihood ratio {text!r}")
    return parsed


def strength_for_lr(lr: float) -> tuple[str | None, str | None, int]:
    """Apply the exact ENIGMA VCEP v1.2 PP4/BP5 thresholds."""
    if lr >= 350:
        return "PP4", "Very Strong", 8
    if lr >= 18.7:
        return "PP4", "Strong", 4
    if lr >= 4.3:
        return "PP4", "Moderate", 2
    if lr >= 2.08:
        return "PP4", "Supporting", 1
    if lr <= 0.00285:
        return "BP5", "Very Strong", -8
    if lr <= 0.05:
        return "BP5", "Strong", -4
    if lr <= 0.23:
        return "BP5", "Moderate", -2
    if lr <= 0.48:
        return "BP5", "Supporting", -1
    return None, None, 0


def load_indel_reference() -> tuple[dict[str, dict], dict[str, str], dict[str, str]]:
    metadata = _load_json_object(INDEL_METADATA, "Normalized BRCA indel metadata")
    records = _load_json_object(INDEL_INDEX, "Normalized BRCA indel snapshot")
    if metadata.get("status") != "validated_reference_snapshot":
        raise RuntimeError("Normalized BRCA indel snapshot is not validated")
    if metadata.get("index_sha256") != sha256(INDEL_INDEX):
        raise RuntimeError("Normalized BRCA indel snapshot checksum mismatch")
    if metadata.get("records") != len(records):
        raise RuntimeError("Normalized BRCA indel snapshot record count mismatch")

    aliases: dict[str, str] = {key: key for key in records}
    ambiguous: set[str] = set()
    for canonical_key, record in records.items():
        for notation in record.get("input_c_notations", []):
            alias = f"{record['gene']}:{notation}"
            previous = aliases.get(alias)
            if previous and previous != canonical_key:
                ambiguous.add(alias)
            else:
                aliases[alias] = canonical_key
    if ambiguous:
        raise RuntimeError(
            f"Normalized BRCA indel snapshot has ambiguous aliases: {len(ambiguous)}"
        )
    return records, aliases, {
        "index_sha256": sha256(INDEL_INDEX),
        "metadata_sha256": sha256(INDEL_METADATA),
        "source_release": str(metadata.get("source_release", "")),
    }


def canonicalize_source_variant(
    gene: str,
    c_notation: str,
    indel_records: dict[str, dict],
    indel_aliases: dict[str, str],
) -> tuple[str, dict[str, str], bool]:
    from backend.modules.hgvs_engine import derive_protein_consequence

    normalized = derive_protein_consequence(gene, c_notation)
    canonical_c = normalized.canonical_c_notation
    matched_indel = False
    if re.search(r"delins|del|dup|ins", c_notation, re.IGNORECASE):
        submitted_key = f"{gene}:{c_notation}"
        normalized_key = f"{gene}:{canonical_c}"
        indel_key = indel_aliases.get(submitted_key) or indel_aliases.get(normalized_key)
        if indel_key:
            reference_c = indel_records[indel_key]["canonical_c_notation"]
            if reference_c != canonical_c:
                raise RuntimeError(
                    "HGVS normalization conflicts with normalized indel snapshot: "
                    f"{gene} {c_notation} -> {canonical_c}, snapshot -> {reference_c}"
                )
            matched_indel = True
    return canonical_c, normalized.provenance, matched_indel


def _publication_components(row: dict[str, str]) -> list[dict]:
    components = []
    for publication in PUBLICATIONS:
        raw = str(row[publication["field"]] or "").strip()
        if not raw or raw.upper() == "NULL":
            continue
        components.append({
            "pmid": publication["pmid"],
            "citation": publication["citation"],
            "reported_values": raw,
        })
    return components


def _source_component(
    row: dict[str, str],
    combined_lr: float,
    source_context: SourceContext,
) -> dict:
    clinical_data = []
    for data_type, field in CLINICAL_TYPE_COLUMNS.items():
        value = _optional_float(row[field])
        if value is not None:
            clinical_data.append({"data_type": data_type, "lr": value})
    publications = _publication_components(row)
    return {
        "source_id": source_context.dataset_id,
        "citation": (
            "UCSC ENIGMA BRCA1/BRCA2 likelihood track, "
            f"{source_context.release_date}"
        ),
        "pmid": ", ".join(item["pmid"] for item in publications),
        "pmids": [item["pmid"] for item in publications],
        "clinical_data": clinical_data,
        "component_lr": combined_lr,
        "evidence_family": "publisher_combined_clinical_lr",
        "independence_group": source_context.dataset_id,
        "source_bundle_id": source_context.dataset_id,
        "source_dataset": (
            "UCSC ENIGMA BRCAmfa track released "
            f"{source_context.release_date}"
        ),
        "publication_components": publications,
        "zanti_dataset_lrs": {
            "BRIDGES": _optional_float(row["bridgesLR"]),
            "CARRIERS": _optional_float(row["carriersLR"]),
            "UK_Biobank": _optional_float(row["ukbLR"]),
        },
        "zanti_suggested_code": str(row["zantiSuggestedCode"] or "").strip(),
    }


def _record_from_row(
    row: dict[str, str],
    canonical_c: str,
    submitted_c: str,
    source_context: SourceContext,
) -> dict:
    transcript = row["name"].split(":", 1)[0]
    gene = TRANSCRIPTS[transcript]
    combined_lr = _optional_float(row["combinedLR"])
    if combined_lr is None:
        raise RuntimeError(f"Published combined LR is missing for {row['name']}")
    code, strength, points = strength_for_lr(combined_lr)
    component = _source_component(row, combined_lr, source_context)
    return {
        "gene": gene,
        "reference_transcript": transcript,
        "canonical_c_notation": canonical_c,
        "input_c_notations": sorted({submitted_c, canonical_c}),
        "source_grch38_intervals": [{
            "chrom": row["chrom"].removeprefix("chr"),
            "start_0_based": int(row["chromStart"]),
            "end_0_based": int(row["chromEnd"]),
        }],
        "source_components": [component],
        "source_bundle_ids": [source_context.dataset_id],
        "source_bundle_count": 1,
        "independent_source_group_count": 1,
        "candidate_combined_lr": combined_lr,
        "combined_lr": combined_lr,
        "likelihood_ratio_status": (
            "source_reported_zero" if combined_lr == 0 else "available"
        ),
        "log10_combined_lr": math.log10(combined_lr) if combined_lr > 0 else None,
        "criterion": code,
        "strength": strength,
        "points": points,
        "informative": code is not None,
        "automatic_application_status": "eligible",
        "overlap_status": "source_curated_combination",
        "double_counting_risk": False,
        "source_reported_overlap_caveat": True,
        "automatic_combination_allowed": True,
        "assessment_note": (
            "ARIANE uses the combined LR published by the UCSC ENIGMA track as one "
            "evidence item and does not multiply its publication components again. "
            "The source replaced overlapping Parsons iCOGS case-control values with "
            "Zanti 2025 ccLR values. UCSC reports that some residual Parsons/Easton "
            "overlap may remain and recommends reviewing the component values."
        ),
        "assessment_sources": [
            source_context.description_url,
            source_context.build_url,
        ],
        "source": {
            "dataset": "UCSC ENIGMA BRCA1/BRCA2 likelihood for PP4 and BP5",
            "track_release": source_context.release_date,
            "track_url": source_context.track_url,
            "description_url": source_context.description_url,
        },
        "source_acmg_label": row["ACMGcode"],
        "source_hgvsp": row["HGVSp"],
    }


def build(
    source: Path | None = None,
    output: Path = DEFAULT_OUTPUT,
    metadata_path: Path = DEFAULT_METADATA,
    source_manifest_path: Path = DEFAULT_SOURCE_MANIFEST,
) -> dict:
    source_context = load_source_manifest(source_manifest_path, source)
    source_manifest = source_context.manifest
    source = source_context.path
    dataset_contract = source_context.dataset
    indel_records, indel_aliases, indel_dependency = load_indel_reference()

    records: dict[str, dict] = {}
    rows_seen = 0
    normalization_provenance: dict[str, str] | None = None
    normalization_counts = Counter()
    excluded = Counter()
    excluded_records = []
    normalization_failures = []

    with source.open("r", encoding="utf-8", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            rows_seen += 1
            if len(values) != len(BED_COLUMNS):
                raise RuntimeError(
                    f"Unexpected ENIGMA track column count on row {rows_seen}: "
                    f"expected {len(BED_COLUMNS)}, found {len(values)}"
                )
            row = dict(zip(BED_COLUMNS, values))
            name = row["name"].strip()
            transcript, separator, submitted_c = name.partition(":")
            if not separator or transcript not in TRANSCRIPTS or not submitted_c.startswith("c."):
                excluded["no_reference_transcript_hgvsc"] += 1
                excluded_records.append({
                    "source_name": name,
                    "reason": "The source row has genomic coordinates but no reference-transcript c. HGVS.",
                })
                continue
            gene = TRANSCRIPTS[transcript]
            try:
                canonical_c, provenance, matched_indel = canonicalize_source_variant(
                    gene, submitted_c, indel_records, indel_aliases
                )
            except (ValueError, RuntimeError) as exc:
                excluded["reference_transcript_hgvs_not_validated"] += 1
                failure = {
                    "source_name": name,
                    "reason": str(exc),
                    "error_code": str(getattr(exc, "code", type(exc).__name__)),
                }
                excluded_records.append(failure)
                normalization_failures.append(failure)
                continue
            if normalization_provenance is None:
                normalization_provenance = provenance
            elif normalization_provenance != provenance:
                raise RuntimeError("HGVS normalization provenance changed during snapshot build")
            normalization_counts["source_records_normalized"] += 1
            if canonical_c != submitted_c:
                normalization_counts["notations_canonicalized"] += 1
            if re.search(r"delins|del|dup|ins", submitted_c, re.IGNORECASE):
                normalization_counts[
                    "known_indels_cross_checked" if matched_indel
                    else "indels_not_in_reference_snapshot"
                ] += 1

            record = _record_from_row(
                row,
                canonical_c,
                submitted_c,
                source_context,
            )
            key = f"{gene}:{canonical_c}"
            previous = records.get(key)
            if previous is not None:
                if previous["combined_lr"] != record["combined_lr"]:
                    source_rows = previous.setdefault("conflicting_source_rows", [{
                        "source_name": previous["reference_transcript"] + ":" +
                        previous["input_c_notations"][0],
                        "combined_lr": previous["candidate_combined_lr"],
                    }])
                    source_rows.append({
                        "source_name": name,
                        "combined_lr": record["combined_lr"],
                    })
                    first_source_id = previous["source_components"][0]["source_id"]
                    if ":source_row_" not in first_source_id:
                        previous["source_components"][0]["source_id"] += ":source_row_1"
                    incoming_component = record["source_components"][0]
                    incoming_component["source_id"] += f":source_row_{len(source_rows)}"
                    previous["source_components"].append(incoming_component)
                    previous.update({
                        "candidate_combined_lr": None,
                        "combined_lr": None,
                        "likelihood_ratio_status": "unavailable_conflict",
                        "log10_combined_lr": None,
                        "criterion": None,
                        "strength": None,
                        "points": 0,
                        "informative": False,
                        "automatic_application_status": "review_required",
                        "overlap_status": "conflicting_normalized_source_rows",
                        "double_counting_risk": True,
                        "automatic_combination_allowed": False,
                        "independent_source_group_count": 0,
                        "assessment_note": (
                            "Multiple source rows normalize to the same reference-transcript "
                            "allele but report different combined LR values. ARIANE does not "
                            "choose or combine them automatically. Expert source review is required."
                        ),
                    })
                    continue
                previous["input_c_notations"] = sorted(set(
                    previous["input_c_notations"] + record["input_c_notations"]
                ))
                previous["source_grch38_intervals"] = sorted(
                    previous["source_grch38_intervals"] + record["source_grch38_intervals"],
                    key=lambda item: (item["chrom"], item["start_0_based"], item["end_0_based"]),
                )
                normalization_counts["canonical_source_rows_merged"] += 1
            else:
                records[key] = record

    expected_rows = int(dataset_contract["expected_rows"])
    if rows_seen != expected_rows:
        raise RuntimeError(
            f"Unexpected ENIGMA track row count: expected {expected_rows}, found {rows_seen}"
        )
    expected_without_hgvsc = int(dataset_contract["expected_rows_without_transcript_hgvsc"])
    if excluded["no_reference_transcript_hgvsc"] != expected_without_hgvsc:
        raise RuntimeError(
            "Unexpected number of ENIGMA rows without reference-transcript c. HGVS: "
            f"expected {expected_without_hgvsc}, found "
            f"{excluded['no_reference_transcript_hgvsc']}"
        )

    records = dict(sorted(records.items()))
    output.parent.mkdir(parents=True, exist_ok=True)
    _write_json_lf(output, records)
    criteria = Counter(record["criterion"] or "not_informative" for record in records.values())
    application_statuses = Counter(
        record["automatic_application_status"] for record in records.values()
    )
    likelihood_ratio_statuses = Counter(
        record["likelihood_ratio_status"] for record in records.values()
    )
    conflict_records = [
        record
        for record in records.values()
        if record["overlap_status"] == "conflicting_normalized_source_rows"
    ]
    conflict_source_row_count = sum(
        len(record.get("conflicting_source_rows", []))
        for record in conflict_records
    )
    normalization_conflicts = {
        "variant_count": len(conflict_records),
        "source_row_count": conflict_source_row_count,
        "excess_source_row_count": conflict_source_row_count - len(conflict_records),
    }
    evidence_type_counts = Counter(
        item["data_type"]
        for record in records.values()
        for item in record["source_components"][0]["clinical_data"]
    )
    metadata = {
        "dataset": "BRCA1/2 combined clinical likelihood-ratio snapshot",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "validated_derived_snapshot",
        "source_manifest_file": source_manifest_path.name,
        "source_manifest_sha256": sha256(source_manifest_path),
        "source_manifest": source_manifest,
        "source_files": {
            source_context.dataset_id: {
                "file": source.name,
                "sha256": sha256(source),
                "url": source_context.track_url,
                "description_url": source_context.description_url,
                "track_release": source_context.release_date,
                "combination_status": "publisher_combined",
            }
        },
        "target_rule_version": (
            "ENIGMA BRCA1/2 VCEP "
            f"{source_manifest['target_specification']['version']} PP4/BP5 thresholds"
        ),
        "clinical_lr_data_release": (
            f"UCSC ENIGMA BRCAmfa {source_context.release_date}"
        ),
        "reference_transcripts": TRANSCRIPTS,
        "rows_seen": rows_seen,
        "records": len(records),
        "criteria": dict(sorted(criteria.items())),
        "automatic_application_statuses": dict(sorted(application_statuses.items())),
        "likelihood_ratio_statuses": dict(sorted(likelihood_ratio_statuses.items())),
        "clinical_evidence_type_records": dict(sorted(evidence_type_counts.items())),
        "normalization_conflicts": normalization_conflicts,
        "combination_policy": {
            "method": "use publisher-combined LR without recomputing publication components",
            "single_final_code": True,
            "component_codes_scored_separately": False,
            "external_source_multiplication": "not allowed",
            "source_reported_residual_overlap_caveat": True,
        },
        "normalization": {
            "method": "biocommons.hgvs with checksum-pinned cdot panel provider",
            "provenance": normalization_provenance or {},
            "counts": dict(sorted(normalization_counts.items())),
            "normalized_indel_dependency": indel_dependency,
            "failures": normalization_failures,
        },
        "excluded": dict(sorted(excluded.items())),
        "excluded_source_records": excluded_records,
        "index_sha256": sha256(output),
    }
    _write_json_lf(metadata_path, metadata)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--source-manifest", type=Path, default=DEFAULT_SOURCE_MANIFEST)
    args = parser.parse_args()
    print(json.dumps(build(
        source=args.source,
        output=args.output,
        metadata_path=args.metadata,
        source_manifest_path=args.source_manifest,
    ), indent=2))


if __name__ == "__main__":
    main()
