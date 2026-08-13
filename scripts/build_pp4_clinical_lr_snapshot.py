"""Build the local BRCA1/2 PP4/BP5 clinical likelihood-ratio snapshot.

The source is the UCSC ENIGMA BRCAmfa hg38 track. Only publications named in
ENIGMA Appendix B are admitted. Caputo et al. values present in the UCSC track
are deliberately excluded from this snapshot until separately reviewed.
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
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DEFAULT_SOURCE = ROOT / "data/sources/enigma/BRCAmfa.hg38.v1.1.bed"
DEFAULT_OUTPUT = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.index.json"
DEFAULT_METADATA = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json"
INDEL_INDEX = ROOT / "data/precomputed/brca_normalized_indel_snapshot.index.json"
INDEL_METADATA = ROOT / "data/precomputed/brca_normalized_indel_snapshot.metadata.json"
TRANSCRIPTS = {"NM_007294.4": "BRCA1", "NM_000059.4": "BRCA2"}
TRACK_URL = "https://hgdownload.soe.ucsc.edu/hubs/enigma/hg38/BRCAmfa.bb"
TRACK_DESCRIPTION_URL = "https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html"

# Field order is Family history, Co-occurrence, Segregation, Pathology,
# Case-control, as documented by the UCSC item schema/description.
SOURCES = {
    "parsonsLRs": {
        "citation": "Parsons et al. 2019",
        "pmid": "31131967",
        "data_types": ["family_history", "cooccurrence", "segregation", "pathology", "case_control"],
    },
    "liLRs": {
        "citation": "Li et al. 2020",
        "pmid": "31853058",
        "data_types": ["personal_and_family_history"],
    },
    "eastonLRs": {
        "citation": "Easton et al. 2007",
        "pmid": "17924331",
        "data_types": ["family_history", "cooccurrence", "segregation", "pathology", "case_control"],
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_values(raw: str) -> list[float | None]:
    if not raw.strip():
        return []
    values = []
    for item in raw.strip().split(","):
        item = item.strip()
        values.append(None if item in {"", "NULL", "NA"} else float(item))
    return values


def strength_for_lr(lr: float) -> tuple[str | None, str | None, int]:
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


def apply_combined_evidence(record: dict) -> None:
    """Recompute the combined LR after independently sourced rows are merged."""
    components = sorted(record["source_components"], key=lambda item: item["pmid"])
    combined_lr = math.prod(component["component_lr"] for component in components)
    code, strength, points = strength_for_lr(combined_lr)
    record.update({
        "source_components": components,
        "combined_lr": combined_lr,
        "log10_combined_lr": math.log10(combined_lr) if combined_lr > 0 else None,
        "criterion": code,
        "strength": strength,
        "points": points,
        "informative": code is not None,
    })


def load_indel_reference() -> tuple[dict[str, dict], dict[str, str], dict[str, str]]:
    """Load the indel snapshot as a required, checksum-verified build input."""
    if not INDEL_INDEX.is_file() or not INDEL_METADATA.is_file():
        raise RuntimeError("Normalized BRCA indel snapshot or metadata is missing")
    metadata = json.loads(INDEL_METADATA.read_text(encoding="utf-8"))
    if metadata.get("status") != "validated_reference_snapshot":
        raise RuntimeError("Normalized BRCA indel snapshot is not validated")
    if metadata.get("index_sha256") != sha256(INDEL_INDEX):
        raise RuntimeError("Normalized BRCA indel snapshot checksum mismatch")
    records = json.loads(INDEL_INDEX.read_text(encoding="utf-8"))
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
    """Normalize against the local transcript and cross-check known indels.

    The HGVS engine verifies any explicitly supplied deleted or duplicated
    sequence against the checksum-pinned transcript. The indel snapshot is an
    independent cross-check, not a fallback source of an unverified alias.
    """
    from backend.modules.hgvs_engine import derive_protein_consequence

    normalized = derive_protein_consequence(gene, c_notation)
    canonical_c = normalized.canonical_c_notation
    matched_indel = False
    if re.search(r"delins|del|dup|ins", c_notation, re.IGNORECASE):
        submitted_key = f"{gene}:{c_notation}"
        normalized_key = f"{gene}:{canonical_c}"
        indel_key = indel_aliases.get(submitted_key) or indel_aliases.get(normalized_key)
        if indel_key:
            reference_record = indel_records[indel_key]
            reference_c = reference_record["canonical_c_notation"]
            if reference_c != canonical_c:
                raise RuntimeError(
                    "HGVS normalization conflicts with normalized indel snapshot: "
                    f"{gene} {c_notation} -> {canonical_c}, snapshot -> {reference_c}"
                )
            matched_indel = True
    return canonical_c, normalized.provenance, matched_indel


def build(source: Path, output: Path, metadata_path: Path) -> dict:
    indel_records, indel_aliases, indel_dependency = load_indel_reference()
    records: dict[str, dict] = {}
    conflicting_keys: set[str] = set()
    conflicting_records: dict[str, dict] = {}
    excluded: Counter[str] = Counter()
    normalization_counts: Counter[str] = Counter()
    normalization_failures: list[dict[str, str]] = []
    normalization_provenance: dict[str, str] | None = None
    rows_seen = 0

    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            rows_seen += 1
            name = row.get("name", "")
            if ":c." not in name:
                excluded["invalid_hgvs"] += 1
                continue
            transcript, c_notation = name.split(":", 1)
            gene = TRANSCRIPTS.get(transcript)
            if not gene:
                excluded["unsupported_transcript"] += 1
                continue

            try:
                canonical_c, provenance, matched_indel = canonicalize_source_variant(
                    gene, c_notation, indel_records, indel_aliases
                )
            except ValueError as exc:
                excluded[f"hgvs_normalization_failed:{type(exc).__name__}"] += 1
                normalization_failures.append({
                    "variant": name,
                    "error_code": str(getattr(exc, "code", type(exc).__name__)),
                    "reason": str(exc),
                })
                continue
            if normalization_provenance is None:
                normalization_provenance = provenance
            elif normalization_provenance != provenance:
                raise RuntimeError("HGVS normalization provenance changed during snapshot build")
            normalization_counts["source_records_normalized"] += 1
            if canonical_c != c_notation:
                normalization_counts["notations_canonicalized"] += 1
            if re.search(r"delins|del|dup|ins", c_notation, re.IGNORECASE):
                normalization_counts[
                    "known_indels_cross_checked" if matched_indel else "indels_not_in_reference_snapshot"
                ] += 1

            components = []
            all_values = []
            for field, definition in SOURCES.items():
                values = parse_values(row.get(field, ""))
                typed_values = [
                    {"data_type": definition["data_types"][idx], "lr": value}
                    for idx, value in enumerate(values[:len(definition["data_types"])])
                    if value is not None
                ]
                if not typed_values:
                    continue
                component_lr = math.prod(value["lr"] for value in typed_values)
                all_values.extend(value["lr"] for value in typed_values)
                components.append({
                    "citation": definition["citation"],
                    "pmid": definition["pmid"],
                    "clinical_data": typed_values,
                    "component_lr": component_lr,
                    "appendix_b_source": True,
                })
            if not all_values:
                excluded["no_appendix_b_lr"] += 1
                continue

            input_notations = {c_notation, canonical_c}
            source_interval = {
                "chrom": row["chrom"].removeprefix("chr"),
                "start_0_based": int(row["chromStart"]),
                "end_0_based": int(row["chromEnd"]),
            }

            key = f"{gene}:{canonical_c}"
            if key in conflicting_keys:
                excluded["conflicting_canonical_record"] += 1
                continue
            record = {
                "gene": gene,
                "reference_transcript": transcript,
                "canonical_c_notation": canonical_c,
                "input_c_notations": sorted(input_notations),
                "source_grch38_intervals": [source_interval],
                "source_components": components,
                "source": {
                    "dataset": "UCSC ENIGMA BRCAmfa track",
                    "track_version": "ENIGMA specifications 1.1.0",
                    "track_url": TRACK_URL,
                    "description_url": TRACK_DESCRIPTION_URL,
                    "excluded_source": "Caputo et al. 2021 (not listed in Appendix B v1.2)",
                },
            }
            apply_combined_evidence(record)
            previous = records.get(key)
            if previous:
                existing_by_pmid = {
                    component["pmid"]: component
                    for component in previous["source_components"]
                }
                component_conflict = any(
                    component["pmid"] in existing_by_pmid
                    and existing_by_pmid[component["pmid"]] != component
                    for component in record["source_components"]
                )
                if component_conflict:
                    excluded["conflicting_canonical_record"] += 1
                    conflicting_keys.add(key)
                    conflicting_records[key] = {
                        "reason": "different clinical LR components under the same PMID",
                        "existing_input_c_notations": previous["input_c_notations"],
                        "incoming_input_c_notations": record["input_c_notations"],
                        "existing_source_components": previous["source_components"],
                        "incoming_source_components": record["source_components"],
                    }
                    records.pop(key, None)
                    continue
                normalization_counts["canonical_source_rows_merged"] += 1
                previous["input_c_notations"] = sorted(set(
                    previous["input_c_notations"] + record["input_c_notations"]
                ))
                previous["source_grch38_intervals"] = sorted(
                    {
                        (item["chrom"], item["start_0_based"], item["end_0_based"])
                        for item in previous["source_grch38_intervals"]
                        + record["source_grch38_intervals"]
                    }
                )
                previous["source_grch38_intervals"] = [
                    {"chrom": chrom, "start_0_based": start, "end_0_based": end}
                    for chrom, start, end in previous["source_grch38_intervals"]
                ]
                new_components = [
                    component for component in record["source_components"]
                    if component["pmid"] not in existing_by_pmid
                ]
                normalization_counts["independent_source_components_merged"] += len(new_components)
                normalization_counts["duplicate_source_components_deduplicated"] += (
                    len(record["source_components"]) - len(new_components)
                )
                previous["source_components"].extend(new_components)
                apply_combined_evidence(previous)
                continue
            records[key] = record

    records = dict(sorted(records.items()))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes((json.dumps(records, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    criteria = Counter(record["criterion"] or "not_informative" for record in records.values())
    evidence_counts = Counter(
        component["pmid"]
        for record in records.values()
        for component in record["source_components"]
    )
    metadata = {
        "dataset": "BRCA1/2 Appendix B clinical likelihood-ratio snapshot",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "validated_derived_snapshot",
        "source_file": source.name,
        "source_sha256": sha256(source),
        "source_url": TRACK_URL,
        "source_description_url": TRACK_DESCRIPTION_URL,
        "source_track_version": "ENIGMA specifications 1.1.0",
        "target_rule_version": "ENIGMA BRCA1/2 VCEP 1.2 PP4/BP5 thresholds",
        "reference_transcripts": TRANSCRIPTS,
        "rows_seen": rows_seen,
        "records": len(records),
        "criteria": dict(sorted(criteria.items())),
        "records_by_appendix_b_pmid": dict(sorted(evidence_counts.items())),
        "included_appendix_b_sources": SOURCES,
        "normalization": {
            "method": "biocommons.hgvs with checksum-pinned cdot panel provider",
            "provenance": normalization_provenance or {},
            "counts": dict(sorted(normalization_counts.items())),
            "normalized_indel_dependency": indel_dependency,
            "failures": normalization_failures,
        },
        "excluded": dict(sorted(excluded.items())),
        "conflicting_canonical_keys": sorted(conflicting_keys),
        "conflicting_canonical_records": dict(sorted(conflicting_records.items())),
        "index_sha256": sha256(output),
    }
    metadata_path.write_bytes((json.dumps(metadata, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    args = parser.parse_args()
    print(json.dumps(build(args.source, args.output, args.metadata), indent=2))


if __name__ == "__main__":
    main()
