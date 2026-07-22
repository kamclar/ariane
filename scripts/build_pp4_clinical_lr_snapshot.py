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
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data/sources/enigma/BRCAmfa.hg38.v1.1.bed"
DEFAULT_OUTPUT = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.index.json"
DEFAULT_METADATA = ROOT / "data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json"
INDEL_INDEX = ROOT / "data/precomputed/brca_normalized_indel_snapshot.index.json"
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


def indel_aliases() -> dict[str, set[str]]:
    aliases: dict[str, set[str]] = {}
    if not INDEL_INDEX.is_file():
        return aliases
    records = json.loads(INDEL_INDEX.read_text(encoding="utf-8"))
    for key, record in records.items():
        gene = record["gene"]
        names = set(record.get("input_c_notations", [])) | {record["canonical_c_notation"]}
        for notation in names:
            aliases.setdefault(f"{gene}:{notation}", set()).add(key)
    return aliases


def build(source: Path, output: Path, metadata_path: Path) -> dict:
    alias_to_indel = indel_aliases()
    records: dict[str, dict] = {}
    conflicting_keys: set[str] = set()
    excluded: Counter[str] = Counter()
    evidence_counts: Counter[str] = Counter()
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
                evidence_counts[definition["pmid"]] += 1

            if not all_values:
                excluded["no_appendix_b_lr"] += 1
                continue
            combined_lr = math.prod(all_values)
            code, strength, points = strength_for_lr(combined_lr)

            input_notations = {c_notation}
            alias_key = f"{gene}:{c_notation}"
            indel_targets = alias_to_indel.get(alias_key, set())
            canonical_c = c_notation
            if len(indel_targets) == 1:
                canonical_key = next(iter(indel_targets))
                canonical_c = canonical_key.split(":", 1)[1]
                input_notations.add(canonical_c)
            elif len(indel_targets) > 1:
                excluded["ambiguous_indel_alias"] += 1
                continue

            key = f"{gene}:{canonical_c}"
            if key in conflicting_keys:
                excluded["conflicting_canonical_record"] += 1
                continue
            record = {
                "gene": gene,
                "reference_transcript": transcript,
                "canonical_c_notation": canonical_c,
                "input_c_notations": sorted(input_notations),
                "grch38": {
                    "chrom": row["chrom"].removeprefix("chr"),
                    "start_0_based": int(row["chromStart"]),
                    "end_0_based": int(row["chromEnd"]),
                },
                "combined_lr": combined_lr,
                "log10_combined_lr": math.log10(combined_lr) if combined_lr > 0 else None,
                "criterion": code,
                "strength": strength,
                "points": points,
                "informative": code is not None,
                "source_components": components,
                "source": {
                    "dataset": "UCSC ENIGMA BRCAmfa track",
                    "track_version": "ENIGMA specifications 1.1.0",
                    "track_url": TRACK_URL,
                    "description_url": TRACK_DESCRIPTION_URL,
                    "excluded_source": "Caputo et al. 2021 (not listed in Appendix B v1.2)",
                },
            }
            previous = records.get(key)
            if previous:
                identity_fields = (
                    "grch38", "combined_lr", "criterion", "strength",
                    "points", "source_components",
                )
                if any(previous[field] != record[field] for field in identity_fields):
                    excluded["conflicting_canonical_record"] += 1
                    conflicting_keys.add(key)
                    records.pop(key, None)
                    continue
                previous["input_c_notations"] = sorted(set(
                    previous["input_c_notations"] + record["input_c_notations"]
                ))
                continue
            records[key] = record

    records = dict(sorted(records.items()))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes((json.dumps(records, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    criteria = Counter(record["criterion"] or "not_informative" for record in records.values())
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
        "excluded": dict(sorted(excluded.items())),
        "conflicting_canonical_keys": sorted(conflicting_keys),
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
