"""Build generic BRCA exon-CNV population data from official sources.

No variant or criterion is enrolled by hand. Exons come from the lossless
ENIGMA Table 4 snapshot, GRCh37 coordinates come from the versioned local
coordinate snapshot, and population observations come from gnomAD-SV v2.1.
The runtime module applies the ENIGMA Appendix G decision path.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "sources" / "enigma" / "exon_cnv_evidence_manifest.json"
DEFAULT_OUTPUT = ROOT / "backend" / "data" / "exon_cnv_evidence.json"


def canonical_sha256(value) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def obtain_source(url: str, source: Path | None):
    if source is not None:
        source = source.resolve()
        if not source.is_file():
            raise RuntimeError(f"gnomAD-SV source does not exist: {source}")
        return source, {"url": url, "retrieval": "provided_local_file"}, None

    temporary = tempfile.TemporaryDirectory(prefix="ariane-gnomad-sv-")
    target = Path(temporary.name) / "gnomad_sv.bed.gz"
    request = urllib.request.Request(url, headers={"User-Agent": "ARIANE-data-builder/2"})
    with urllib.request.urlopen(request, timeout=120) as response, target.open("wb") as output:
        headers = {
            "url": response.geturl(),
            "retrieval": "downloaded",
            "etag": response.headers.get("ETag"),
            "x_goog_hash": response.headers.get("x-goog-hash"),
            "last_modified": response.headers.get("Last-Modified"),
            "content_length": response.headers.get("Content-Length"),
        }
        shutil.copyfileobj(response, output, length=1024 * 1024)
    return target, headers, temporary


def spans_interval(row: dict, interval: dict) -> bool:
    if row.get("svtype") != "DEL":
        return False
    chrom = str(row.get("#chrom", "")).removeprefix("chr")
    return (
        chrom == str(interval["chrom"]).removeprefix("chr")
        and int(row["start"]) <= int(interval["start_1_based_inclusive"]) - 1
        and int(row["end"]) >= int(interval["end_1_based_inclusive"])
    )


def spans_all_exons(row: dict, intervals: list[dict]) -> bool:
    return bool(intervals) and all(spans_interval(row, interval) for interval in intervals)


def compact_match(row: dict) -> dict:
    fields = ("#chrom", "start", "end", "name", "svtype", "SVLEN", "AN", "AC", "AF", "FILTER")
    return {field: row.get(field) for field in fields}


def _coordinate_by_c_position(coordinates: dict, gene: str) -> dict[int, tuple[str, int]]:
    result: dict[int, tuple[str, int]] = {}
    pattern = re.compile(r"^c\.(\d+)[ACGT]>[ACGT]$")
    prefix = f"{gene}:"
    for key, record in coordinates.items():
        if not key.startswith(prefix):
            continue
        match = pattern.match(record.get("c_notation", ""))
        grch37 = record.get("grch37") or ""
        if not match or not grch37:
            continue
        parts = grch37.split(":")
        if len(parts) != 3:
            continue
        result.setdefault(int(match.group(1)), (parts[0], int(parts[1])))
    return result


def derive_exons(manifest: dict) -> dict:
    exon_source = manifest["exon_source"]
    coordinate_source = manifest["coordinate_source"]
    table4_path = ROOT / exon_source["path"]
    coordinate_path = ROOT / coordinate_source["path"]
    for label, path, expected in (
        ("ENIGMA Table 4", table4_path, exon_source["sha256"]),
        ("coordinate snapshot", coordinate_path, coordinate_source["sha256"]),
    ):
        if not path.is_file() or file_sha256(path) != expected:
            raise RuntimeError(f"{label} is missing or its checksum changed")

    table4 = json.loads(table4_path.read_text(encoding="utf-8"))
    coordinates = json.loads(coordinate_path.read_text(encoding="utf-8"))
    exons = {}
    for gene, ranges in table4[exon_source["section"]].items():
        position_map = _coordinate_by_c_position(coordinates, gene)
        for exon, (c_start, c_end) in ranges.items():
            coding_start = max(1, int(c_start))
            coding_end = int(c_end)
            key = f"{gene}:{exon}"
            record = {
                "gene": gene,
                "exon": exon,
                "c_range": [c_start, c_end],
                "coding_length_bp": max(0, coding_end - coding_start + 1),
                "coordinate_status": "unavailable",
                "grch37_coding_interval": None,
                "all_matching_deletions": [],
                "pass_matching_deletions": [],
            }
            start_coordinate = position_map.get(coding_start)
            end_coordinate = position_map.get(coding_end)
            if coding_start <= coding_end and start_coordinate and end_coordinate:
                if start_coordinate[0] != end_coordinate[0]:
                    raise RuntimeError(f"Exon endpoints map to different chromosomes: {key}")
                record["coordinate_status"] = "ok"
                record["grch37_coding_interval"] = {
                    "chrom": start_coordinate[0],
                    "start_1_based_inclusive": min(start_coordinate[1], end_coordinate[1]),
                    "end_1_based_inclusive": max(start_coordinate[1], end_coordinate[1]),
                }
            exons[key] = record
    return exons


def build(manifest_path: Path, output_path: Path, source: Path | None) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 2:
        raise RuntimeError("Exon-CNV population manifest is missing or invalid")
    exons = derive_exons(manifest)
    population = manifest["population_source"]
    source_path, source_identity, temporary = obtain_source(population["url"], source)
    try:
        source_identity["sha256"] = file_sha256(source_path)
        source_identity["bytes"] = source_path.stat().st_size
        expected_source = population["expected_source_identity"]
        for field in ("bytes", "sha256"):
            if source_identity.get(field) != expected_source.get(field):
                raise RuntimeError(
                    f"gnomAD-SV source {field} mismatch: expected "
                    f"{expected_source.get(field)!r}, found {source_identity.get(field)!r}"
                )
        if source_identity.get("etag") and expected_source.get("etag") and (
            source_identity["etag"] != expected_source["etag"]
        ):
            raise RuntimeError("gnomAD-SV source ETag mismatch")

        rows_scanned = 0
        deletion_rows_scanned = 0
        available = [record for record in exons.values() if record["coordinate_status"] == "ok"]
        with gzip.open(source_path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"#chrom", "start", "end", "name", "svtype", "FILTER"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise RuntimeError("gnomAD-SV BED has an unexpected schema")
            for row in reader:
                rows_scanned += 1
                if row.get("svtype") != "DEL":
                    continue
                deletion_rows_scanned += 1
                for record in available:
                    if spans_interval(row, record["grch37_coding_interval"]):
                        match = compact_match(row)
                        record["all_matching_deletions"].append(match)
                        if row.get("FILTER") == "PASS":
                            record["pass_matching_deletions"].append(match)

        payload = {
            "schema_version": 2,
            "snapshot_id": manifest["manifest_id"],
            "rule_set": manifest["rule_set"],
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "builder": "scripts/build_exon_cnv_evidence_snapshot.py",
            "manifest_sha256": file_sha256(manifest_path),
            "source_identity": source_identity,
            "rows_scanned": rows_scanned,
            "deletion_rows_scanned": deletion_rows_scanned,
            "pm2_policy": population["pm2_policy"],
            "exons": exons,
        }
        payload["exons_sha256"] = canonical_sha256(exons)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        return payload
    finally:
        if temporary is not None:
            temporary.cleanup()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--source", type=Path, help="Existing official gnomAD-SV BED.gz")
    args = parser.parse_args()
    payload = build(args.manifest, args.output, args.source)
    print(
        f"Wrote {args.output}: {len(payload['exons'])} exon record(s), "
        f"source sha256={payload['source_identity']['sha256']}"
    )


if __name__ == "__main__":
    main()
