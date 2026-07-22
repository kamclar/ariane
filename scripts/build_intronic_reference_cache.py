"""Build versioned BRCA1/2 intronic SNV coordinate and SpliceAI caches.

The manifest covers every SNV within ``--window`` bases of each internal exon
boundary on the reference transcripts. Exon boundaries are derived from the
versioned coding-SNV snapshot; reference bases come from the UCSC sequence API.
SpliceAI scoring is resumable and uses the same Broad-compatible API contract
as the application.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import re
import threading
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "data/precomputed/brca_module1_snv_classification_snapshot.index.json"
COORDINATES = ROOT / "data/coordinates/brca_intronic_snv_coordinates.json"
COORDINATE_METADATA = ROOT / "data/coordinates/brca_intronic_snv_coordinates.metadata.json"
SPLICEAI = ROOT / "data/spliceai/spliceai_brca_intronic_snv_reference_cache.json"
SPLICEAI_METADATA = ROOT / "data/spliceai/spliceai_brca_intronic_snv_reference_cache.metadata.json"
DEFAULT_API = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/"
TRANSCRIPTS = {
    "BRCA1": {"refseq": "NM_007294.4", "ensembl": "ENST00000357654.9"},
    "BRCA2": {"refseq": "NM_000059.4", "ensembl": "ENST00000380152.8"},
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class CodingPosition:
    c_pos: int
    chrom: str
    pos: int
    ref: str


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_coord(value: str) -> tuple[str, int, str, str]:
    match = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([ACGT]+)>([ACGT]+)", value or "")
    if not match:
        raise ValueError(f"Invalid snapshot coordinate: {value!r}")
    return match.group(1), int(match.group(2)), match.group(3), match.group(4)


def coding_positions(snapshot: dict, gene: str, assembly_key: str) -> list[CodingPosition]:
    positions: dict[int, CodingPosition] = {}
    pattern = re.compile(rf"^{gene}:c\.(\d+)([ACGT])>[ACGT]$")
    for key, record in snapshot.items():
        match = pattern.match(key)
        if not match:
            continue
        c_pos = int(match.group(1))
        chrom, pos, ref, _ = _parse_coord(record[assembly_key])
        item = CodingPosition(c_pos, chrom, pos, ref)
        previous = positions.setdefault(c_pos, item)
        if previous != item:
            raise ValueError(f"Inconsistent coordinates for {gene} c.{c_pos}")
    result = sorted(positions.values(), key=lambda item: item.c_pos)
    if not result:
        raise ValueError(f"No coding positions for {gene}/{assembly_key}")
    return result


def exon_blocks(positions: list[CodingPosition]) -> tuple[int, list[list[CodingPosition]]]:
    steps = [b.pos - a.pos for a, b in zip(positions, positions[1:]) if b.c_pos == a.c_pos + 1]
    strand = 1 if steps.count(1) > steps.count(-1) else -1
    blocks: list[list[CodingPosition]] = [[positions[0]]]
    for item in positions[1:]:
        previous = blocks[-1][-1]
        if item.c_pos == previous.c_pos + 1 and item.pos == previous.pos + strand:
            blocks[-1].append(item)
        else:
            blocks.append([item])
    return strand, blocks


def _fetch_sequence(assembly: str, chrom: str, start_1: int, end_1: int) -> str:
    genome = "hg19" if assembly == "GRCh37" else "hg38"
    query = urllib.parse.urlencode({"genome": genome, "chrom": f"chr{chrom}", "start": start_1 - 1, "end": end_1})
    url = f"https://api.genome.ucsc.edu/getData/sequence?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "ARIANE-intronic-cache-builder/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.loads(response.read())
    sequence = str(payload.get("dna", "")).upper()
    if len(sequence) != end_1 - start_1 + 1 or set(sequence) - set("ACGT"):
        raise ValueError(f"Invalid reference sequence from {url}")
    return sequence


def _reference_lookup(assembly: str, chrom: str, positions: set[int]) -> dict[int, str]:
    start, end = min(positions), max(positions)
    sequence = _fetch_sequence(assembly, chrom, start, end)
    return {pos: sequence[pos - start] for pos in positions}


def build_coordinates(window: int) -> tuple[dict, dict]:
    snapshot = _read_json(SNAPSHOT)
    assembly_maps: dict[str, dict[tuple[str, int], tuple[str, int, str, int]]] = {}
    summary: dict[str, dict] = {}
    for assembly, field in (("GRCh37", "grch37"), ("GRCh38", "grch38")):
        assembly_maps[assembly] = {}
        for gene in TRANSCRIPTS:
            strand, blocks = exon_blocks(coding_positions(snapshot, gene, field))
            requested: list[tuple[str, int, str, int]] = []
            for left, right in zip(blocks, blocks[1:]):
                for offset in range(1, window + 1):
                    requested.append((f"c.{left[-1].c_pos}+{offset}", left[-1].pos + strand * offset, left[-1].chrom, strand))
                    requested.append((f"c.{right[0].c_pos}-{offset}", right[0].pos - strand * offset, right[0].chrom, strand))
            chrom = requested[0][2]
            refs = _reference_lookup(assembly, chrom, {item[1] for item in requested})
            for notation, pos, chrom, strand in requested:
                assembly_maps[assembly][(gene, notation)] = (chrom, pos, refs[pos], strand)
            summary.setdefault(gene, {})[assembly] = {"strand": strand, "coding_exon_blocks": len(blocks)}

    output: dict[str, dict] = {}
    for gene, notation in sorted(assembly_maps["GRCh38"]):
        chrom38, pos38, genomic_ref38, strand = assembly_maps["GRCh38"][(gene, notation)]
        chrom37, pos37, genomic_ref37, strand37 = assembly_maps["GRCh37"][(gene, notation)]
        if strand != strand37:
            raise ValueError(f"Strand mismatch for {gene}")
        transcript_ref = genomic_ref38 if strand == 1 else genomic_ref38.translate(COMPLEMENT)
        transcript_ref37 = genomic_ref37 if strand == 1 else genomic_ref37.translate(COMPLEMENT)
        if transcript_ref != transcript_ref37:
            raise ValueError(f"Reference-build mismatch for {gene}:{notation}")
        for transcript_alt in "ACGT":
            if transcript_alt == transcript_ref:
                continue
            genomic_alt38 = transcript_alt if strand == 1 else transcript_alt.translate(COMPLEMENT)
            genomic_alt37 = transcript_alt if strand == 1 else transcript_alt.translate(COMPLEMENT)
            c_notation = f"{notation}{transcript_ref}>{transcript_alt}"
            output[f"{gene}:{c_notation}"] = {
                "gene": gene,
                "transcript": TRANSCRIPTS[gene]["refseq"],
                "c_notation": c_notation,
                "status": "ok",
                "source": "versioned_intronic_coordinate_map",
                "grch37": {"chrom": chrom37, "pos": pos37, "ref": genomic_ref37, "alt": genomic_alt37, "assembly": "GRCh37"},
                "grch38": {"chrom": chrom38, "pos": pos38, "ref": genomic_ref38, "alt": genomic_alt38, "assembly": "GRCh38"},
            }
    _atomic_json(COORDINATES, output)
    metadata = {
        "dataset": "BRCA1/BRCA2 intronic SNV coordinate map",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "window_bp_per_exon_boundary": window,
        "variants": len(output),
        "source_snapshot": str(SNAPSHOT.relative_to(ROOT)).replace("\\", "/"),
        "source_snapshot_sha256": _sha256(SNAPSHOT),
        "reference_sequence_source": "UCSC Genome Browser sequence API (hg19/hg38)",
        "reference_transcripts": TRANSCRIPTS,
        "genes": summary,
        "sha256": _sha256(COORDINATES),
    }
    _atomic_json(COORDINATE_METADATA, metadata)
    return output, metadata


def _score_variant(api_url: str, gene: str, entry: dict, timeout: int, distance: int) -> dict:
    grch38 = entry["grch38"]
    variant = f"chr{grch38['chrom']}-{grch38['pos']}-{grch38['ref']}-{grch38['alt']}"
    url = api_url.rstrip("/") + "/?" + urllib.parse.urlencode({"hg": 38, "variant": variant, "distance": distance, "mask": 0})
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ARIANE-intronic-cache-builder/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read())
    rows = payload.get("scores") or []
    reference = TRANSCRIPTS[gene]
    matching = []
    for row in rows:
        tid = str(row.get("t_id") or "")
        refseq_ids = [str(value) for value in row.get("t_refseq_ids") or []]
        if tid.split(".")[0] == reference["ensembl"].split(".")[0] or any(value.split(".")[0] == reference["refseq"].split(".")[0] for value in refseq_ids):
            matching.append(row)
    if not matching:
        raise ValueError("Reference transcript absent from SpliceAI response")
    row = matching[-1]
    fields = {name: float(row.get(name, 0) or 0) for name in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")}
    max_field = max(fields, key=fields.get)
    return {
        "status": "ok", "score": fields[max_field], "max_delta_field": max_field,
        "delta_scores": fields, "selected_transcript": str(row.get("t_id") or ""),
        "selected_transcript_policy": "reference_transcript", "source": "Broad-compatible SpliceAI API",
        "variant": variant, "grch38": f"{grch38['chrom']}:{grch38['pos']}:{grch38['ref']}>{grch38['alt']}",
    }


def build_spliceai(api_url: str, workers: int, timeout: int, delay: float, distance: int) -> tuple[dict, dict]:
    coordinates = _read_json(COORDINATES)
    cache = _read_json(SPLICEAI) if SPLICEAI.exists() else {}
    pending = [(key, value) for key, value in coordinates.items() if cache.get(key, {}).get("status") != "ok"]
    lock = threading.Lock()
    api_urls = [value.strip() for value in api_url.split(",") if value.strip()]
    if not api_urls:
        raise ValueError("At least one SpliceAI API URL is required")

    def task(item: tuple[str, dict]) -> tuple[str, dict]:
        key, entry = item
        if delay:
            time.sleep(delay)
        endpoint_index = int(hashlib.sha256(key.encode("utf-8")).hexdigest(), 16) % len(api_urls)
        return key, _score_variant(api_urls[endpoint_index], entry["gene"], entry, timeout, distance)

    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(task, item): item[0] for item in pending}
        for future in concurrent.futures.as_completed(futures):
            key = futures[future]
            try:
                result_key, result = future.result()
                cache[result_key] = result
            except Exception as exc:
                cache[key] = {"status": "error", "error": f"{type(exc).__name__}: {exc}"}
            completed += 1
            if completed % 100 == 0 or completed == len(pending):
                with lock:
                    _atomic_json(SPLICEAI, cache)
                print(f"SpliceAI: {completed}/{len(pending)} attempted; {sum(v.get('status') == 'ok' for v in cache.values())}/{len(coordinates)} complete", flush=True)
    _atomic_json(SPLICEAI, cache)
    ok = sum(value.get("status") == "ok" for value in cache.values())
    metadata = {
        "dataset": "BRCA1/BRCA2 intronic SNV reference-transcript SpliceAI cache",
        "created_utc": datetime.now(timezone.utc).isoformat(), "coordinate_cache_sha256": _sha256(COORDINATES),
        "coordinate_variants": len(coordinates), "cache_entries": len(cache), "status_ok": ok,
        "status_error": len(cache) - ok, "api_urls": api_urls, "distance": distance,
        "mask": 0, "transcript_policy": "reference_transcript",
        "reference_transcripts": TRANSCRIPTS, "sha256": _sha256(SPLICEAI),
    }
    _atomic_json(SPLICEAI_METADATA, metadata)
    return cache, metadata


def validate() -> None:
    coordinates = _read_json(COORDINATES)
    coordinate_metadata = _read_json(COORDINATE_METADATA)
    if coordinate_metadata["sha256"] != _sha256(COORDINATES):
        raise SystemExit("Coordinate cache checksum mismatch")
    if coordinate_metadata["variants"] != len(coordinates):
        raise SystemExit("Coordinate cache count mismatch")
    if SPLICEAI.exists():
        spliceai = _read_json(SPLICEAI)
        missing = sorted(set(coordinates) - {key for key, value in spliceai.items() if value.get("status") == "ok"})
        if missing:
            raise SystemExit(f"SpliceAI cache incomplete: {len(missing)} variants missing/error")
        print(f"Validated {len(coordinates)} coordinates and {len(spliceai)} SpliceAI entries")
    else:
        print(f"Validated {len(coordinates)} coordinates; SpliceAI cache not built")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("coordinates", "spliceai", "all", "validate"))
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument("--api-url", default=DEFAULT_API)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=90)
    parser.add_argument("--distance", type=int, default=50)
    parser.add_argument("--delay", type=float, default=1.5, help="Per-request delay; use 0 only for a local endpoint")
    args = parser.parse_args()
    if args.command in {"coordinates", "all"}:
        _, metadata = build_coordinates(args.window)
        print(json.dumps(metadata, indent=2))
    if args.command in {"spliceai", "all"}:
        _, metadata = build_spliceai(args.api_url, args.workers, args.timeout, args.delay, args.distance)
        print(json.dumps(metadata, indent=2))
    if args.command == "validate":
        validate()


if __name__ == "__main__":
    main()
