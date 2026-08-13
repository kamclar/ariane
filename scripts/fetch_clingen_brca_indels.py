"""Fetch a versioned BRCA1/2 indel overlay from ClinGen Allele Registry.

This script deliberately writes a separate overlay. It never modifies the
existing BRCA Exchange indel snapshot. Only records with an exact configured
RefSeq transcript, a protein consequence, and GRCh37 plus GRCh38 mappings are
retained. Runtime integration is a separate, explicitly validated step.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/precomputed/brca_clingen_indel_overlay.index.json"
DEFAULT_METADATA = ROOT / "data/precomputed/brca_clingen_indel_overlay.metadata.json"
DEFAULT_CACHE_DIR = ROOT / "data/sources/clingen_allele_registry/chunks"
API = "https://reg.genome.network/alleles"
TRANSCRIPTS = {"BRCA1": "NM_007294.4", "BRCA2": "NM_000059.4"}
# Querying beyond the transcript end is harmless and avoids hard-coding a
# sequence length as biological data. The registry returns only overlaps.
QUERY_END = 12_000
INDEL_RE = re.compile(r"delins|del|dup|ins", re.IGNORECASE)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _request(transcript: str, begin: int, end: int, timeout: int) -> list[dict]:
    params = urllib.parse.urlencode(
        {
            "refseq": transcript,
            "begin": begin,
            "end": end,
            "fields": "none @id genomicAlleles transcriptAlleles",
        }
    )
    request = urllib.request.Request(
        f"{API}?{params}",
        headers={"Accept": "application/json", "User-Agent": "ARIANE-indel-map-builder/1.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        value = json.load(response)
    if not isinstance(value, list):
        raise RuntimeError(f"ClinGen returned {type(value).__name__}, expected a list")
    return value


def _exact_transcript_allele(record: dict, transcript: str) -> dict | None:
    prefix = transcript + ":"
    for allele in record.get("transcriptAlleles") or []:
        if any(str(item).startswith(prefix) for item in allele.get("hgvs") or []):
            return allele
    return None


def _c_notation(allele: dict, transcript: str) -> str:
    prefix = transcript + ":"
    values = [str(item)[len(prefix):] for item in allele.get("hgvs") or [] if str(item).startswith(prefix)]
    values = sorted(set(item for item in values if item.startswith("c.") and INDEL_RE.search(item)))
    return values[0] if len(values) == 1 else ""


def _p_notation(allele: dict) -> str:
    value = str((allele.get("proteinEffect") or {}).get("hgvs") or "")
    if ":p." not in value:
        return ""
    protein = value.split(":", 1)[1]
    if not protein.startswith("p.("):
        protein = f"p.({protein[2:]})"
    return re.sub(r"\*(?=\d|\))", "Ter", protein)


def _genomic(record: dict, assembly: str) -> dict | None:
    matches = [item for item in record.get("genomicAlleles") or [] if item.get("referenceGenome") == assembly]
    if len(matches) != 1:
        return None
    allele = matches[0]
    coordinates = allele.get("coordinates") or []
    if len(coordinates) != 1:
        return None
    coordinate = coordinates[0]
    ref = str(coordinate.get("referenceAllele") or "").upper()
    alt = str(coordinate.get("allele") or "").upper()
    if (ref and not re.fullmatch(r"[ACGT]+", ref)) or (alt and not re.fullmatch(r"[ACGT]+", alt)):
        return None
    return {
        "assembly": assembly,
        "chrom": str(allele.get("chromosome") or ""),
        "start_0based": int(coordinate["start"]),
        "end_0based": int(coordinate["end"]),
        "ref": ref,
        "alt": alt,
        "hgvs": next((str(item) for item in allele.get("hgvs") or [] if item.startswith("NC_")), ""),
    }


def _aliases(c_notation: str, transcript_allele: dict) -> list[str]:
    found = {c_notation}
    coordinates = transcript_allele.get("coordinates") or []
    if len(coordinates) == 1:
        ref = str(coordinates[0].get("referenceAllele") or "").upper()
        alt = str(coordinates[0].get("allele") or "").upper()
        if ref and not alt and c_notation.endswith("del"):
            found.add(c_notation + ref)
        if not ref and alt and c_notation.endswith(("dup", "ins")):
            found.add(c_notation + alt)
    return sorted(found)


def _kind(c_notation: str, p_notation: str) -> str:
    c, p = c_notation.lower(), p_notation.lower()
    if "fs" in p:
        return "frameshift"
    if "delins" in c:
        return "inframe_delins"
    if "dup" in c:
        return "inframe_duplication"
    if "ins" in c:
        return "inframe_insertion"
    if "del" in c:
        return "inframe_deletion"
    return ""


def build(
    output: Path,
    metadata_path: Path,
    *,
    query_start: int,
    query_end: int,
    chunk_size: int,
    timeout: int,
    retries: int,
    workers: int,
    cache_dir: Path,
) -> dict:
    records_by_caid: dict[str, dict] = {}
    stats: Counter[str] = Counter()
    request_count = 0
    cache_hits = 0
    cache_dir.mkdir(parents=True, exist_ok=True)

    jobs = [
        (gene, transcript, begin, min(begin + chunk_size, query_end))
        for gene, transcript in TRANSCRIPTS.items()
        for begin in range(query_start, query_end, chunk_size)
    ]

    def fetch(job: tuple[str, str, int, int]) -> tuple[str, str, list[dict], bool]:
        gene, transcript, begin, end = job
        cache_path = cache_dir / f"{transcript}_{begin}_{end}.json"
        if cache_path.is_file():
            with cache_path.open(encoding="utf-8") as handle:
                cached = json.load(handle)
            if not isinstance(cached, list):
                raise RuntimeError(f"Invalid cached ClinGen response: {cache_path}")
            return gene, transcript, cached, True
        for attempt in range(retries + 1):
            try:
                batch = _request(transcript, begin, end, timeout)
                temporary = cache_path.with_suffix(".json.tmp")
                temporary.write_bytes((json.dumps(batch, separators=(",", ":")) + "\n").encode())
                temporary.replace(cache_path)
                return gene, transcript, batch, False
            except Exception:
                if attempt >= retries:
                    raise
                time.sleep(min(2 ** attempt, 10))
        raise AssertionError("unreachable")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        for completed, (gene, transcript, batch, was_cached) in enumerate(executor.map(fetch, jobs), start=1):
            request_count += 1
            cache_hits += int(was_cached)
            print(
                f"ClinGen chunk {completed}/{len(jobs)}: {transcript}; "
                f"records={len(batch)}; source={'cache' if was_cached else 'api'}",
                flush=True,
            )
            stats["api_records_seen"] += len(batch)
            for raw in batch:
                caid = str(raw.get("@id") or "").rsplit("/", 1)[-1]
                if not caid.startswith("CA"):
                    stats["missing_caid"] += 1
                    continue
                allele = _exact_transcript_allele(raw, transcript)
                if not allele:
                    stats["missing_exact_transcript"] += 1
                    continue
                c_notation = _c_notation(allele, transcript)
                if not c_notation:
                    stats["not_single_coding_indel"] += 1
                    continue
                p_notation = _p_notation(allele)
                if not p_notation:
                    stats["missing_protein"] += 1
                    continue
                grch37 = _genomic(raw, "GRCh37")
                grch38 = _genomic(raw, "GRCh38")
                if not grch37 or not grch38:
                    stats["missing_genomic_mapping"] += 1
                    continue
                record = {
                    "gene": gene,
                    "reference_transcript": transcript,
                    "canonical_c_notation": c_notation,
                    "input_c_notations": _aliases(c_notation, allele),
                    "p_notation": p_notation,
                    "variant_type": _kind(c_notation, p_notation),
                    "grch37_unanchored": grch37,
                    "grch38_unanchored": grch38,
                    "source": {"dataset": "ClinGen Allele Registry", "ca_id": caid, "api": API},
                }
                previous = records_by_caid.get(caid)
                if previous and previous != record:
                    raise RuntimeError(f"ClinGen CAID changed across overlapping queries: {caid}")
                records_by_caid[caid] = record

    by_key: dict[str, dict] = {}
    conflicts: dict[str, list[str]] = {}
    for caid, record in records_by_caid.items():
        key = f"{record['gene']}:{record['canonical_c_notation']}"
        previous = by_key.get(key)
        if previous and previous != record:
            conflicts.setdefault(key, [previous["source"]["ca_id"]]).append(caid)
            continue
        by_key[key] = record
    for key in conflicts:
        by_key.pop(key, None)

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes((json.dumps(dict(sorted(by_key.items())), indent=2, sort_keys=True) + "\n").encode())
    metadata = {
        "dataset": "ClinGen Allele Registry BRCA1/2 coding-indel overlay",
        "status": "comparison_only_not_runtime",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": API,
        "query_transcripts": TRANSCRIPTS,
        "query_start": query_start,
        "query_end": query_end,
        "query_chunk_size": chunk_size,
        "parallel_workers": workers,
        "chunks_processed": request_count,
        "api_requests": request_count - cache_hits,
        "cache_hits": cache_hits,
        "source_chunks": {
            path.name: sha256(path)
            for path in sorted(cache_dir.glob("*.json"))
            if re.fullmatch(r"NM_\d+\.\d+_\d+_\d+\.json", path.name)
        },
        "unique_caids": len(records_by_caid),
        "records": len(by_key),
        "conflicts": conflicts,
        "excluded": dict(sorted(stats.items())),
        "index_sha256": sha256(output),
    }
    metadata_path.write_bytes((json.dumps(metadata, indent=2, sort_keys=True) + "\n").encode())
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--query-start", type=int, default=0)
    parser.add_argument("--query-end", type=int, default=QUERY_END)
    parser.add_argument("--chunk-size", type=int, default=250)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    args = parser.parse_args()
    print(json.dumps(build(
        args.output,
        args.metadata,
        query_start=args.query_start,
        query_end=args.query_end,
        chunk_size=args.chunk_size,
        timeout=args.timeout,
        retries=args.retries,
        workers=args.workers,
        cache_dir=args.cache_dir,
    ), indent=2))


if __name__ == "__main__":
    main()
