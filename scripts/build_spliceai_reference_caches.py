"""Build ENIGMA Appendix J compatible BRCA SpliceAI reference caches.

Only local, pinned Broad-compatible SpliceAI servers are accepted. Builds use
separate checkpoint files and replace an active cache only after every expected
variant has a complete, validated result.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import re
import sys
import threading
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from backend.spliceai_profile import SPLICEAI_PROFILE, scoring_profile_metadata
from scripts.build_intronic_reference_cache import (
    ROOT,
    SPLICEAI as INTRONIC_CACHE,
    SPLICEAI_METADATA as INTRONIC_METADATA,
    TRANSCRIPTS,
    _atomic_json,
    _read_json,
    _score_variant,
    _sha256,
    build_spliceai as build_intronic_spliceai,
)


CODING_SOURCE = ROOT / "data/precomputed/brca_module1_snv_classification_snapshot.index.json"
CODING_CACHE = ROOT / "data/spliceai/spliceai_brca_snv_reference_cache.json"
CODING_METADATA = ROOT / "data/spliceai/spliceai_brca_snv_reference_cache.metadata.json"
CODING_WORK = ROOT / "data/spliceai/build/spliceai_brca_snv_reference_cache.building.json"
CODING_WORK_METADATA = CODING_WORK.with_name(CODING_WORK.stem + ".metadata.json")
DEFAULT_API = "http://127.0.0.1:8080/spliceai/"


def _parse_grch38(value: str) -> dict:
    match = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([ACGT]+)>([ACGT]+)", value or "")
    if not match:
        raise ValueError(f"Invalid GRCh38 coordinate: {value!r}")
    return {
        "chrom": match.group(1),
        "pos": int(match.group(2)),
        "ref": match.group(3),
        "alt": match.group(4),
        "assembly": "GRCh38",
    }


def coding_source_records() -> dict[str, dict]:
    snapshot = _read_json(CODING_SOURCE)
    records = {}
    for key, entry in snapshot.items():
        match = re.fullmatch(r"(BRCA[12]):(c\.\d+[ACGT]>[ACGT])", key)
        if not match:
            raise ValueError(f"Unexpected coding SNV key: {key}")
        records[key] = {
            "gene": match.group(1),
            "c_notation": match.group(2),
            "grch38": _parse_grch38(entry.get("grch38", "")),
        }
    if len(records) != 47547:
        raise ValueError(f"Expected 47,547 coding SNVs, found {len(records):,}")
    return records


def _local_api_urls(value: str) -> list[str]:
    urls = [item.strip() for item in value.split(",") if item.strip()]
    if not urls:
        raise ValueError("At least one local SpliceAI API URL is required")
    for url in urls:
        if urllib.parse.urlparse(url).hostname not in {"127.0.0.1", "localhost", "::1"}:
            raise ValueError(f"Remote SpliceAI endpoint rejected for reference build: {url}")
    return urls


def _checkpoint_metadata(records: dict, cache: dict, api_urls: list[str]) -> dict:
    ok = sum(entry.get("status") == "ok" for entry in cache.values())
    metadata = {
        **scoring_profile_metadata(),
        "build_status": "building",
        "source_snapshot_sha256": _sha256(CODING_SOURCE),
        "expected_variants": len(records),
        "cache_entries": len(cache),
        "status_ok": ok,
        "status_error": len(cache) - ok,
        "api_urls": api_urls,
    }
    if CODING_WORK.exists():
        metadata["sha256"] = _sha256(CODING_WORK)
    return metadata


def _load_work(records: dict, api_urls: list[str]) -> dict:
    if not CODING_WORK.exists() and not CODING_WORK_METADATA.exists():
        return {}
    if not CODING_WORK.exists() or not CODING_WORK_METADATA.exists():
        raise RuntimeError(
            f"Incomplete build checkpoint at {CODING_WORK.parent}; preserve or move it aside before a new build"
        )
    metadata = _read_json(CODING_WORK_METADATA)
    expected = _checkpoint_metadata(records, {}, api_urls)
    protected_fields = (
        "scoring_profile_id", "scoring_profile_sha256", "genome_assembly",
        "distance", "mask", "annotation_subset", "transcript_policy",
        "aggregation", "delta_score_fields", "reference_score_fields",
        "alternate_score_fields", "source_snapshot_sha256", "expected_variants",
    )
    mismatches = [field for field in protected_fields if metadata.get(field) != expected.get(field)]
    if metadata.get("sha256") != _sha256(CODING_WORK):
        mismatches.append("sha256")
    if mismatches:
        raise RuntimeError(
            "Coding SpliceAI checkpoint is incompatible: " + ", ".join(mismatches)
        )
    return _read_json(CODING_WORK)


def _write_checkpoint(records: dict, cache: dict, api_urls: list[str]) -> None:
    _atomic_json(CODING_WORK, cache)
    _atomic_json(CODING_WORK_METADATA, _checkpoint_metadata(records, cache, api_urls))


def build_coding(api_url: str, workers: int, timeout: int, delay: float) -> tuple[dict, dict]:
    records = coding_source_records()
    api_urls = _local_api_urls(api_url)
    cache = _load_work(records, api_urls)
    pending = [(key, entry) for key, entry in records.items() if cache.get(key, {}).get("status") != "ok"]
    lock = threading.Lock()

    def task(index: int, item: tuple[str, dict]) -> tuple[str, dict]:
        key, entry = item
        if delay:
            time.sleep(delay)
        endpoint = api_urls[index % len(api_urls)]
        return key, _score_variant(endpoint, entry["gene"], entry, timeout)

    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(task, index, item): item[0]
            for index, item in enumerate(pending)
        }
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
                    _write_checkpoint(records, cache, api_urls)
                ok = sum(value.get("status") == "ok" for value in cache.values())
                print(f"Coding SpliceAI: {completed}/{len(pending)} attempted; {ok}/{len(records)} complete", flush=True)

    _write_checkpoint(records, cache, api_urls)
    ok = sum(value.get("status") == "ok" for value in cache.values())
    if ok != len(records) or len(cache) != len(records):
        raise RuntimeError(
            f"Coding SpliceAI build incomplete: {ok}/{len(records)} successful. Production cache was not replaced."
        )
    _atomic_json(CODING_CACHE, cache)
    metadata = {
        **scoring_profile_metadata(),
        "dataset": "BRCA1/BRCA2 coding SNV reference-transcript SpliceAI cache",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_snapshot": str(CODING_SOURCE.relative_to(ROOT)).replace("\\", "/"),
        "source_snapshot_sha256": _sha256(CODING_SOURCE),
        "expected_variants": len(records),
        "cache_entries": len(cache),
        "status_ok": ok,
        "status_error": 0,
        "api_urls": api_urls,
        "reference_transcripts": TRANSCRIPTS,
        "approved_engine": SPLICEAI_PROFILE["approved_engine"],
        "sha256": _sha256(CODING_CACHE),
    }
    _atomic_json(CODING_METADATA, metadata)
    return cache, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", choices=("coding", "intronic", "all"))
    parser.add_argument("--api-url", default=DEFAULT_API, help="Comma-separated local pinned server URLs")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--delay", type=float, default=0.0)
    args = parser.parse_args()
    if args.dataset in {"coding", "all"}:
        _, metadata = build_coding(args.api_url, args.workers, args.timeout, args.delay)
        print(json.dumps(metadata, indent=2))
    if args.dataset in {"intronic", "all"}:
        _, metadata = build_intronic_spliceai(args.api_url, args.workers, args.timeout, args.delay)
        print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
