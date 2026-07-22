"""Build the authoritative gnomAD v3 BRCA frequency and coverage snapshot.

Requires requests and pysam. The script reads the official gnomAD v3.1.2
regional VCF through its tabix index and obtains per-base genome coverage from
the official gnomAD browser API. Existing v2.1 records are retained.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
GNOMAD_DIR = ROOT / "backend" / "data" / "gnomad"
VARIANT_PATH = GNOMAD_DIR / "gnomad_brca_region_cache_by_variant.with_real_coverage.json"
COVERAGE_PATH = GNOMAD_DIR / "gnomad_brca_coverage_cache.json"
COVERAGE_PROGRESS_PATH = GNOMAD_DIR / ".gnomad_v3_coverage_build_progress.json"
DATASET = "gnomad_v3_1_2_genomes_grch38"
VCF_TEMPLATE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/"
    "vcf/genomes/gnomad.genomes.v3.1.2.sites.chr{chrom}.vcf.bgz"
)
API_URL = "https://gnomad.broadinstitute.org/api"
COVERAGE_QUERY = """
query Coverage($chrom: String!, $start: Int!, $stop: Int!, $dataset: DatasetId!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    coverage(dataset: $dataset) {
      genome { pos mean median over_20 over_25 }
    }
  }
}
"""

# GRCh38 gene intervals, padded to match the v2 regional snapshot policy.
REGIONS = {
    "BRCA1": {"chrom": "17", "start": 43044295, "end": 43125364},
    "BRCA2": {"chrom": "13", "start": 32315474, "end": 32400266},
}
PADDING = 1000
NON_FOUNDER_ANCESTRIES = ("afr", "amr", "eas", "nfe", "sas")


def _atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, separators=(",", ":"))
    os.replace(temporary, path)


def _alt_value(value, index: int):
    if isinstance(value, tuple):
        return value[index] if index < len(value) else None
    return value


def _extract_variants() -> tuple[dict, list]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - build-time dependency
        raise SystemExit("pysam is required to build the snapshot") from exc
    variants = {}
    extraction_log = []
    for gene, region in REGIONS.items():
        chrom = region["chrom"]
        start = region["start"] - PADDING
        end = region["end"] + PADDING
        url = VCF_TEMPLATE.format(chrom=chrom)
        source = pysam.VariantFile(url)
        count = 0
        for record in source.fetch(f"chr{chrom}", start - 1, end):
            for alt_index, alt in enumerate(record.alts or ()):
                ac = _alt_value(record.info.get("AC_non_cancer"), alt_index)
                an = record.info.get("AN_non_cancer")
                af = _alt_value(record.info.get("AF_non_cancer"), alt_index)
                if not ac or not an:
                    continue
                ancestry_values = []
                for ancestry in NON_FOUNDER_ANCESTRIES:
                    value = _alt_value(record.info.get(f"AF_non_cancer_{ancestry}"), alt_index)
                    if value is not None:
                        ancestry_values.append((float(value), ancestry))
                popmax_af, popmax_pop = max(ancestry_values, default=(None, None))
                variant_id = f"{chrom}-{record.pos}-{record.ref}-{alt}"
                variants.setdefault(variant_id, []).append({
                    "variant_id": variant_id,
                    "chrom": chrom,
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": alt,
                    "filter": ";".join(record.filter.keys()) or "PASS",
                    "af": float(af) if af is not None else None,
                    "ac": int(ac),
                    "an": int(an),
                    "nhomalt": int(_alt_value(record.info.get("nhomalt_non_cancer"), alt_index) or 0),
                    "faf95_max": None,
                    "popmax_af": popmax_af,
                    "popmax_pop": popmax_pop,
                    "dataset": DATASET,
                    "build": "GRCh38",
                })
                count += 1
        source.close()
        extraction_log.append({
            "dataset": DATASET, "status": "ok", "chrom": chrom,
            "gene": gene, "start": start, "end": end, "records": count,
            "source": url,
        })
    return variants, extraction_log


def _coverage_chunk(
    session: requests.Session, chrom: str, start: int, stop: int, attempts: int = 5
) -> list:
    payload = {
        "query": COVERAGE_QUERY,
        "variables": {"chrom": chrom, "start": start, "stop": stop, "dataset": "gnomad_r3"},
    }
    last_error = None
    for attempt in range(attempts):
        try:
            response = session.post(API_URL, json=payload, timeout=90)
            if response.status_code == 429:
                wait_seconds = int(response.headers.get("Retry-After") or 60)
                time.sleep(max(wait_seconds, 30))
                continue
            if not response.ok:
                raise RuntimeError(f"HTTP {response.status_code}: {response.text[:2000]}")
            body = response.json()
            if body.get("errors"):
                raise RuntimeError(body["errors"])
            return body["data"]["region"]["coverage"]["genome"]
        except Exception as exc:
            last_error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"gnomAD coverage query failed for {chrom}:{start}-{stop}: {last_error}")


def _extract_coverage(chunk_size: int = 500) -> dict:
    if COVERAGE_PROGRESS_PATH.exists():
        progress = json.loads(COVERAGE_PROGRESS_PATH.read_text(encoding="utf-8"))
    else:
        progress = {"coverage": {}, "completed": []}
    result = progress.get("coverage", {})
    completed = set(progress.get("completed", []))
    session = requests.Session()
    session.headers.update({"User-Agent": "ARIANE gnomAD snapshot builder/1.0"})
    for gene, region in REGIONS.items():
        chrom = region["chrom"]
        interval_start = region["start"] - PADDING
        interval_end = region["end"] + PADDING
        for start in range(interval_start, interval_end + 1, chunk_size):
            stop = min(start + chunk_size - 1, interval_end)
            chunk_key = f"{chrom}:{start}-{stop}"
            if chunk_key in completed:
                continue
            for item in _coverage_chunk(session, chrom, start, stop, attempts=10):
                pos = int(item["pos"])
                key = f"{DATASET}|GRCh38|{chrom}|{pos}"
                mean = item.get("mean")
                result[key] = {
                    "chrom": chrom, "pos": pos,
                    "mean_depth": mean, "median_depth": item.get("median"),
                    "over_20": item.get("over_20"), "over_25": item.get("over_25"),
                    "dataset": DATASET, "dataset_key": DATASET,
                    "build": "GRCh38", "gene": gene, "threshold": 25.0,
                    "passes": mean is not None and float(mean) >= 25.0,
                    "source": "gnomAD browser coverage API (gnomad_r3)",
                    "position_key": key,
                }
            completed.add(chunk_key)
            _atomic_json(
                COVERAGE_PROGRESS_PATH,
                {"coverage": result, "completed": sorted(completed)},
            )
            time.sleep(0.5)
    return result


def build() -> None:
    variant_payload = json.loads(VARIANT_PATH.read_text(encoding="utf-8"))
    coverage_payload = json.loads(COVERAGE_PATH.read_text(encoding="utf-8"))
    new_variants, extraction_log = _extract_variants()
    new_coverage = _extract_coverage()

    old_variants = variant_payload.get("variants", {})
    merged_variants = {}
    for key, records in old_variants.items():
        kept = [record for record in records if record.get("dataset") != DATASET]
        if kept:
            merged_variants[key] = kept
    for key, records in new_variants.items():
        merged_variants.setdefault(key, []).extend(records)

    old_coverage = coverage_payload.get("coverage_by_position", {})
    merged_coverage = {
        key: value for key, value in old_coverage.items()
        if value.get("dataset_key") != DATASET
    }
    merged_coverage.update(new_coverage)

    metadata = variant_payload.setdefault("metadata", {})
    metadata.setdefault("regions", {})["GRCh38"] = REGIONS
    old_log = [item for item in metadata.get("extraction_log", []) if item.get("dataset") != DATASET]
    metadata["extraction_log"] = old_log + extraction_log
    metadata["region_padding_bp"] = PADDING
    metadata["updated_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["v3_source"] = VCF_TEMPLATE
    metadata["v3_records"] = sum(len(records) for records in new_variants.values())
    variant_payload["variants"] = merged_variants

    coverage_metadata = coverage_payload.setdefault("metadata", {})
    coverage_metadata.update({
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "records": len(merged_coverage),
        "v3_records": len(new_coverage),
        "v3_source": API_URL,
        "v3_dataset": "gnomad_r3 genome coverage",
    })
    coverage_payload["coverage_by_position"] = merged_coverage

    _atomic_json(VARIANT_PATH, variant_payload)
    _atomic_json(COVERAGE_PATH, coverage_payload)
    COVERAGE_PROGRESS_PATH.unlink(missing_ok=True)
    for path in (VARIANT_PATH, COVERAGE_PATH):
        print(path, hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_size)
    print(f"v3 variants={metadata['v3_records']} coverage_positions={len(new_coverage)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe-coverage", action="store_true")
    args = parser.parse_args()
    if args.probe_coverage:
        values = _coverage_chunk(requests.Session(), "17", 43043295, 43044294, attempts=1)
        print(len(values), values[0], values[-1])
    else:
        build()
