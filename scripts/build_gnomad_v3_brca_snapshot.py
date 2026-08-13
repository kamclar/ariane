"""Build the authoritative gnomAD BRCA frequency and coverage snapshot.

Requires requests, pysam, and scipy. The script reads official gnomAD v2.1.1
exome and v3.1.2 genome regional VCFs through their tabix indexes and obtains
per-base genome coverage from the official gnomAD browser API.

gnomAD v3.1.2 publishes non-cancer AC/AN values but not non-cancer FAF95 in its
VCF. This builder reproduces Hail's official filtering_allele_frequency
algorithm for each non-founder ancestry and stores their maximum. BA1/BS1 must
never be calculated from raw AF or popmax AF as a fallback.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
from scipy.optimize import brentq
from scipy.stats import poisson

ROOT = Path(__file__).resolve().parents[1]
GNOMAD_DIR = ROOT / "backend" / "data" / "gnomad"
VARIANT_PATH = GNOMAD_DIR / "gnomad_brca_region_cache_by_variant.with_real_coverage.json"
COVERAGE_PATH = GNOMAD_DIR / "gnomad_brca_coverage_cache.json"
COVERAGE_PROGRESS_PATH = GNOMAD_DIR / ".gnomad_v3_coverage_build_progress.json"
DATASET = "gnomad_v3_1_2_genomes_grch38"  # Coverage cache compatibility.
V3_DATASET = DATASET
V2_DATASET = "gnomad_v2_1_1_exomes_grch37"
V3_VCF_TEMPLATE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/"
    "vcf/genomes/gnomad.genomes.v3.1.2.sites.chr{chrom}.vcf.bgz"
)
V2_VCF_TEMPLATE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/"
    "vcf/exomes/gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz"
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
V3_REGIONS = {
    "BRCA1": {"chrom": "17", "start": 43044295, "end": 43125364},
    "BRCA2": {"chrom": "13", "start": 32315474, "end": 32400266},
}
V2_REGIONS = {
    "BRCA1": {"chrom": "17", "start": 41196312, "end": 41277500},
    "BRCA2": {"chrom": "13", "start": 32889617, "end": 32973809},
}
PADDING = 1000
NON_FOUNDER_ANCESTRIES = ("afr", "amr", "eas", "nfe", "sas")
FAF_CONFIDENCE = 0.95
FAF_LOWER = 1e-10
FAF_UPPER = 2.0
FAF_TOLERANCE = 1e-7
FAF_PRECISION = 1e-6
HAIL_FAF_COMMIT = "27c50a725f36169fd24e2fcde5256a9cb37ad998"
HAIL_FAF_SOURCE = (
    f"https://github.com/hail-is/hail/blob/{HAIL_FAF_COMMIT}/hail/hail/src/is/hail/"
    "experimental/package.scala"
)
ENIGMA_FREQUENCY_THRESHOLDS = (0.00002, 0.0001, 0.001)


def _atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, separators=(",", ":"))
    os.replace(temporary, path)


def _alt_value(value, index: int):
    if isinstance(value, tuple):
        return value[index] if index < len(value) else None
    return value


def filtering_allele_frequency(
    ac: int,
    an: int,
    confidence: float = FAF_CONFIDENCE,
) -> float:
    """Reproduce Hail calcFilterAlleleFreq used by gnomAD.

    The control flow, bounds, tolerance, rounding, and singleton handling match
    the Hail implementation. scipy supplies the Poisson quantile and root
    solver; no statistical approximation is used here.
    """
    ac = int(ac)
    an = int(an)
    if ac <= 1 or an == 0:
        return 0.0

    def difference(allele_frequency: float) -> float:
        return float(ac - 1 - poisson.ppf(confidence, an * allele_frequency))

    root = brentq(
        difference,
        FAF_LOWER,
        FAF_UPPER,
        xtol=FAF_TOLERANCE,
    )
    rounder = 1.0 / (FAF_PRECISION / 100.0)
    max_af = int(root * rounder + 0.5) / rounder
    while int(poisson.ppf(confidence, an * max_af)) < ac:
        max_af += FAF_PRECISION
    return max(0.0, max_af - FAF_PRECISION)


def _non_cancer_faf95(record, alt_index: int) -> tuple[float | None, str | None, dict]:
    """Return maximum non-cancer FAF95 across ENIGMA non-founder ancestries."""
    ancestry_values = {}
    for ancestry in NON_FOUNDER_ANCESTRIES:
        ac = _alt_value(record.info.get(f"AC_non_cancer_{ancestry}"), alt_index)
        an = record.info.get(f"AN_non_cancer_{ancestry}")
        if ac is None or an in (None, 0):
            continue
        ancestry_values[ancestry] = filtering_allele_frequency(int(ac), int(an))
    if not ancestry_values:
        return None, None, {}
    population, value = max(ancestry_values.items(), key=lambda item: item[1])
    return value, population, ancestry_values


def _assert_faf95_matches_hail(
    record,
    alt_index: int,
    *,
    ac_prefix: str,
    an_prefix: str,
    faf_prefix: str,
) -> int:
    """Cross-check published VCF FAF95 values against the reproduced algorithm."""
    comparisons = 0
    for ancestry in NON_FOUNDER_ANCESTRIES:
        ac = _alt_value(record.info.get(f"{ac_prefix}{ancestry}"), alt_index)
        an = record.info.get(f"{an_prefix}{ancestry}")
        published = _alt_value(record.info.get(f"{faf_prefix}{ancestry}"), alt_index)
        if ac is None or an in (None, 0) or published is None:
            continue
        calculated = filtering_allele_frequency(int(ac), int(an))
        published_value = float(published)
        calculated_band = sum(
            calculated > threshold for threshold in ENIGMA_FREQUENCY_THRESHOLDS
        )
        published_band = sum(
            published_value > threshold for threshold in ENIGMA_FREQUENCY_THRESHOLDS
        )
        if calculated_band != published_band or not math.isclose(
            calculated,
            published_value,
            # scipy and the Hail/Breeze Poisson quantile differ slightly for
            # very large lambda; this tolerance is still far below 0.01%.
            rel_tol=2e-5,
            abs_tol=FAF_PRECISION * 1.01,
        ):
            raise RuntimeError(
                "Hail FAF95 cross-check failed for "
                f"{record.contig}:{record.pos}:{record.ref}:{alt_index} "
                f"{ancestry}: calculated={calculated}, published={published}"
            )
        comparisons += 1
    return comparisons


def _source_identity(url: str) -> dict:
    response = requests.head(url, timeout=90)
    response.raise_for_status()
    return {
        "url": url,
        "etag": response.headers.get("ETag"),
        "x_goog_hash": response.headers.get("x-goog-hash"),
        "content_length": response.headers.get("Content-Length"),
        "last_modified": response.headers.get("Last-Modified"),
    }


def _extract_v3_variants() -> tuple[dict, list]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - build-time dependency
        raise SystemExit("pysam is required to build the snapshot") from exc
    variants = {}
    extraction_log = []
    for gene, region in V3_REGIONS.items():
        chrom = region["chrom"]
        start = region["start"] - PADDING
        end = region["end"] + PADDING
        url = V3_VCF_TEMPLATE.format(chrom=chrom)
        source = pysam.VariantFile(url, index_filename=url + ".tbi")
        source_identity = _source_identity(url)
        count = 0
        faf_validation_count = 0
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
                faf95_max, faf95_pop, faf95_by_ancestry = _non_cancer_faf95(
                    record, alt_index
                )
                faf_validation_count += _assert_faf95_matches_hail(
                    record,
                    alt_index,
                    ac_prefix="AC_",
                    an_prefix="AN_",
                    faf_prefix="faf95_",
                )
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
                    "faf95_max": faf95_max,
                    "faf95_pop": faf95_pop,
                    "faf95_by_ancestry": faf95_by_ancestry,
                    "faf95_scope": "non_cancer_non_founder_ancestries",
                    "faf95_method": "hail_calcFilterAlleleFreq_poisson_95_ci",
                    "popmax_af": popmax_af,
                    "popmax_pop": popmax_pop,
                    "dataset": V3_DATASET,
                    "build": "GRCh38",
                })
                count += 1
        source.close()
        extraction_log.append({
            "dataset": V3_DATASET, "status": "ok", "chrom": chrom,
            "gene": gene, "start": start, "end": end, "records": count,
            "source": url, "source_identity": source_identity,
            "faf_algorithm_validation_comparisons": faf_validation_count,
        })
    return variants, extraction_log


def _extract_v2_variants() -> tuple[dict, list]:
    """Extract published non-cancer FAF95 fields from gnomAD v2.1.1 exomes."""
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - build-time dependency
        raise SystemExit("pysam is required to build the snapshot") from exc
    variants = {}
    extraction_log = []
    for gene, region in V2_REGIONS.items():
        chrom = region["chrom"]
        start = region["start"] - PADDING
        end = region["end"] + PADDING
        url = V2_VCF_TEMPLATE.format(chrom=chrom)
        source = pysam.VariantFile(url, index_filename=url + ".tbi")
        source_identity = _source_identity(url)
        count = 0
        faf_validation_count = 0
        for record in source.fetch(chrom, start - 1, end):
            for alt_index, alt in enumerate(record.alts or ()):
                ac = _alt_value(record.info.get("non_cancer_AC"), alt_index)
                an = record.info.get("non_cancer_AN")
                af = _alt_value(record.info.get("non_cancer_AF"), alt_index)
                if not ac or not an:
                    continue
                af_by_ancestry = {}
                faf95_by_ancestry = {}
                for ancestry in NON_FOUNDER_ANCESTRIES:
                    ancestry_af = _alt_value(
                        record.info.get(f"non_cancer_AF_{ancestry}"), alt_index
                    )
                    ancestry_faf95 = _alt_value(
                        record.info.get(f"non_cancer_faf95_{ancestry}"), alt_index
                    )
                    if ancestry_af is not None:
                        af_by_ancestry[ancestry] = float(ancestry_af)
                    if ancestry_faf95 is not None:
                        faf95_by_ancestry[ancestry] = float(ancestry_faf95)
                if not faf95_by_ancestry:
                    raise RuntimeError(
                        f"Published non-cancer FAF95 missing for {chrom}:{record.pos}:{alt}"
                    )
                faf_validation_count += _assert_faf95_matches_hail(
                    record,
                    alt_index,
                    ac_prefix="non_cancer_AC_",
                    an_prefix="non_cancer_AN_",
                    faf_prefix="non_cancer_faf95_",
                )
                popmax_pop, popmax_af = max(
                    af_by_ancestry.items(), key=lambda item: item[1]
                )
                faf95_pop, faf95_max = max(
                    faf95_by_ancestry.items(), key=lambda item: item[1]
                )
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
                    "nhomalt": int(
                        _alt_value(record.info.get("non_cancer_nhomalt"), alt_index) or 0
                    ),
                    "faf95_max": faf95_max,
                    "faf95_pop": faf95_pop,
                    "faf95_by_ancestry": faf95_by_ancestry,
                    "faf95_scope": "non_cancer_non_founder_ancestries",
                    "faf95_method": "official_vcf_non_cancer_faf95_poisson_95_ci",
                    "popmax_af": popmax_af,
                    "popmax_pop": popmax_pop,
                    "dataset": V2_DATASET,
                    "build": "GRCh37",
                })
                count += 1
        source.close()
        extraction_log.append({
            "dataset": V2_DATASET, "status": "ok", "chrom": chrom,
            "gene": gene, "start": start, "end": end, "records": count,
            "source": url, "source_identity": source_identity,
            "faf_algorithm_validation_comparisons": faf_validation_count,
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
    for gene, region in V3_REGIONS.items():
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


def build(*, refresh_coverage: bool = True) -> None:
    variant_payload = json.loads(VARIANT_PATH.read_text(encoding="utf-8"))
    coverage_payload = json.loads(COVERAGE_PATH.read_text(encoding="utf-8"))
    v2_variants, v2_extraction_log = _extract_v2_variants()
    v3_variants, v3_extraction_log = _extract_v3_variants()
    new_variants = {}
    for extracted in (v2_variants, v3_variants):
        for key, records in extracted.items():
            new_variants.setdefault(key, []).extend(records)
    extraction_log = v2_extraction_log + v3_extraction_log
    new_coverage = _extract_coverage() if refresh_coverage else None

    old_variants = variant_payload.get("variants", {})
    merged_variants = {}
    for key, records in old_variants.items():
        kept = [
            record for record in records
            if record.get("dataset") not in (V2_DATASET, V3_DATASET)
        ]
        if kept:
            merged_variants[key] = kept
    for key, records in new_variants.items():
        merged_variants.setdefault(key, []).extend(records)

    metadata = variant_payload.setdefault("metadata", {})
    metadata.setdefault("regions", {})["GRCh37"] = V2_REGIONS
    metadata.setdefault("regions", {})["GRCh38"] = V3_REGIONS
    old_log = [
        item for item in metadata.get("extraction_log", [])
        if item.get("dataset") not in (V2_DATASET, V3_DATASET)
    ]
    metadata["extraction_log"] = old_log + extraction_log
    metadata["region_padding_bp"] = PADDING
    metadata["updated_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["source"] = "official gnomAD v2.1.1 and v3.1.2 regional VCFs"
    metadata["v2_source"] = V2_VCF_TEMPLATE
    metadata["v3_source"] = V3_VCF_TEMPLATE
    metadata["v2_records"] = sum(len(records) for records in v2_variants.values())
    metadata["v3_records"] = sum(len(records) for records in v3_variants.values())
    metadata["variants_count"] = sum(len(records) for records in merged_variants.values())
    metadata["unique_ids"] = len(merged_variants)
    metadata["v2_faf95"] = {
        "scope": "non-cancer subset; maximum across afr, amr, eas, nfe, sas",
        "source_fields": "non_cancer_faf95_{ancestry}",
        "confidence": FAF_CONFIDENCE,
        "method": "published gnomAD v2.1.1 VCF annotation",
        "published_values_crosschecked_against_hail": sum(
            item.get("faf_algorithm_validation_comparisons", 0)
            for item in v2_extraction_log
        ),
        "raw_af_fallback_allowed": False,
    }
    metadata["v3_faf95"] = {
        "scope": "non-cancer subset; maximum across afr, amr, eas, nfe, sas",
        "source_fields": "AC_non_cancer_{ancestry}, AN_non_cancer_{ancestry}",
        "confidence": FAF_CONFIDENCE,
        "method": "Hail calcFilterAlleleFreq",
        "method_source": HAIL_FAF_SOURCE,
        "method_commit": HAIL_FAF_COMMIT,
        "published_all_samples_crosscheck": sum(
            item.get("faf_algorithm_validation_comparisons", 0)
            for item in v3_extraction_log
        ),
        "singletons": 0.0,
        "raw_af_fallback_allowed": False,
    }
    variant_payload["variants"] = merged_variants

    _atomic_json(VARIANT_PATH, variant_payload)
    output_paths = [VARIANT_PATH]
    if refresh_coverage:
        old_coverage = coverage_payload.get("coverage_by_position", {})
        merged_coverage = {
            key: value for key, value in old_coverage.items()
            if value.get("dataset_key") != DATASET
        }
        merged_coverage.update(new_coverage or {})
        coverage_metadata = coverage_payload.setdefault("metadata", {})
        coverage_metadata.update({
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "records": len(merged_coverage),
            "v3_records": len(new_coverage or {}),
            "v3_source": API_URL,
            "v3_dataset": "gnomad_r3 genome coverage",
        })
        coverage_payload["coverage_by_position"] = merged_coverage
        _atomic_json(COVERAGE_PATH, coverage_payload)
        COVERAGE_PROGRESS_PATH.unlink(missing_ok=True)
        output_paths.append(COVERAGE_PATH)
    for path in output_paths:
        print(path, hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_size)
    coverage_count = len(new_coverage or {}) if refresh_coverage else "preserved"
    print(
        f"v2 variants={metadata['v2_records']} v3 variants={metadata['v3_records']} "
        f"coverage_positions={coverage_count}"
    )


if __name__ == "__main__":
    raise SystemExit(
        "This legacy VCF/SciPy builder is audit-only and must not publish runtime data. "
        "Use: python scripts/refresh_gnomad_panel_snapshot.py refresh"
    )
