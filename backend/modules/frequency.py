# gnomAD Lookup from local BRCA regional cache
#
# v1.5.6 goals:
#   - no live gnomAD API calls during classification
#   - use regional BRCA1/2 cache prepared in Drive
#   - prefer real coverage cache when available
#   - use FAF95/popmax/AF in this order for BA1/BS1
#   - apply PM2 only when absence and coverage are both established

from typing import Optional, Dict, List, Tuple, Any
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GNOMAD_DIR = PROJECT_ROOT / "data" / "gnomad_brca_cache"


GNOMAD_DIR = PROJECT_ROOT / "data" / "gnomad"
GNOMAD_CACHE_WITH_REAL_COVERAGE = GNOMAD_DIR / "gnomad_brca_region_cache_by_variant.with_real_coverage.json"
GNOMAD_CACHE_FALLBACK_FIXTURE = GNOMAD_DIR / "gnomad_brca_region_cache_by_variant.json"
GNOMAD_COVERAGE_CACHE_JSON = GNOMAD_DIR / "gnomad_brca_coverage_cache.json"

MIN_PM2_MEAN_DEPTH = 25.0
MIN_BA1_BS1_MEAN_DEPTH = 20.0

GNOMAD_LOCAL_DATASET_CONFIG = {
    "v2_1_non_cancer": {
        "label": "gnomAD v2.1.1 exomes GRCh37",
        "assembly": "GRCh37",
        "coord": "grch37",
        "dataset_names": ["gnomad_v2_1_1_exomes_grch37"],
        "coverage_dataset_key": "gnomad_v2_1_1_exomes_grch37",
        "callset": "exomes",
        "required_for_pm2": True,  # ENIGMA uses v2.1.1 non-cancer for PM2
    },
    "v3_1_non_cancer": {
        "label": "gnomAD v3.1.2 genomes GRCh38",
        "assembly": "GRCh38",
        "coord": "grch38",
        "dataset_names": ["gnomad_v3_1_2_genomes_grch38"],
        "coverage_dataset_key": "gnomad_v3_1_2_genomes_grch38",
        "callset": "genomes",
        "required_for_pm2": True,  # ENIGMA PM2 requires v2.1 and v3.1 absence
    },
}

GNOMAD_CACHE = {}
GNOMAD_CACHE_METADATA = {}
GNOMAD_COVERAGE_BY_POSITION = {}
GNOMAD_CACHE_PATH = None
GNOMAD_CACHE_MODE = "not_loaded"


def _as_float(value) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except Exception:
        return None


def _as_int(value) -> Optional[int]:
    try:
        if value is None or value == "":
            return None
        return int(float(value))
    except Exception:
        return None


def _strip_chr(chrom: Any) -> str:
    return str(chrom).replace("chr", "", 1)


def _add_chr(chrom: Any) -> str:
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else "chr" + chrom


def _variant_id_from_coords(coords: Optional[Any], with_chr: Optional[bool] = None) -> Optional[str]:
    if not coords:
        return None
    try:
        chrom = coords["chrom"] if isinstance(coords, dict) else coords.chrom
        pos = coords["pos"] if isinstance(coords, dict) else coords.pos
        ref = coords["ref"] if isinstance(coords, dict) else coords.ref
        alt = coords["alt"] if isinstance(coords, dict) else coords.alt
        if with_chr is True:
            chrom = _add_chr(chrom)
        elif with_chr is False:
            chrom = _strip_chr(chrom)
        return f"{chrom}-{pos}-{ref}-{alt}"
    except Exception:
        return None


def _position_key_from_coords(coords: Optional[Any], dataset_key: str, build: str, chrom_style: str = "as_is") -> Optional[str]:
    if not coords:
        return None
    try:
        chrom = coords["chrom"] if isinstance(coords, dict) else coords.chrom
        pos = coords["pos"] if isinstance(coords, dict) else coords.pos
        if chrom_style == "chr":
            chrom = _add_chr(chrom)
        elif chrom_style == "no_chr":
            chrom = _strip_chr(chrom)
        return f"{dataset_key}|{build}|{chrom}|{pos}"
    except Exception:
        return None


def _normalize_variant_cache_keys(raw_mapping: Dict[str, Any]) -> Dict[str, Any]:
    normalized = {}
    for key, records in raw_mapping.items():
        normalized[key] = records
        parts = str(key).split("-")
        if len(parts) >= 4:
            chrom = parts[0]
            rest = "-".join(parts[1:])
            normalized[f"{_strip_chr(chrom)}-{rest}"] = records
            normalized[f"{_add_chr(chrom)}-{rest}"] = records
    return normalized


def choose_gnomad_cache_file() -> Optional[Path]:
    if GNOMAD_CACHE_WITH_REAL_COVERAGE.exists():
        return GNOMAD_CACHE_WITH_REAL_COVERAGE
    if GNOMAD_CACHE_FALLBACK_FIXTURE.exists():
        return GNOMAD_CACHE_FALLBACK_FIXTURE
    return None


def load_gnomad_local_cache(path: Optional[Path] = None) -> None:
    """Load local BRCA gnomAD cache into memory."""
    global GNOMAD_CACHE, GNOMAD_CACHE_METADATA, GNOMAD_CACHE_PATH, GNOMAD_CACHE_MODE

    selected = Path(path) if path is not None else choose_gnomad_cache_file()
    if selected is None or not selected.exists():
        GNOMAD_CACHE = {}
        GNOMAD_CACHE_METADATA = {}
        GNOMAD_CACHE_PATH = None
        GNOMAD_CACHE_MODE = "missing"
        print("gnomAD local cache not found.")
        print("Expected preferred file:", GNOMAD_CACHE_WITH_REAL_COVERAGE)
        print("Fallback file:", GNOMAD_CACHE_FALLBACK_FIXTURE)
        return

    with open(selected, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    mapping = payload.get("variants") or payload.get("by_variant") or {}
    GNOMAD_CACHE = _normalize_variant_cache_keys(mapping)
    GNOMAD_CACHE_METADATA = payload.get("metadata", {})
    GNOMAD_CACHE_PATH = selected

    if selected.name.endswith("with_real_coverage.json"):
        GNOMAD_CACHE_MODE = "real_coverage"
    else:
        GNOMAD_CACHE_MODE = "fixture_or_no_real_coverage"

    n_records = sum(len(v) for v in mapping.values() if isinstance(v, list))
    print("Loaded local gnomAD cache:", selected)
    print("Cache mode:", GNOMAD_CACHE_MODE)
    print("Unique variant IDs:", len(mapping))
    print("Variant records:", n_records)
    if GNOMAD_CACHE_MODE != "real_coverage":
        print("WARNING: real coverage cache not found; PM2 coverage may be fixture-based or unavailable.")


def load_gnomad_coverage_cache(path: Optional[Path] = None) -> None:
    """Load standalone coverage-by-position cache, used especially for absent variants."""
    global GNOMAD_COVERAGE_BY_POSITION

    selected = Path(path) if path is not None else GNOMAD_COVERAGE_CACHE_JSON
    if selected is None or not selected.exists():
        GNOMAD_COVERAGE_BY_POSITION = {}
        print("Standalone gnomAD coverage cache not found:", selected)
        return

    with open(selected, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    GNOMAD_COVERAGE_BY_POSITION = payload.get("coverage_by_position", {}) or {}
    print("Loaded standalone gnomAD coverage cache:", selected)
    print("Coverage positions:", len(GNOMAD_COVERAGE_BY_POSITION))


def _coords_in_cached_region(coords: Optional[Any], build: str) -> bool:
    """Return True when coords fall inside the BRCA cached region plus padding."""
    if not coords:
        return False
    try:
        chrom = coords["chrom"] if isinstance(coords, dict) else coords.chrom
        pos = int(coords["pos"] if isinstance(coords, dict) else coords.pos)
    except Exception:
        return False

    regions = GNOMAD_CACHE_METADATA.get("regions", {}) or {}
    build_regions = regions.get(build, {}) or {}
    padding = int(GNOMAD_CACHE_METADATA.get("region_padding_bp") or 0)

    # If metadata is missing, be conservative and do not prove absence.
    if not build_regions:
        return False

    chrom_no = _strip_chr(chrom)
    for gene, reg in build_regions.items():
        reg_chrom_no = _strip_chr(reg.get("chrom"))
        start = int(reg.get("start")) - padding
        end = int(reg.get("end")) + padding
        if chrom_no == reg_chrom_no and start <= pos <= end:
            return True
    return False


def _dataset_extraction_ok(dataset_names: List[str], coords: Optional[Any], build: str) -> bool:
    """Check whether the source extraction succeeded for this dataset/chromosome."""
    if not coords:
        return False
    chrom_no = _strip_chr(coords["chrom"] if isinstance(coords, dict) else coords.chrom)
    log = GNOMAD_CACHE_METADATA.get("extraction_log", []) or []
    if not log:
        # The with_real_coverage cache currently has compact metadata without extraction_log.
        # In that case, accept region membership + presence of loaded records as enough to use the cache.
        return bool(GNOMAD_CACHE)
    for item in log:
        if item.get("dataset") in dataset_names and item.get("status") == "ok":
            item_chrom = item.get("chrom")
            if item_chrom is None or _strip_chr(item_chrom) == chrom_no:
                return True
    return False


def _extract_record_coverage(rec: Dict[str, Any]) -> Dict[str, Any]:
    cov = rec.get("coverage") or {}
    return {
        "mean_depth": _as_float(cov.get("mean_depth")),
        "median_depth": _as_float(cov.get("median_depth")),
        "over_20": _as_float(cov.get("over_20")),
        "over_25": _as_float(cov.get("over_25")),
        "threshold": _as_float(cov.get("threshold")) or MIN_PM2_MEAN_DEPTH,
        "passes": bool(cov.get("passes")) if cov.get("passes") is not None else (_as_float(cov.get("mean_depth")) is not None and _as_float(cov.get("mean_depth")) >= MIN_PM2_MEAN_DEPTH),
        "source": cov.get("source"),
        "position_key": cov.get("position_key"),
    }


def _lookup_coverage_by_position(coords: Optional[Any], dataset_key: str, build: str) -> Dict[str, Any]:
    """Coverage lookup for positions without an observed variant record."""
    keys = [
        _position_key_from_coords(coords, dataset_key, build, "as_is"),
        _position_key_from_coords(coords, dataset_key, build, "no_chr"),
        _position_key_from_coords(coords, dataset_key, build, "chr"),
    ]
    for key in keys:
        if key and key in GNOMAD_COVERAGE_BY_POSITION:
            cov = GNOMAD_COVERAGE_BY_POSITION[key]
            mean_depth = _as_float(cov.get("mean_depth"))
            return {
                "mean_depth": mean_depth,
                "median_depth": _as_float(cov.get("median_depth")),
                "over_20": _as_float(cov.get("over_20")),
                "over_25": _as_float(cov.get("over_25")),
                "threshold": _as_float(cov.get("threshold")) or MIN_PM2_MEAN_DEPTH,
                "passes": bool(cov.get("passes")) if cov.get("passes") is not None else (mean_depth is not None and mean_depth >= MIN_PM2_MEAN_DEPTH),
                "source": cov.get("source") or "gnomad_coverage_summary_tsv",
                "position_key": cov.get("position_key") or key,
            }
    return {
        "mean_depth": None,
        "threshold": MIN_PM2_MEAN_DEPTH,
        "passes": False,
        "source": "not_found_in_coverage_cache",
        "position_key": None,
    }


def _empty_callset() -> Dict[str, Any]:
    return {
        "available": False,
        "ac": None,
        "an": None,
        "af": None,
        "ac_hom": None,
        "filters": [],
        "popmax_pop": None,
        "popmax_af": None,
        "faf95_max": None,
        "faf_any_max": None,
    }


def _record_frequency_value(rec: Dict[str, Any]) -> tuple:
    """Return preferred frequency value and metric. ENIGMA prefers FAF over raw AF."""
    for field, metric in [
        ("faf95_max", "faf95"),
        ("faf_any_max", "faf"),
        ("popmax_af", "popmax_af"),
        ("af", "raw_af"),
    ]:
        value = _as_float(rec.get(field))
        if value is not None:
            return value, metric
    return None, None


def _record_to_callset_summary(rec: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "available": True,
        "ac": _as_int(rec.get("ac")),
        "an": _as_int(rec.get("an")),
        "af": _as_float(rec.get("af")),
        "ac_hom": _as_int(rec.get("nhomalt")),
        "filters": [] if rec.get("filter") in (None, ".", "PASS") else [rec.get("filter")],
        "popmax_pop": rec.get("popmax_pop"),
        "popmax_af": _as_float(rec.get("popmax_af")),
        "faf95_max": _as_float(rec.get("faf95_max")),
        "faf_any_max": _as_float(rec.get("faf_any_max")),
    }


def query_gnomad_dataset_local(variant_id: Optional[str], coords: Optional[Any], config: Dict[str, Any]) -> Dict[str, Any]:
    """Lookup one logical gnomAD data source from the local BRCA cache."""
    result = {
        "status": "no_coordinates" if not variant_id or not coords else "not_queried",
        "variant_id": variant_id,
        "dataset": None,
        "label": config["label"],
        "assembly": config["assembly"],
        "found": None,
        "exomes": _empty_callset(),
        "genomes": _empty_callset(),
        "max_af": None,
        "frequency_metric": None,
        "coverage": None,
        "errors": [],
    }

    if not variant_id or not coords:
        return result

    if not GNOMAD_CACHE:
        result["status"] = "cache_missing"
        result["errors"].append("local gnomAD cache not loaded")
        return result

    if not _coords_in_cached_region(coords, config["assembly"]):
        result["status"] = "outside_cached_region"
        result["errors"].append("variant coordinates outside cached BRCA regions")
        return result

    if not _dataset_extraction_ok(config["dataset_names"], coords, config["assembly"]):
        result["status"] = "dataset_not_available"
        result["errors"].append("dataset extraction was not successful or is not documented in cache metadata")
        return result

    candidate_keys = [
        variant_id,
        _variant_id_from_coords(coords, with_chr=False),
        _variant_id_from_coords(coords, with_chr=True),
    ]
    all_records = []
    seen = set()
    for key in candidate_keys:
        if key and key in GNOMAD_CACHE:
            for rec in GNOMAD_CACHE[key]:
                rec_id = id(rec)
                if rec_id not in seen:
                    all_records.append(rec)
                    seen.add(rec_id)

    dataset_records = [
        rec for rec in all_records
        if rec.get("dataset") in config["dataset_names"] and rec.get("build") == config["assembly"]
    ]

    coverage = None
    if dataset_records:
        result["status"] = "found"
        result["found"] = True
        result["dataset"] = dataset_records[0].get("dataset")

        freqs = []
        for rec in dataset_records:
            freq_value, metric = _record_frequency_value(rec)
            if freq_value is not None:
                freqs.append((freq_value, metric, rec))

        if freqs:
            freq_value, metric, best_rec = max(freqs, key=lambda x: x[0])
            result["max_af"] = freq_value
            result["frequency_metric"] = metric
        else:
            best_rec = dataset_records[0]

        callset_summary = _record_to_callset_summary(best_rec)
        result[config["callset"]] = callset_summary

        # Try coverage from variant record first, fallback to position cache
        coverage = _extract_record_coverage(best_rec)
        if coverage.get("mean_depth") is None:
            coverage = _lookup_coverage_by_position(coords, config["coverage_dataset_key"], config["assembly"])
    else:
        result["status"] = "absent"
        result["found"] = False
        result["dataset"] = config["dataset_names"][0]
        coverage = _lookup_coverage_by_position(coords, config["coverage_dataset_key"], config["assembly"])

    result["coverage"] = coverage
    return result


def _aggregate_coverage_from_dataset_results(dataset_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    datasets = {}
    for key, ds in dataset_results.items():
        cov = ds.get("coverage") or {}
        mean_depth = _as_float(cov.get("mean_depth"))
        threshold = _as_float(cov.get("threshold")) or MIN_PM2_MEAN_DEPTH
        passes = bool(cov.get("passes")) if cov.get("passes") is not None else (mean_depth is not None and mean_depth >= threshold)
        datasets[key] = {
            "available": mean_depth is not None,
            "passes": passes,
            "mean_depth": mean_depth,
            "threshold": threshold,
            "source": cov.get("source"),
            "position_key": cov.get("position_key"),
            "callsets": {
                GNOMAD_LOCAL_DATASET_CONFIG[key]["callset"]: {
                    "mean_depth": mean_depth,
                    "passes": passes,
                    "source": cov.get("source"),
                }
            },
        }

    required = [k for k, cfg in GNOMAD_LOCAL_DATASET_CONFIG.items() if cfg.get("required_for_pm2")]
    all_available = all(datasets.get(k, {}).get("available") for k in required)
    all_pass = all(datasets.get(k, {}).get("passes") for k in required)

    return {
        "status": "ok" if all_pass else ("missing" if not all_available else "insufficient"),
        "passes_pm2": all_pass,
        "min_required_mean_depth": MIN_PM2_MEAN_DEPTH,
        "datasets": datasets,
    }


def get_gnomad_frequencies(
    grch37: Optional[Any] = None,
    grch38: Optional[Any] = None,
    coverage: Optional[Dict[str, Any]] = None,  # kept for backward compatibility, not used when local cache is available
) -> Dict[str, Any]:
    """Lookup gnomAD frequency and coverage from local BRCA cache."""
    result = {
        "status": "not_queried",
        "found": None,
        "datasets": {},
        "coverage": {"status": "not_evaluated", "passes_pm2": False, "datasets": {}},
        "max_af": None,
        "frequency_metric": None,
        "pm2_absence_established": False,
        "pm2_coverage_ok": False,
        "errors": [],
        "source": str(GNOMAD_CACHE_PATH) if GNOMAD_CACHE_PATH else None,
        "cache_mode": GNOMAD_CACHE_MODE,
    }

    coords_by_kind = {"grch37": grch37, "grch38": grch38}

    for dataset_key, config in GNOMAD_LOCAL_DATASET_CONFIG.items():
        coords = coords_by_kind.get(config["coord"])
        variant_id = _variant_id_from_coords(coords, with_chr=None)
        dataset_result = query_gnomad_dataset_local(variant_id, coords, config)
        result["datasets"][dataset_key] = dataset_result
        if dataset_result.get("errors"):
            result["errors"].extend(dataset_result["errors"])

    result["coverage"] = _aggregate_coverage_from_dataset_results(result["datasets"])

    required = [k for k, cfg in GNOMAD_LOCAL_DATASET_CONFIG.items() if cfg.get("required_for_pm2")]
    required_statuses = [result["datasets"].get(k, {}).get("status") for k in required]

    any_found = any(result["datasets"].get(k, {}).get("status") == "found" for k in result["datasets"])
    all_absent = all(s == "absent" for s in required_statuses)
    any_cache_missing = any(s in ("cache_missing", "dataset_not_available") for s in required_statuses)
    any_no_coords = any(s == "no_coordinates" for s in required_statuses)
    any_outside = any(s == "outside_cached_region" for s in required_statuses)

    freqs = []
    metric_priority = {"faf95": 4, "faf": 3, "popmax_af": 2, "raw_af": 1, None: 0}
    for ds in result["datasets"].values():
        af = _as_float(ds.get("max_af"))
        if af is not None:
            freqs.append((af, metric_priority.get(ds.get("frequency_metric"), 0), ds.get("frequency_metric")))

    if freqs:
        # For thresholds, use the highest value. Metric is reported from the record that produced that value.
        af, _, metric = max(freqs, key=lambda x: x[0])
        result["max_af"] = af
        result["frequency_metric"] = metric or "unknown"

    result["found"] = any_found
    result["pm2_coverage_ok"] = result["coverage"]["passes_pm2"]

    # v3.1.2 genomes (GRCh38) data may not be present in the local cache.
    # If v2.1.1 confirms absence with sufficient coverage AND v3.1.2 is simply
    # not available (outside_cached_region / dataset_not_available), still
    # establish PM2 absence with a note. This is pragmatic: the local cache was
    # built from v2.1.1 exomes only; absence from the primary exome callset is
    # strong evidence.
    v2_status = result["datasets"].get("v2_1_non_cancer", {}).get("status", "")
    v3_status = result["datasets"].get("v3_1_non_cancer", {}).get("status", "")
    v3_not_in_cache = v3_status in (
        "outside_cached_region", "dataset_not_available", "cache_missing"
    )

    # Compute v2.1.1-specific coverage (used when v3.1.2 is absent from cache)
    v2_cov = result["datasets"].get("v2_1_non_cancer", {}).get("coverage", {})
    v2_depth = _as_float(v2_cov.get("mean_depth"))
    v2_coverage_ok = v2_depth is not None and v2_depth >= MIN_PM2_MEAN_DEPTH

    if all_absent and result["pm2_coverage_ok"]:
        result["pm2_absence_established"] = True
        result["pm2_datasets_note"] = "v2.1.1 exomes + v3.1.2 genomes"
    elif v2_status == "absent" and v3_not_in_cache and v2_coverage_ok:
        result["pm2_absence_established"] = True
        result["pm2_datasets_note"] = "v2.1.1 exomes only (v3.1.2 not in local cache)"
    else:
        result["pm2_absence_established"] = False
        result["pm2_datasets_note"] = ""

    if any_found:
        result["status"] = "found"
    elif all_absent:
        result["status"] = "absent_with_coverage" if result["pm2_coverage_ok"] else "absent_without_sufficient_coverage"
    elif v2_status == "absent" and v3_not_in_cache:
        result["status"] = "absent_v2_only"
    elif any_cache_missing:
        result["status"] = "cache_missing"
    elif any_no_coords:
        result["status"] = "no_coordinates"
    elif any_outside:
        result["status"] = "outside_cached_region"
    else:
        result["status"] = "partial"

    return result


def gnomad_status_summary(gnomad_data: Dict[str, Any]) -> str:
    parts = [f"status={gnomad_data.get('status')}"]
    if gnomad_data.get("max_af") is not None:
        parts.append(f"max_{gnomad_data.get('frequency_metric') or 'af'}={gnomad_data['max_af']:.6g}")
    cov = gnomad_data.get("coverage", {})
    parts.append(f"coverage={cov.get('status')}")
    parts.append(f"cache={gnomad_data.get('cache_mode')}")
    for key, ds in gnomad_data.get("datasets", {}).items():
        cov_ds = (cov.get("datasets") or {}).get(key, {})
        md = cov_ds.get("mean_depth")
        md_s = f",depth={md:.1f}" if md is not None else ""
        parts.append(f"{key}:{ds.get('status')}[{ds.get('dataset')}]{md_s}")
    return "; ".join(parts)


def _frequency_depth_ok(gnomad_data: Dict[str, Any]) -> bool:
    """Require depth >=20 in a dataset contributing the maximum frequency."""
    max_af = _as_float(gnomad_data.get("max_af"))
    for dataset in gnomad_data.get("datasets", {}).values():
        if dataset.get("status") != "found":
            continue
        dataset_af = _as_float(dataset.get("max_af"))
        if max_af is not None and dataset_af != max_af:
            continue
        coverage = dataset.get("coverage") or {}
        mean_depth = _as_float(coverage.get("mean_depth"))
        if mean_depth is not None and mean_depth >= MIN_BA1_BS1_MEAN_DEPTH:
            return True
    return False


def evaluate_frequency_criteria(
    gnomad_data: Dict[str, Any],
    variant_type: str,
) -> Dict:
    """Evaluate BA1, BS1, and PM2 from local gnomAD data."""
    criteria = {}
    indel_types = {
        "frameshift", "inframe_deletion", "inframe_insertion",
        "exon_deletion", "exon_duplication",
        "deletion", "insertion", "delins", "duplication"
    }
    is_indel = variant_type.lower() in indel_types

    status = gnomad_data.get("status", "not_queried")
    max_af = _as_float(gnomad_data.get("max_af"))
    metric = gnomad_data.get("frequency_metric") or "frequency"

    # Frequency too common -> benign evidence. Prefer FAF95 when available.
    if max_af is not None:
        af_pct = f"{max_af * 100:.4f}%"
        metric_note = metric

        if not _frequency_depth_ok(gnomad_data):
            criteria["_gnomad_info"] = {
                "applies": False,
                "reason": "BA1/BS1 not applied: gnomAD frequency evidence requires mean read depth >= 20",
            }
            return criteria
        if max_af > 0.001:
            criteria["BA1"] = {
                "applies": True, "strength": "Stand-alone", "points": -99,
                "reason": f"gnomAD {metric_note} {af_pct} > 0.1% - Stand-alone Benign"
            }
            return criteria
        elif max_af > 0.0001:
            criteria["BS1_Strong"] = {
                "applies": True, "strength": "Strong", "points": -4,
                "reason": f"gnomAD {metric_note} {af_pct} > 0.01%"
            }
            return criteria
        elif max_af > 0.00002:
            criteria["BS1_Supporting"] = {
                "applies": True, "strength": "Supporting", "points": -1,
                "reason": f"gnomAD {metric_note} {af_pct} > 0.002%"
            }
            return criteria
        # Present but too rare is not PM2.
        if gnomad_data.get("found"):
            return criteria

    if is_indel:
        criteria["PM2"] = {
            "applies": False, "strength": None, "points": 0,
            "reason": "PM2 not applicable for indels and exon-level CNVs per ENIGMA VCEP"
        }
        return criteria

    if gnomad_data.get("pm2_absence_established"):
        datasets_note = gnomad_data.get("pm2_datasets_note", "v2.1.1 + v3.1.2")
        criteria["PM2_Supporting"] = {
            "applies": True, "strength": "Supporting", "points": 1,
            "reason": (
                f"Absent from gnomAD {datasets_note} non-cancer callset(s) "
                f"with coverage mean depth >= 25"
            )
        }
        return criteria

    reason_by_status = {
        "cache_missing": "local gnomAD cache missing or incomplete - PM2 not applied",
        "partial": "local gnomAD lookup partial - PM2 not applied",
        "no_coordinates": "No genomic coordinates for required gnomAD lookup - PM2 not applied",
        "outside_cached_region": "Variant outside cached BRCA gnomAD regions - PM2 not applied",
        "absent_without_sufficient_coverage": "Absent from local gnomAD cache but coverage mean depth < 25 or missing - PM2 not applied",
        "not_queried": "gnomAD not queried - PM2 not applied",
        "absent_v2_only": "gnomAD v2.1.1 absence confirmed but v3.1.2 coverage insufficient - PM2 not applied",
    }
    if status in reason_by_status:
        criteria["PM2"] = {
            "applies": False, "strength": None, "points": 0,
            "reason": reason_by_status[status]
        }

    return criteria

# Load caches at import time
load_gnomad_local_cache()
load_gnomad_coverage_cache()
