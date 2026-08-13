# gnomAD lookup from a local gene-panel regional cache
#
# v1.5.6 goals:
#   - no live gnomAD API calls during classification
#   - use the checksum-verified regional panel cache
#   - require the checksum-verified frequency and coverage snapshots
#   - use only non-cancer FAF95 for BA1/BS1; never substitute raw/popmax AF
#   - apply PM2 only when absence and coverage are both established

from typing import Optional, Dict, List, Tuple, Any
from pathlib import Path
import hashlib
import json
import re
from backend.data_health import clear_issue, register_issue

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GNOMAD_DIR = PROJECT_ROOT / "data" / "gnomad"
# The historical filename is retained for deployment compatibility. Its
# contents and accepted regions are governed entirely by the panel manifest.
GNOMAD_CACHE_WITH_REAL_COVERAGE = GNOMAD_DIR / "gnomad_brca_region_cache_by_variant.with_real_coverage.json"
GNOMAD_COVERAGE_CACHE_JSON = GNOMAD_DIR / "gnomad_brca_coverage_cache.json"
GNOMAD_PANEL_MANIFEST_JSON = GNOMAD_DIR / "gnomad_panel_manifest.json"

GNOMAD_LOCAL_DATASET_CONFIG = {
    "v2_1_non_cancer": {
        "label": "gnomAD v2.1.1 exomes GRCh37",
        "assembly": "GRCh37",
        "coord": "grch37",
        "dataset_names": ["gnomad_v2_1_1_exomes_grch37"],
        "coverage_dataset_key": "gnomad_v2_1_1_exomes_grch37",
        "callset": "exomes",
    },
    "v3_1_non_cancer": {
        "label": "gnomAD v3.1.2 genomes GRCh38",
        "assembly": "GRCh38",
        "coord": "grch38",
        "dataset_names": ["gnomad_v3_1_2_genomes_grch38"],
        "coverage_dataset_key": "gnomad_v3_1_2_genomes_grch38",
        "callset": "genomes",
    },
}

GNOMAD_CACHE = {}
GNOMAD_CACHE_METADATA = {}
GNOMAD_COVERAGE_BY_POSITION = {}
GNOMAD_CACHE_PATH = None
GNOMAD_CACHE_MODE = "not_loaded"


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _approved_manifest_sha256() -> Optional[str]:
    try:
        manifest = json.loads(GNOMAD_PANEL_MANIFEST_JSON.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return _canonical_sha256(manifest)


def _approved_manifest() -> Optional[Dict[str, Any]]:
    try:
        manifest = json.loads(GNOMAD_PANEL_MANIFEST_JSON.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return manifest if isinstance(manifest, dict) else None


def _classification_policy_for_gene(
    gene: Optional[str], manifest: Optional[Dict[str, Any]] = None
) -> Optional[Dict[str, Any]]:
    """Return the explicitly assigned active policy; never inherit one."""
    if not gene:
        return None
    manifest = manifest or _approved_manifest()
    if not manifest or manifest.get("schema_version") != 2:
        return None
    target = next(
        (
            item
            for item in manifest.get("targets", [])
            if str(item.get("gene", "")).upper() == str(gene).upper()
        ),
        None,
    )
    if not target:
        return None
    if target.get("activation_status") != "active":
        return None
    policy_id = target.get("classification_policy_id")
    policy = (manifest.get("classification_policies") or {}).get(policy_id)
    if not policy or policy.get("status") != "active":
        return None
    if str(gene).upper() not in {
        str(item).upper() for item in policy.get("applicable_genes", [])
    }:
        return None
    return {"policy_id": policy_id, **policy}


def _manifest_scored_ancestries(manifest: Dict[str, Any]) -> set:
    return {
        ancestry
        for policy in (manifest.get("classification_policies") or {}).values()
        if policy.get("status") == "active"
        for ancestry in (
            policy.get("frequency_criteria", {}).get(
                "scored_non_founder_ancestries", []
            )
        )
    }


def _runtime_dataset_binding_error(manifest: Dict[str, Any]) -> Optional[str]:
    """Verify that policy runtime keys resolve to the pinned manifest products."""
    manifest_by_runtime = {
        item.get("runtime_key"): item for item in manifest.get("datasets", [])
    }
    required = {
        runtime_key
        for policy in (manifest.get("classification_policies") or {}).values()
        if policy.get("status") == "active"
        for runtime_key in (
            policy.get("frequency_criteria", {}).get(
                "required_dataset_runtime_keys", []
            )
        )
    }
    for runtime_key in required:
        item = manifest_by_runtime.get(runtime_key)
        config = GNOMAD_LOCAL_DATASET_CONFIG.get(runtime_key)
        if item is None or config is None:
            return f"runtime dataset {runtime_key!r} is not implemented"
        if (
            item.get("dataset_key") not in config["dataset_names"]
            or item.get("assembly") != config["assembly"]
            or item.get("callset") != config["callset"]
            or item.get("subset") != "non_cancer"
            or item.get("dataset_key") != config["coverage_dataset_key"]
        ):
            return f"runtime dataset {runtime_key!r} differs from the manifest"
    return None


def _validate_gnomad_cache_payload(payload: Dict[str, Any]) -> Optional[str]:
    """Validate that every scored gnomAD record has auditable non-cancer FAF95."""
    mapping = payload.get("variants") or payload.get("by_variant") or {}
    if not isinstance(mapping, dict) or not mapping:
        return "variant mapping is missing or empty"

    metadata = payload.get("metadata") or {}
    if metadata.get("schema_version") != 2:
        return "snapshot schema_version must be 2"
    approved_manifest = _approved_manifest()
    if approved_manifest is None:
        return f"approved panel manifest is missing or invalid: {GNOMAD_PANEL_MANIFEST_JSON}"
    binding_error = _runtime_dataset_binding_error(approved_manifest)
    if binding_error:
        return binding_error
    approved_manifest_hash = _canonical_sha256(approved_manifest)
    if metadata.get("manifest_sha256") != approved_manifest_hash:
        return "snapshot was built from a different panel/source manifest"
    if metadata.get("automatic_release_activation") is not False:
        return "automatic release activation must be disabled"
    if metadata.get("records_sha256") != _canonical_sha256(mapping):
        return "variant records checksum mismatch"
    if metadata.get("classification_policies") != approved_manifest.get(
        "classification_policies"
    ):
        return "snapshot classification policies differ from the approved manifest"

    expected_datasets = {
        name
        for config in GNOMAD_LOCAL_DATASET_CONFIG.values()
        for name in config["dataset_names"]
    }
    counts = {dataset: 0 for dataset in expected_datasets}
    missing_faf95 = []
    expected_methods = {
        "gnomad_v2_1_1_exomes_grch37": "official_gnomad_hail_table_non_cancer_faf95",
        "gnomad_v3_1_2_genomes_grch38": (
            "hail.experimental.filtering_allele_frequency_from_official_non_cancer_ac_an"
        ),
    }
    manifest_datasets = {
        item["dataset_key"]: item for item in approved_manifest.get("datasets", [])
    }
    scored_ancestries = _manifest_scored_ancestries(approved_manifest)
    if not scored_ancestries:
        return "approved manifest contains no active scored population groups"
    for variant_id, records in mapping.items():
        if not isinstance(records, list):
            return f"record list is invalid for {variant_id}"
        for record in records:
            dataset = record.get("dataset")
            if dataset not in expected_datasets:
                continue
            counts[dataset] += 1
            if (
                _as_float(record.get("faf95_max")) is None
                or record.get("faf95_scope") != "non_cancer_non_founder_ancestries"
                or record.get("faf95_method") != expected_methods.get(dataset)
            ):
                missing_faf95.append(variant_id)
                if len(missing_faf95) >= 3:
                    break
                continue
            if set(record.get("non_founder_ac_by_ancestry") or {}) != scored_ancestries:
                return f"non-founder population counts are incomplete for {variant_id}"
            observed = any(
                (_as_int(value) or 0) > 0
                for value in record["non_founder_ac_by_ancestry"].values()
            )
            if record.get("non_founder_observed") is not observed:
                return f"non-founder presence flag is inconsistent for {variant_id}"
            expected_context = {
                item["code"]
                for item in manifest_datasets.get(dataset, {}).get(
                    "excluded_population_context", []
                )
            }
            context = record.get("excluded_population_context") or {}
            if set(context) != expected_context or any(
                item.get("used_for_ba1_bs1") is not False
                or item.get("used_for_pm2_presence") is not False
                for item in context.values()
            ):
                return f"excluded population context is invalid for {variant_id}"
        if len(missing_faf95) >= 3:
            break

    absent_datasets = sorted(dataset for dataset, count in counts.items() if count == 0)
    if absent_datasets:
        return f"required dataset records are missing: {', '.join(absent_datasets)}"
    if missing_faf95:
        return (
            "records lack ENIGMA-compatible non-cancer FAF95 provenance: "
            + ", ".join(missing_faf95)
        )

    for release in ("v2_faf95", "v3_faf95"):
        item = metadata.get(release) or {}
        if item.get("raw_af_fallback_allowed") is not False:
            return f"metadata {release}.raw_af_fallback_allowed must be false"

    logged_datasets = {
        item.get("dataset")
        for item in metadata.get("extraction_log", [])
        if item.get("status") == "ok"
        and (item.get("source_identity") or {}).get("x_goog_hash")
        and (item.get("source_identity") or {}).get("etag")
    }
    if not expected_datasets.issubset(logged_datasets):
        return "official source identity or successful extraction log is incomplete"
    return None


def _validate_gnomad_coverage_payload(payload: Dict[str, Any]) -> Optional[str]:
    mapping = payload.get("coverage_by_position") or {}
    if not isinstance(mapping, dict) or not mapping:
        return "coverage mapping is missing or empty"
    metadata = payload.get("metadata") or {}
    if metadata.get("schema_version") != 2:
        return "coverage schema_version must be 2"
    approved_manifest_hash = _approved_manifest_sha256()
    if approved_manifest_hash is None:
        return f"approved panel manifest is missing or invalid: {GNOMAD_PANEL_MANIFEST_JSON}"
    if metadata.get("manifest_sha256") != approved_manifest_hash:
        return "coverage was built from a different panel/source manifest"
    approved_manifest = _approved_manifest()
    if approved_manifest is not None:
        binding_error = _runtime_dataset_binding_error(approved_manifest)
        if binding_error:
            return binding_error
    if approved_manifest is None or metadata.get(
        "classification_policies"
    ) != approved_manifest.get("classification_policies"):
        return "coverage classification policies differ from the approved manifest"
    if metadata.get("records") != len(mapping):
        return "coverage record count does not match metadata"
    if metadata.get("records_sha256") != _canonical_sha256(mapping):
        return "coverage records checksum mismatch"

    expected_datasets = {
        name
        for config in GNOMAD_LOCAL_DATASET_CONFIG.values()
        for name in config["dataset_names"]
    }
    observed_datasets = {item.get("dataset_key") for item in mapping.values()}
    if observed_datasets != expected_datasets:
        return "coverage does not contain exactly the required gnomAD datasets"
    documented_datasets = {
        item.get("dataset")
        for item in metadata.get("datasets", [])
        if (item.get("coverage_source_identity") or {}).get("etag")
        and (item.get("coverage_source_identity") or {}).get("x_goog_hash")
        and item.get("coverage_hail_uri")
    }
    if documented_datasets != expected_datasets:
        return "official coverage source identity is incomplete"
    return None


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


def _reference_span_from_coords(coords: Optional[Any]) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """Return the genomic span represented by the VCF REF allele.

    ENIGMA v1.2 specifies an average depth for the region around a variant but
    does not define a flanking-window size. ARIANE therefore uses the only
    non-arbitrary reproducible region available from the variant itself: every
    genomic base represented by REF. For an SNV this is exactly one base.
    """
    if not coords:
        return None, None, None
    try:
        chrom = coords["chrom"] if isinstance(coords, dict) else coords.chrom
        pos = int(coords["pos"] if isinstance(coords, dict) else coords.pos)
        ref = coords.get("ref") if isinstance(coords, dict) else coords.ref
        ref_length = max(len(str(ref or "")), 1)
        return str(chrom), pos, pos + ref_length - 1
    except Exception:
        return None, None, None


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
        print("Expected approved file:", GNOMAD_CACHE_WITH_REAL_COVERAGE)
        register_issue(
            "gnomAD variant cache",
            f"approved real-coverage cache is missing: {GNOMAD_CACHE_WITH_REAL_COVERAGE}",
        )
        return

    try:
        with open(selected, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        GNOMAD_CACHE = {}
        GNOMAD_CACHE_METADATA = {}
        GNOMAD_CACHE_PATH = None
        GNOMAD_CACHE_MODE = "load_failed"
        register_issue(
            "gnomAD variant cache",
            f"could not load {selected}: {type(exc).__name__}: {exc}",
        )
        return

    validation_error = _validate_gnomad_cache_payload(payload)
    if validation_error:
        GNOMAD_CACHE = {}
        GNOMAD_CACHE_METADATA = payload.get("metadata", {})
        GNOMAD_CACHE_PATH = selected
        GNOMAD_CACHE_MODE = "invalid_faf95"
        register_issue(
            "gnomAD variant cache",
            f"cache failed non-cancer FAF95 validation: {validation_error}",
        )
        return

    mapping = payload.get("variants") or payload.get("by_variant") or {}
    GNOMAD_CACHE = _normalize_variant_cache_keys(mapping)
    GNOMAD_CACHE_METADATA = payload.get("metadata", {})
    GNOMAD_CACHE_PATH = selected

    GNOMAD_CACHE_MODE = "approved_snapshot"

    n_records = sum(len(v) for v in mapping.values() if isinstance(v, list))
    print("Loaded local gnomAD cache:", selected)
    print("Cache mode:", GNOMAD_CACHE_MODE)
    print("Unique variant IDs:", len(mapping))
    print("Variant records:", n_records)
    clear_issue("gnomAD variant cache")


def load_gnomad_coverage_cache(path: Optional[Path] = None) -> None:
    """Load standalone coverage-by-position cache, used especially for absent variants."""
    global GNOMAD_COVERAGE_BY_POSITION

    selected = Path(path) if path is not None else GNOMAD_COVERAGE_CACHE_JSON
    if selected is None or not selected.exists():
        GNOMAD_COVERAGE_BY_POSITION = {}
        print("Standalone gnomAD coverage cache not found:", selected)
        register_issue("gnomAD coverage cache", f"coverage cache is missing: {selected}")
        return

    try:
        with open(selected, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        GNOMAD_COVERAGE_BY_POSITION = {}
        register_issue(
            "gnomAD coverage cache",
            f"could not load {selected}: {type(exc).__name__}: {exc}",
        )
        return

    validation_error = _validate_gnomad_coverage_payload(payload)
    if validation_error:
        GNOMAD_COVERAGE_BY_POSITION = {}
        register_issue(
            "gnomAD coverage cache",
            f"coverage cache validation failed: {validation_error}",
        )
        return

    GNOMAD_COVERAGE_BY_POSITION = payload.get("coverage_by_position", {}) or {}
    print("Loaded standalone gnomAD coverage cache:", selected)
    print("Coverage positions:", len(GNOMAD_COVERAGE_BY_POSITION))
    clear_issue("gnomAD coverage cache")


def _coords_in_cached_region(coords: Optional[Any], build: str) -> bool:
    """Return True when coords fall inside a manifest panel region plus padding."""
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
    for item in log:
        if item.get("dataset") in dataset_names and item.get("status") == "ok":
            item_chrom = item.get("chrom")
            if item_chrom is None or _strip_chr(item_chrom) == chrom_no:
                return True

    return False


def _lookup_coverage_by_position(
    coords: Optional[Any], dataset_key: str, build: str, threshold: float
) -> Dict[str, Any]:
    """Average gnomAD depth across the complete reference-allele span."""
    chrom, span_start, span_end = _reference_span_from_coords(coords)
    if chrom is None or span_start is None or span_end is None:
        return {
            "mean_depth": None,
            "threshold": threshold,
            "passes": False,
            "source": "no_coordinates",
            "position_key": None,
            "position_keys": [],
            "coverage_scope": "variant_reference_span",
            "span_start": None,
            "span_end": None,
            "positions_expected": 0,
            "positions_available": 0,
        }

    found = []
    missing_positions = []
    for pos in range(span_start, span_end + 1):
        keys = [
            f"{dataset_key}|{build}|{chrom}|{pos}",
            f"{dataset_key}|{build}|{_strip_chr(chrom)}|{pos}",
            f"{dataset_key}|{build}|{_add_chr(chrom)}|{pos}",
        ]
        match = None
        match_key = None
        for key in dict.fromkeys(keys):
            if key in GNOMAD_COVERAGE_BY_POSITION:
                match = GNOMAD_COVERAGE_BY_POSITION[key]
                match_key = key
                break
        if match is None or _as_float(match.get("mean_depth")) is None:
            missing_positions.append(pos)
        else:
            found.append((match_key, match))

    expected = span_end - span_start + 1
    complete = len(found) == expected
    if complete:
        mean_depth = sum(_as_float(cov.get("mean_depth")) for _, cov in found) / expected
        numeric_medians = [_as_float(cov.get("median_depth")) for _, cov in found]
        numeric_over_20 = [_as_float(cov.get("over_20")) for _, cov in found]
        numeric_over_25 = [_as_float(cov.get("over_25")) for _, cov in found]
        median_depth = sum(numeric_medians) / expected if all(v is not None for v in numeric_medians) else None
        over_20 = sum(numeric_over_20) / expected if all(v is not None for v in numeric_over_20) else None
        over_25 = sum(numeric_over_25) / expected if all(v is not None for v in numeric_over_25) else None
        sources = sorted({cov.get("source") or "gnomad_coverage_summary_tsv" for _, cov in found})
        return {
            "mean_depth": mean_depth,
            "median_depth": median_depth,
            "over_20": over_20,
            "over_25": over_25,
            "threshold": threshold,
            "passes": mean_depth >= threshold,
            "source": ", ".join(sources),
            "position_key": found[0][1].get("position_key") or found[0][0],
            "position_keys": [cov.get("position_key") or key for key, cov in found],
            "coverage_scope": "variant_reference_span",
            "span_start": span_start,
            "span_end": span_end,
            "positions_expected": expected,
            "positions_available": len(found),
            "missing_positions": [],
        }
    return {
        "mean_depth": None,
        "threshold": threshold,
        "passes": False,
        "source": "not_found_in_coverage_cache",
        "position_key": None,
        "position_keys": [cov.get("position_key") or key for key, cov in found],
        "coverage_scope": "variant_reference_span",
        "span_start": span_start,
        "span_end": span_end,
        "positions_expected": expected,
        "positions_available": len(found),
        "missing_positions": missing_positions,
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


def _record_frequency_value(
    rec: Dict[str, Any], scored_ancestries: List[str]
) -> tuple:
    """Return ENIGMA-compatible non-cancer FAF95, without fallback metrics."""
    values = {
        code: _as_float((rec.get("faf95_by_ancestry") or {}).get(code))
        for code in scored_ancestries
    }
    if not values or any(value is None for value in values.values()):
        return None, None, None
    population, value = max(values.items(), key=lambda item: item[1])
    return value, "faf95", population


def _record_to_callset_summary(
    rec: Dict[str, Any], scored_ancestries: List[str]
) -> Dict[str, Any]:
    frequency_value, _, population = _record_frequency_value(
        rec, scored_ancestries
    )
    non_founder_ac = rec.get("non_founder_ac_by_ancestry") or {}
    return {
        "available": True,
        "ac": _as_int(rec.get("ac")),
        "an": _as_int(rec.get("an")),
        "af": _as_float(rec.get("af")),
        "ac_hom": _as_int(rec.get("nhomalt")),
        "filters": [] if rec.get("filter") in (None, ".", "PASS") else [rec.get("filter")],
        "popmax_pop": rec.get("popmax_pop"),
        "popmax_af": _as_float(rec.get("popmax_af")),
        "faf95_max": frequency_value,
        "faf95_pop": population,
        "faf95_scope": rec.get("faf95_scope"),
        "faf95_method": rec.get("faf95_method"),
        "faf95_by_ancestry": rec.get("faf95_by_ancestry") or {},
        "non_founder_ac_by_ancestry": non_founder_ac,
        "non_founder_allele_count": sum(
            _as_int(value) or 0 for value in non_founder_ac.values()
        ),
        "non_founder_an_by_ancestry": rec.get("non_founder_an_by_ancestry") or {},
        "non_founder_observed": rec.get("non_founder_observed"),
        "excluded_population_context": rec.get("excluded_population_context") or {},
        "faf_any_max": _as_float(rec.get("faf_any_max")),
    }


def _population_frequency_audit(gnomad_data: Dict[str, Any]) -> Dict[str, Any]:
    policy = gnomad_data.get("frequency_policy") or {}
    scored_codes = list(policy.get("scored_non_founder_ancestries") or [])
    datasets = []
    for dataset_key, result in gnomad_data.get("datasets", {}).items():
        config = GNOMAD_LOCAL_DATASET_CONFIG.get(dataset_key, {})
        callset = result.get(config.get("callset", "")) or {}
        scored_faf = callset.get("faf95_by_ancestry") or {}
        scored_ac = callset.get("non_founder_ac_by_ancestry") or {}
        scored_an = callset.get("non_founder_an_by_ancestry") or {}
        excluded = callset.get("excluded_population_context") or {}
        datasets.append(
            {
                "dataset_key": dataset_key,
                "label": result.get("label") or config.get("label", ""),
                "status": result.get("status"),
                "scored_max_faf95": result.get("max_af"),
                "scored_max_population": callset.get("faf95_pop"),
                "scored_non_founder_populations": [
                    {
                        "code": code,
                        "ac": scored_ac.get(code),
                        "an": scored_an.get(code),
                        "faf95": scored_faf.get(code),
                        "used_for_scoring": True,
                    }
                    for code in scored_codes
                ],
                "excluded_population_context": [
                    {"code": code, **value}
                    for code, value in excluded.items()
                ],
            }
        )
    return {
        "status": gnomad_data.get("status"),
        "policy_id": gnomad_data.get("policy_id"),
        "policy_source": policy.get("source", ""),
        "policy_source_url": policy.get("source_url", ""),
        "scoring_rule": policy.get("scoring_rule", ""),
        "scored_non_founder_population_codes": scored_codes,
        "excluded_populations_are_context_only": True,
        "founder_context_only_observed": any(
            item.get("status") == "absent_in_non_founder_populations"
            for item in gnomad_data.get("datasets", {}).values()
        ),
        "datasets": datasets,
    }


def _scored_frequency_label(gnomad_data: Dict[str, Any]) -> str:
    maximum = _as_float(gnomad_data.get("max_af"))
    if maximum is None:
        return "non-cancer FAF95"
    matches = []
    for dataset_key, result in gnomad_data.get("datasets", {}).items():
        if _as_float(result.get("max_af")) != maximum:
            continue
        config = GNOMAD_LOCAL_DATASET_CONFIG.get(dataset_key, {})
        callset = result.get(config.get("callset", "")) or {}
        population = str(callset.get("faf95_pop") or "unknown").upper()
        label = result.get("label") or config.get("label", dataset_key)
        matches.append(f"{population}, {label}")
    location = "; ".join(matches) or "ENIGMA non-founder populations"
    return f"non-cancer FAF95 in {location}"


def query_gnomad_dataset_local(
    variant_id: Optional[str],
    coords: Optional[Any],
    config: Dict[str, Any],
    coverage_threshold: float,
    scored_ancestries: List[str],
) -> Dict[str, Any]:
    """Lookup one logical gnomAD data source from the local panel cache."""
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
        "quality_filter_passed": None,
        "filtered_records": [],
        "coverage": None,
        "errors": [],
    }

    if not variant_id or not coords:
        return result

    if not GNOMAD_CACHE:
        result["status"] = "cache_missing"
        result["errors"].append("local gnomAD cache not loaded")
        return result

    if GNOMAD_CACHE_MODE != "approved_snapshot":
        result["status"] = "cache_untrusted"
        result["errors"].append(
            f"gnomAD cache mode {GNOMAD_CACHE_MODE!r} is not approved for classification"
        )
        return result

    if not _coords_in_cached_region(coords, config["assembly"]):
        result["status"] = "outside_cached_region"
        result["errors"].append("variant coordinates outside cached panel regions")
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
        result["database_record_found"] = True
        result["dataset"] = dataset_records[0].get("dataset")

        result["filtered_records"] = [
            {
                "variant_id": rec.get("variant_id"),
                "filter": rec.get("filter"),
            }
            for rec in dataset_records
            if rec.get("filter") not in (None, ".", "PASS")
        ]
        passing_records = [
            rec for rec in dataset_records
            if rec.get("filter") in (None, ".", "PASS")
        ]
        result["quality_filter_passed"] = bool(passing_records)
        passing_non_founder_records = [
            rec for rec in passing_records if rec.get("non_founder_observed") is True
        ]
        if not passing_records:
            result["status"] = "filtered_record"
            result["found"] = None
        elif passing_non_founder_records:
            result["status"] = "found"
            result["found"] = True
        else:
            result["status"] = "absent_in_non_founder_populations"
            result["found"] = False

        freqs = []
        for rec in passing_records:
            freq_value, metric, population = _record_frequency_value(
                rec, scored_ancestries
            )
            if freq_value is not None:
                freqs.append((freq_value, metric, population, rec))

        if freqs:
            freq_value, metric, population, best_rec = max(
                freqs, key=lambda x: x[0]
            )
            result["max_af"] = freq_value
            result["frequency_metric"] = metric
        else:
            best_rec = passing_records[0] if passing_records else dataset_records[0]

        callset_summary = _record_to_callset_summary(
            best_rec, scored_ancestries
        )
        result[config["callset"]] = callset_summary
        result["non_founder_allele_count"] = callset_summary[
            "non_founder_allele_count"
        ]

        # Coverage always comes from the complete standalone positional cache.
        # A record-embedded single-locus value cannot prove a multi-base span.
        coverage = _lookup_coverage_by_position(
            coords,
            config["coverage_dataset_key"],
            config["assembly"],
            coverage_threshold,
        )
    else:
        result["status"] = "absent"
        result["found"] = False
        result["dataset"] = config["dataset_names"][0]
        coverage = _lookup_coverage_by_position(
            coords,
            config["coverage_dataset_key"],
            config["assembly"],
            coverage_threshold,
        )

    result["coverage"] = coverage
    return result


def _aggregate_coverage_from_dataset_results(
    dataset_results: Dict[str, Dict[str, Any]],
    required_dataset_keys: List[str],
    threshold: float,
) -> Dict[str, Any]:
    datasets = {}
    for key, ds in dataset_results.items():
        cov = ds.get("coverage") or {}
        mean_depth = _as_float(cov.get("mean_depth"))
        dataset_threshold = _as_float(cov.get("threshold"))
        passes = (
            bool(cov.get("passes"))
            if cov.get("passes") is not None
            else (
                mean_depth is not None
                and dataset_threshold is not None
                and mean_depth >= dataset_threshold
            )
        )
        datasets[key] = {
            "available": mean_depth is not None,
            "passes": passes,
            "mean_depth": mean_depth,
            "threshold": dataset_threshold,
            "source": cov.get("source"),
            "position_key": cov.get("position_key"),
            "position_keys": cov.get("position_keys") or [],
            "coverage_scope": cov.get("coverage_scope"),
            "span_start": cov.get("span_start"),
            "span_end": cov.get("span_end"),
            "positions_expected": cov.get("positions_expected"),
            "positions_available": cov.get("positions_available"),
            "missing_positions": cov.get("missing_positions") or [],
            "callsets": {
                GNOMAD_LOCAL_DATASET_CONFIG[key]["callset"]: {
                    "mean_depth": mean_depth,
                    "passes": passes,
                    "source": cov.get("source"),
                }
            },
        }

    required = list(required_dataset_keys)
    all_available = all(datasets.get(k, {}).get("available") for k in required)
    all_pass = all(datasets.get(k, {}).get("passes") for k in required)

    return {
        "status": "ok" if all_pass else ("missing" if not all_available else "insufficient"),
        "passes_pm2": all_pass,
        "min_required_mean_depth": threshold,
        "datasets": datasets,
    }


def get_gnomad_frequencies(
    gene: Optional[str],
    grch37: Optional[Any] = None,
    grch38: Optional[Any] = None,
) -> Dict[str, Any]:
    """Lookup gnomAD frequency and coverage from the local panel cache."""
    result = {
        "status": "not_queried",
        "gene": gene,
        "policy_id": None,
        "frequency_policy": {},
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

    policy = _classification_policy_for_gene(gene)
    if policy is None:
        result["status"] = "policy_unavailable"
        result["errors"].append(
            f"no active gene-specific gnomAD classification policy for {gene or 'unspecified gene'}"
        )
        return result
    frequency_policy = policy.get("frequency_criteria") or {}
    pm2_policy = frequency_policy.get("pm2") or {}
    required = list(
        pm2_policy.get("required_absence_dataset_runtime_keys") or []
    )
    scored_ancestries = list(
        frequency_policy.get("scored_non_founder_ancestries") or []
    )
    pm2_threshold = _as_float(pm2_policy.get("minimum_mean_depth"))
    if (
        not required
        or not scored_ancestries
        or pm2_threshold is None
        or any(key not in GNOMAD_LOCAL_DATASET_CONFIG for key in required)
    ):
        result["status"] = "policy_unavailable"
        result["errors"].append(
            f"active gnomAD policy {policy.get('policy_id')} is incomplete"
        )
        return result
    result["policy_id"] = policy["policy_id"]
    result["frequency_policy"] = frequency_policy

    coords_by_kind = {"grch37": grch37, "grch38": grch38}

    for dataset_key in required:
        config = GNOMAD_LOCAL_DATASET_CONFIG[dataset_key]
        coords = coords_by_kind.get(config["coord"])
        variant_id = _variant_id_from_coords(coords, with_chr=None)
        dataset_result = query_gnomad_dataset_local(
            variant_id,
            coords,
            config,
            pm2_threshold,
            scored_ancestries,
        )
        result["datasets"][dataset_key] = dataset_result
        if dataset_result.get("errors"):
            result["errors"].extend(dataset_result["errors"])

    result["coverage"] = _aggregate_coverage_from_dataset_results(
        result["datasets"], required, pm2_threshold
    )
    required_statuses = [result["datasets"].get(k, {}).get("status") for k in required]

    any_found = any(
        result["datasets"].get(k, {}).get("status") == "found"
        for k in required
    )
    all_absent = all(
        s in {"absent", "absent_in_non_founder_populations"}
        for s in required_statuses
    )
    any_cache_missing = any(
        s in ("cache_missing", "cache_untrusted", "dataset_not_available")
        for s in required_statuses
    )
    any_no_coords = any(s == "no_coordinates" for s in required_statuses)
    any_outside = any(s == "outside_cached_region" for s in required_statuses)

    freqs = []
    for key in required:
        ds = result["datasets"].get(key, {})
        af = _as_float(ds.get("max_af"))
        if af is not None and ds.get("frequency_metric") == "faf95":
            freqs.append((af, ds.get("frequency_metric")))

    if freqs:
        af, metric = max(freqs, key=lambda x: x[0])
        result["max_af"] = af
        result["frequency_metric"] = metric

    result["found"] = any_found
    result["pm2_coverage_ok"] = result["coverage"]["passes_pm2"]

    # Both required datasets must establish absence. Missing v3 data is not
    # evidence of absence and must never be converted into PM2.
    v2_status = result["datasets"].get("v2_1_non_cancer", {}).get("status", "")
    v3_status = result["datasets"].get("v3_1_non_cancer", {}).get("status", "")
    v3_not_in_cache = v3_status in (
        "outside_cached_region", "dataset_not_available", "cache_missing"
    )

    if all_absent and result["pm2_coverage_ok"]:
        result["pm2_absence_established"] = True
        result["pm2_datasets_note"] = " + ".join(
            GNOMAD_LOCAL_DATASET_CONFIG[key]["label"] for key in required
        )
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

    result["population_frequency_audit"] = _population_frequency_audit(result)
    result["founder_context_only_observed"] = result[
        "population_frequency_audit"
    ]["founder_context_only_observed"]
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


def _frequency_depth_ok(
    gnomad_data: Dict[str, Any], minimum_mean_depth: float
) -> bool:
    """Check policy depth in a dataset contributing the maximum frequency."""
    max_af = _as_float(gnomad_data.get("max_af"))
    for dataset in gnomad_data.get("datasets", {}).values():
        if dataset.get("status") != "found":
            continue
        dataset_af = _as_float(dataset.get("max_af"))
        if max_af is not None and dataset_af != max_af:
            continue
        coverage = dataset.get("coverage") or {}
        mean_depth = _as_float(coverage.get("mean_depth"))
        if mean_depth is not None and mean_depth >= minimum_mean_depth:
            return True
    return False


def _frequency_qc_ok(gnomad_data: Dict[str, Any]) -> bool:
    """Require a PASS record in each dataset contributing the scored FAF95."""
    max_af = _as_float(gnomad_data.get("max_af"))
    for dataset in gnomad_data.get("datasets", {}).values():
        if dataset.get("status") != "found":
            continue
        dataset_af = _as_float(dataset.get("max_af"))
        if max_af is not None and dataset_af != max_af:
            continue
        if dataset.get("quality_filter_passed") is True:
            return True
        for callset_name in ("exomes", "genomes"):
            callset = dataset.get(callset_name) or {}
            if callset.get("available") and not callset.get("filters"):
                return True
    return False


def evaluate_frequency_criteria(
    gnomad_data: Dict[str, Any],
    variant_type: str,
    gene: Optional[str] = None,
    c_notation: Optional[str] = None,
) -> Dict:
    """Evaluate BA1, BS1, and PM2 from local gnomAD data."""
    criteria = {}
    policy = _classification_policy_for_gene(gene)
    if policy is None or policy.get("policy_id") != gnomad_data.get("policy_id"):
        return {
            "_gnomad_info": {
                "applies": False,
                "reason": (
                    "Frequency criteria not applied: no matching active "
                    f"gene-specific gnomAD policy for {gene or 'unspecified gene'}"
                ),
            }
        }
    frequency_policy = policy.get("frequency_criteria") or {}
    ba1_policy = frequency_policy.get("ba1") or {}
    bs1_policy = frequency_policy.get("bs1") or {}
    bs1_strong_policy = bs1_policy.get("strong") or {}
    bs1_supporting_policy = bs1_policy.get("supporting") or {}
    pm2_policy = frequency_policy.get("pm2") or {}
    observation_policy = frequency_policy.get("outbred_observation_count") or {}
    ba1_threshold = _as_float(ba1_policy.get("threshold"))
    bs1_strong_threshold = _as_float(bs1_strong_policy.get("threshold"))
    bs1_supporting_threshold = _as_float(
        bs1_supporting_policy.get("lower_threshold")
    )
    frequency_depth_threshold = _as_float(
        ba1_policy.get("minimum_mean_depth")
    )
    pm2_depth_threshold = _as_float(pm2_policy.get("minimum_mean_depth"))
    minimum_ba1_bs1_observations = _as_int(
        observation_policy.get("minimum_observations_for_ba1_bs1")
    )
    required_values = (
        ba1_threshold,
        bs1_strong_threshold,
        bs1_supporting_threshold,
        frequency_depth_threshold,
        pm2_depth_threshold,
        minimum_ba1_bs1_observations,
    )
    if any(value is None for value in required_values):
        return {
            "_gnomad_info": {
                "applies": False,
                "reason": "Frequency criteria not applied: active policy is incomplete",
            }
        }
    pm2_excluded_types = {
        str(value).lower() for value in pm2_policy.get("excluded_variant_types", [])
    }
    # The protein consequence can route an indel into the PTC/nonsense branch.
    # For example, c.5533_5534insG has p.(Tyr1845Ter), so its consequence type
    # is "nonsense" even though the underlying allele is an insertion.  PM2
    # applicability is defined by the allele class, therefore the c. HGVS
    # operation must independently exclude indels and exon CNVs.
    c_allele_is_indel = bool(
        re.search(r"(?:delins|del|dup|ins)", (c_notation or "").lower())
    )
    pm2_excluded = (
        variant_type.lower() in pm2_excluded_types
        or c_allele_is_indel
    )

    status = gnomad_data.get("status", "not_queried")
    max_af = _as_float(gnomad_data.get("max_af"))
    metric = gnomad_data.get("frequency_metric") or "frequency"

    # Frequency too common -> benign evidence. ENIGMA requires filtering AF;
    # raw AF and popmax AF are displayed for context but are not substitutes.
    if max_af is not None:
        if metric != "faf95":
            criteria["_gnomad_info"] = {
                "applies": False,
                "reason": (
                    "BA1/BS1 not applied: ENIGMA-compatible non-cancer FAF95 "
                    f"is unavailable; {metric} cannot be used as a fallback"
                ),
            }
            max_af = None

    if max_af is not None:
        af_pct = f"{max_af * 100:.6g}%"
        metric_note = _scored_frequency_label(gnomad_data)
        founder_snapshot_note = ""
        population_policy_note = (
            "; founder and other non-scoring population groups excluded per "
            f"{frequency_policy.get('source') or 'the active gene-specific policy'}"
        )

        contributing_observation_counts = [
            _as_int(dataset.get("non_founder_allele_count"))
            for dataset in gnomad_data.get("datasets", {}).values()
            if dataset.get("status") == "found"
            and _as_float(dataset.get("max_af")) == max_af
        ]
        if max_af > bs1_supporting_threshold and (
            not contributing_observation_counts
            or any(count is None for count in contributing_observation_counts)
        ):
            criteria["_gnomad_info"] = {
                "applies": False,
                "reason": (
                    "BA1/BS1 not applied: the outbred observation count is "
                    "missing from the approved gnomAD result"
                ),
            }
            return criteria
        if max_af > bs1_supporting_threshold and all(
            count < minimum_ba1_bs1_observations
            for count in contributing_observation_counts
        ):
            criteria["_gnomad_info"] = {
                "applies": False,
                "reason": (
                    "BA1/BS1 not applied: a single observation in an ENIGMA "
                    "outbred population is not informative"
                ),
            }
            return criteria

        if max_af > bs1_supporting_threshold:
            if not _frequency_qc_ok(gnomad_data):
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": "BA1/BS1 not applied: the contributing gnomAD record did not pass dataset QC filters",
                }
                return criteria

            if not _frequency_depth_ok(
                gnomad_data, frequency_depth_threshold
            ):
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": (
                        "BA1/BS1 not applied: active policy requires mean read "
                        f"depth >= {frequency_depth_threshold:g}"
                    ),
                }
                return criteria

            from backend.lookups.founder_variants import lookup_pathogenic_founder_variant

            founder = lookup_pathogenic_founder_variant(gene or "", c_notation or "")
            if founder.get("is_pathogenic_founder") is True:
                if max_af > ba1_threshold:
                    excluded_code = "BA1"
                    excluded_policy = ba1_policy
                    threshold = ba1_threshold
                elif max_af > bs1_strong_threshold:
                    excluded_code = "BS1_Strong"
                    excluded_policy = bs1_strong_policy
                    threshold = bs1_strong_threshold
                else:
                    excluded_code = "BS1_Supporting"
                    excluded_policy = bs1_supporting_policy
                    threshold = bs1_supporting_threshold
                exclusion_reason = (
                    f"gnomAD {metric_note} {af_pct} exceeds the "
                    f"{excluded_code.replace('_', ' ')} threshold "
                    f"{threshold * 100:g}%, but {excluded_code.split('_')[0]} was "
                    "not applied and added no points: ENIGMA v1.2 excludes "
                    "well-established pathogenic founder variants; "
                    f"{founder.get('reason')}"
                )
                criteria["_excluded_criteria"] = {
                    excluded_code: {
                        "applies": False,
                        "strength": excluded_policy.get("strength"),
                        "points": 0,
                        "reason": exclusion_reason,
                        "source": frequency_policy.get("source_url", ""),
                    }
                }
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": exclusion_reason,
                    "founder_exception": founder,
                }
                return criteria
            if founder.get("status") == "unavailable":
                criteria["_gnomad_info"] = {
                    "applies": False,
                    "reason": (
                        "BA1/BS1 not applied: the pathogenic founder exception "
                        f"could not be checked; {founder.get('reason')}"
                    ),
                    "founder_exception": founder,
                }
                return criteria
            founder_snapshot_note = (
                f"; pathogenic-founder exception checked against snapshot "
                f"{founder.get('snapshot_version') or 'unknown'}"
            )
            if max_af > ba1_threshold:
                criteria["BA1"] = {
                    "applies": True,
                    "strength": ba1_policy["strength"],
                    "points": ba1_policy["points"],
                    "reason": (
                        f"gnomAD {metric_note} {af_pct} > {ba1_threshold * 100:g}% "
                        "- Stand-alone Benign"
                        f"{population_policy_note}{founder_snapshot_note}"
                    )
                }
                return criteria
            if max_af > bs1_strong_threshold:
                criteria["BS1_Strong"] = {
                    "applies": True,
                    "strength": bs1_strong_policy["strength"],
                    "points": bs1_strong_policy["points"],
                    "reason": (
                        f"gnomAD {metric_note} {af_pct} > "
                        f"{bs1_strong_threshold * 100:g}%"
                        f"{population_policy_note}{founder_snapshot_note}"
                    )
                }
                return criteria
            criteria["BS1_Supporting"] = {
                "applies": True,
                "strength": bs1_supporting_policy["strength"],
                "points": bs1_supporting_policy["points"],
                "reason": (
                    f"gnomAD {metric_note} {af_pct} > "
                    f"{bs1_supporting_threshold * 100:g}%"
                    f"{population_policy_note}{founder_snapshot_note}"
                )
            }
            return criteria
        # Present but too rare is not PM2.
        if gnomad_data.get("found"):
            return criteria

    if gnomad_data.get("found") and max_af is None:
        failed_qc = any(
            dataset.get("status") == "found"
            and dataset.get("quality_filter_passed") is False
            for dataset in gnomad_data.get("datasets", {}).values()
        )
        criteria.setdefault(
            "_gnomad_info",
            {
                "applies": False,
                "reason": (
                    "BA1/BS1 not applied: the variant is present in gnomAD, but "
                    + (
                        "the record failed dataset QC filters"
                        if failed_qc
                        else "ENIGMA-compatible non-cancer FAF95 is unavailable"
                    )
                ),
            },
        )

    if pm2_excluded:
        exclusion_basis = (
            f"c. HGVS {c_notation} describes an indel"
            if c_allele_is_indel
            else f"variant type {variant_type} is excluded"
        )
        criteria["PM2"] = {
            "applies": False, "strength": None, "points": 0,
            "reason": (
                f"PM2 not applicable because {exclusion_basis} under "
                f"policy {policy['policy_id']}"
            )
        }
        return criteria

    if gnomad_data.get("pm2_absence_established"):
        datasets_note = gnomad_data.get("pm2_datasets_note", "v2.1.1 + v3.1.2")
        founder_only_note = (
            "; observations confined to excluded founder/non-scoring populations "
            "were not treated as presence in an ENIGMA outbred population"
            if gnomad_data.get("founder_context_only_observed")
            else ""
        )
        criteria["PM2_Supporting"] = {
            "applies": True,
            "strength": pm2_policy["strength"],
            "points": pm2_policy["points"],
            "reason": (
                f"Absent from gnomAD {datasets_note} non-cancer callset(s) "
                "in all policy-defined non-founder populations with mean depth "
                f">= {pm2_depth_threshold:g} "
                f"across the variant reference span{founder_only_note}"
            )
        }
        return criteria

    reason_by_status = {
        "cache_missing": "local gnomAD cache missing or incomplete - PM2 not applied",
        "cache_untrusted": "local gnomAD cache is not an approved snapshot - frequency criteria not applied",
        "partial": "local gnomAD lookup partial - PM2 not applied",
        "no_coordinates": "No genomic coordinates for required gnomAD lookup - PM2 not applied",
        "outside_cached_region": "Variant outside cached panel gnomAD regions - PM2 not applied",
        "absent_without_sufficient_coverage": (
            "Absent from local gnomAD cache but coverage mean depth is below "
            f"{pm2_depth_threshold:g} or missing - PM2 not applied"
        ),
        "not_queried": "gnomAD not queried - PM2 not applied",
        "absent_v2_only": "gnomAD v2.1.1 absence confirmed but v3.1.2 coverage insufficient - PM2 not applied",
        "policy_unavailable": "Gene-specific gnomAD policy is unavailable - frequency criteria not applied",
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
