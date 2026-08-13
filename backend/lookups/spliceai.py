# ============================================================
# SpliceAI lookup via Broad API and profile-pinned runtime cache
#
# Until the Appendix J reference caches are completely rebuilt, production
# uses the Broad Institute SpliceAI Lookup API and its profile-pinned runtime
# cache. Immutable precomputed caches are opt-in and remain disabled by
# default.
#
# API endpoint: https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/
# Variant format: chr{chrom}-{pos}-{ref}-{alt}
# Rate limit: a few requests per minute - cache prevents repeated calls.
#
# The previous MANE VCF subset approach was removed because the Ensembl MANE v1.0
# file uses an older Gencode version and gives incorrect scores for some variants
# (e.g. BRCA1 c.4185G>A: MANE gives DS_DL=0.01, Broad API gives DS_DL=0.93).
# ============================================================
from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse
import gzip
import json as _json
import os
import shutil
import subprocess
import hashlib
import tempfile
from backend.data_health import clear_issue, register_issue
from backend.spliceai_profile import (
    SPLICEAI_AGGREGATION,
    SPLICEAI_ALTERNATE_FIELDS,
    SPLICEAI_ANNOTATION_SUBSET,
    SPLICEAI_DELTA_FIELDS,
    SPLICEAI_GENOME_ASSEMBLY,
    SPLICEAI_HIGH_THRESHOLD,
    SPLICEAI_LOW_THRESHOLD,
    SPLICEAI_MASK,
    SPLICEAI_MAX_DISTANCE,
    SPLICEAI_PROFILE,
    SPLICEAI_PROFILE_ID,
    SPLICEAI_PROFILE_SHA256,
    SPLICEAI_REFERENCE_FIELDS,
    SPLICEAI_TRANSCRIPT_POLICY_REQUIRED,
    validate_scoring_metadata,
)

from backend.lookups.coordinates import resolve_variant, get_grch38


def choose_project_root() -> Path:
    env_root = os.environ.get("BRCA_ACMG_PROJECT_ROOT")
    if env_root:
        return Path(env_root)
    # The repository root is stable regardless of the process working directory.
    # Deployments with external data must opt in through BRCA_ACMG_PROJECT_ROOT.
    return Path(__file__).resolve().parents[2]

PROJECT_ROOT  = choose_project_root()
SPLICEAI_DIR  = PROJECT_ROOT / "data" / "spliceai"
SPLICEAI_DIR.mkdir(parents=True, exist_ok=True)

def choose_runtime_cache_dir() -> Path:
    """Choose writable runtime storage without mixing it with snapshots."""
    configured = os.environ.get("ARIANE_RUNTIME_CACHE_DIR")
    if configured:
        return Path(configured)
    railway_volume = os.environ.get("RAILWAY_VOLUME_MOUNT_PATH")
    if railway_volume:
        return Path(railway_volume) / "ariane-runtime-cache"
    # Local-development compatibility. Production deployments should provide
    # ARIANE_RUNTIME_CACHE_DIR or attach a Railway volume.
    return SPLICEAI_DIR


RUNTIME_CACHE_DIR = choose_runtime_cache_dir()
RUNTIME_CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Mutable API results are separate from immutable, versioned snapshots.
SPLICEAI_API_CACHE_PATH = RUNTIME_CACHE_DIR / "spliceai_api_cache.json"
SPLICEAI_PRECOMPUTED_CACHE_PATH = Path(os.environ.get(
    "SPLICEAI_PRECOMPUTED_CACHE_PATH",
    SPLICEAI_DIR / "spliceai_brca_snv_reference_cache.json",
))
SPLICEAI_INTRONIC_CACHE_PATH = SPLICEAI_DIR / "spliceai_brca_intronic_snv_reference_cache.json"

# In-memory caches
SPLICEAI_CACHE:        Dict[str, float] = {}   # policy:gene:c_notation -> score
SPLICEAI_STATUS_CACHE: Dict[str, dict]  = {}   # gene:c_notation -> status details
SPLICEAI_PRECOMPUTED_CACHE: Optional[Dict[str, dict]] = None

# Broad API endpoint (Google Cloud Run, hg38) or a local compatible server.
DEFAULT_SPLICEAI_API_URL = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/"


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


def _normalize_api_url(value: str) -> str:
    value = (value or DEFAULT_SPLICEAI_API_URL).strip()
    if not value.endswith("/"):
        value += "/"
    return value


SPLICEAI_API_URL = _normalize_api_url(os.environ.get("SPLICEAI_API_URL", DEFAULT_SPLICEAI_API_URL))
# Interactive requests are proxied by nginx with a 60-second read timeout. A
# non-critical SpliceAI lookup must fail closed early enough for the classifier
# to return a result with an explicit unavailable-source warning. Offline cache
# builders do not use this web-request deadline.
SPLICEAI_API_TIMEOUT = _env_int("SPLICEAI_API_TIMEOUT", 25)
SPLICEAI_API_RATE_SLEEP = _env_float("SPLICEAI_API_RATE_SLEEP", 1.5)
SPLICEAI_API_SOURCE = os.environ.get(
    "SPLICEAI_API_SOURCE",
    "Local Broad SpliceAI API" if "localhost" in SPLICEAI_API_URL else "Broad SpliceAI API",
)

REFERENCE_TRANSCRIPTS = SPLICEAI_PROFILE["reference_transcripts"]

_requested_transcript_policy = os.environ.get(
    "SPLICEAI_TRANSCRIPT_POLICY", SPLICEAI_TRANSCRIPT_POLICY_REQUIRED
).strip().lower()
if _requested_transcript_policy != SPLICEAI_TRANSCRIPT_POLICY_REQUIRED:
    raise RuntimeError(
        "SPLICEAI_TRANSCRIPT_POLICY conflicts with the ENIGMA v1.2 scoring "
        f"profile: {_requested_transcript_policy!r}; expected "
        f"{SPLICEAI_TRANSCRIPT_POLICY_REQUIRED!r}"
    )
SPLICEAI_TRANSCRIPT_POLICY = SPLICEAI_TRANSCRIPT_POLICY_REQUIRED

SPLICEAI_USE_PRECOMPUTED_CACHE = os.environ.get(
    "SPLICEAI_USE_PRECOMPUTED_CACHE",
    "0",
).strip().lower() not in {"0", "false", "no", "off"}

def _load_api_cache() -> dict:
    """Load persisted API cache from Drive."""
    if SPLICEAI_API_CACHE_PATH.exists():
        try:
            with open(SPLICEAI_API_CACHE_PATH) as f:
                result = _json.load(f)
            clear_issue("SpliceAI API cache")
            return result
        except Exception as exc:
            register_issue(
                "SpliceAI API cache",
                f"could not load {SPLICEAI_API_CACHE_PATH}: {type(exc).__name__}: {exc}",
            )
    return {}


def _load_precomputed_cache() -> dict:
    """Load the local BRCA SNV reference-transcript SpliceAI cache once."""
    global SPLICEAI_PRECOMPUTED_CACHE
    if SPLICEAI_PRECOMPUTED_CACHE is not None:
        return SPLICEAI_PRECOMPUTED_CACHE

    SPLICEAI_PRECOMPUTED_CACHE = {}
    if not SPLICEAI_USE_PRECOMPUTED_CACHE:
        return SPLICEAI_PRECOMPUTED_CACHE
    if SPLICEAI_TRANSCRIPT_POLICY != "reference_transcript":
        return SPLICEAI_PRECOMPUTED_CACHE
    for path in (SPLICEAI_PRECOMPUTED_CACHE_PATH, SPLICEAI_INTRONIC_CACHE_PATH):
        component = (
            "SpliceAI intronic cache"
            if path == SPLICEAI_INTRONIC_CACHE_PATH
            else "SpliceAI coding SNV cache"
        )
        if not path.exists():
            register_issue(component, f"cache is missing: {path}")
            continue
        try:
            with open(path, encoding="utf-8") as handle:
                raw = _json.load(handle)
            metadata_path = path.with_name(path.stem + ".metadata.json")
            if not metadata_path.is_file():
                register_issue(
                    component,
                    f"cache build is incomplete: metadata is missing: {metadata_path}",
                )
                continue
            with metadata_path.open(encoding="utf-8") as handle:
                metadata = _json.load(handle)
            if not isinstance(raw, dict) or not isinstance(metadata, dict):
                register_issue(component, "cache or metadata is not a JSON object")
                continue
            profile_errors = validate_scoring_metadata(metadata)
            if profile_errors:
                register_issue(
                    component,
                    "cache scoring profile is not ENIGMA Appendix J compatible: "
                    + "; ".join(profile_errors),
                )
                continue
            actual_sha = hashlib.sha256(path.read_bytes()).hexdigest()
            expected_count = metadata.get("expected_variants")
            if expected_count is None:
                expected_count = metadata.get("coordinate_variants")
            if (
                not isinstance(expected_count, int)
                or metadata.get("status_ok") != expected_count
                or metadata.get("status_error", 0) != 0
                or metadata.get("cache_entries") != expected_count
                or len(raw) != expected_count
                or metadata.get("sha256", "").lower() != actual_sha.lower()
            ):
                register_issue(
                    component,
                    "cache build is incomplete or checksum/count validation failed",
                )
                continue
            invalid_entries = sum(
                1 for entry in raw.values() if not _precomputed_entry_is_complete(entry)
            )
            if invalid_entries:
                register_issue(
                    component,
                    f"cache has {invalid_entries} entries without complete delta, REF, and ALT audit data",
                )
                continue
            SPLICEAI_PRECOMPUTED_CACHE.update(raw)
            clear_issue(component)
        except Exception as exc:
            print(f"Warning: could not load precomputed SpliceAI cache {path}: {exc}")
            register_issue(
                component,
                f"could not load {path}: {type(exc).__name__}: {exc}",
            )
    return SPLICEAI_PRECOMPUTED_CACHE


def _save_api_cache(cache: dict) -> bool:
    """Atomically persist the API cache and report whether it succeeded."""
    temporary_path = None
    try:
        RUNTIME_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=RUNTIME_CACHE_DIR,
            prefix="spliceai_api_cache.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            _json.dump(cache, handle, indent=2)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, SPLICEAI_API_CACHE_PATH)
        clear_issue("SpliceAI API cache")
        return True
    except Exception as e:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError:
                pass
        print(f"Warning: could not save SpliceAI cache: {e}")
        register_issue(
            "SpliceAI API cache",
            f"score was obtained and used, but the runtime cache could not be saved to "
            f"{SPLICEAI_API_CACHE_PATH}; this request is unaffected, but the score may need "
            f"to be fetched again after restart: {type(e).__name__}: {e}",
        )
        return False


def _cache_key(gene: str, c_notation: str) -> str:
    return f"{SPLICEAI_PROFILE_ID}:{SPLICEAI_TRANSCRIPT_POLICY}:{gene}:{c_notation}"


def _precomputed_cache_keys(gene: str, c_notation: str) -> tuple[str, str]:
    raw_key = f"{gene}:{c_notation}"
    return raw_key, f"reference_transcript:{raw_key}"


def _float_score_map(value: object, fields: tuple[str, ...]) -> Optional[dict[str, float]]:
    if not isinstance(value, dict) or set(value) != set(fields):
        return None
    result: dict[str, float] = {}
    try:
        for field in fields:
            score = float(value[field])
            if not 0.0 <= score <= 1.0:
                return None
            result[field] = score
    except (TypeError, ValueError):
        return None
    return result


def _precomputed_entry_is_complete(entry: object) -> bool:
    if not isinstance(entry, dict) or entry.get("status") != "ok":
        return False
    delta_scores = _float_score_map(entry.get("delta_scores"), SPLICEAI_DELTA_FIELDS)
    reference_scores = _float_score_map(
        entry.get("reference_scores"), SPLICEAI_REFERENCE_FIELDS
    )
    alternate_scores = _float_score_map(
        entry.get("alternate_scores"), SPLICEAI_ALTERNATE_FIELDS
    )
    if delta_scores is None or reference_scores is None or alternate_scores is None:
        return False
    max_field = entry.get("max_delta_field")
    if max_field not in SPLICEAI_DELTA_FIELDS:
        return False
    try:
        score = float(entry.get("score"))
    except (TypeError, ValueError):
        return False
    return score == max(delta_scores.values()) and score == delta_scores[max_field]


def _runtime_entry_matches_profile(entry: object) -> bool:
    if not isinstance(entry, dict):
        return False
    if any(
        (
            entry.get("scoring_profile_id") != SPLICEAI_PROFILE_ID,
            entry.get("scoring_profile_sha256") != SPLICEAI_PROFILE_SHA256,
            entry.get("genome_assembly") != SPLICEAI_GENOME_ASSEMBLY,
            entry.get("distance") != SPLICEAI_MAX_DISTANCE,
            entry.get("mask") != SPLICEAI_MASK,
            entry.get("annotation_subset") != SPLICEAI_ANNOTATION_SUBSET,
            entry.get("aggregation") != SPLICEAI_AGGREGATION,
            entry.get("transcript_policy") != SPLICEAI_TRANSCRIPT_POLICY,
        )
    ):
        return False
    return _precomputed_entry_is_complete({**entry, "status": "ok"})


def _entry_score(entry: dict) -> Optional[float]:
    try:
        score = entry.get("score")
        return float(score) if score is not None else None
    except (TypeError, ValueError):
        return None


def _lookup_precomputed_score(gene: str, c_notation: str) -> Optional[dict]:
    if SPLICEAI_TRANSCRIPT_POLICY != "reference_transcript":
        return None
    cache = _load_precomputed_cache()
    for key in _precomputed_cache_keys(gene, c_notation):
        entry = cache.get(key)
        if not isinstance(entry, dict):
            continue
        if not _precomputed_entry_is_complete(entry):
            continue
        score = _entry_score(entry)
        if score is None:
            continue
        return {
            "score": score,
            "max_delta_field": entry.get("max_delta_field", ""),
            "delta_scores": entry.get("delta_scores", {}),
            "reference_scores": entry.get("reference_scores", {}),
            "alternate_scores": entry.get("alternate_scores", {}),
            "selected_transcript": REFERENCE_TRANSCRIPTS.get(gene, {}).get("ensembl", ""),
            "reference_transcript_score": score,
            "max_any_transcript_score": entry.get("max_any_transcript_score"),
            "max_any_transcript": entry.get("max_any_transcript", ""),
            "source": entry.get("source") or "precomputed BRCA SNV SpliceAI cache",
            "cache_key": key,
            "grch38": entry.get("grch38", ""),
            "variant": entry.get("variant", ""),
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
            "distance": SPLICEAI_MAX_DISTANCE,
            "mask": SPLICEAI_MASK,
            "annotation_subset": SPLICEAI_ANNOTATION_SUBSET,
            "genome_assembly": SPLICEAI_GENOME_ASSEMBLY,
            "aggregation": SPLICEAI_AGGREGATION,
        }
    return None


def _score_row(row: dict) -> tuple[float, str]:
    best_score = 0.0
    best_field = ""
    for key in SPLICEAI_DELTA_FIELDS:
        try:
            value = float(row.get(key, 0) or 0)
        except (ValueError, TypeError):
            value = 0.0
        if value > best_score:
            best_score = value
            best_field = key
    return best_score, best_field


def _row_score_audit(row: dict) -> Optional[dict]:
    delta_scores = _float_score_map(
        {field: row.get(field) for field in SPLICEAI_DELTA_FIELDS},
        SPLICEAI_DELTA_FIELDS,
    )
    reference_scores = _float_score_map(
        {field: row.get(field) for field in SPLICEAI_REFERENCE_FIELDS},
        SPLICEAI_REFERENCE_FIELDS,
    )
    alternate_scores = _float_score_map(
        {field: row.get(field) for field in SPLICEAI_ALTERNATE_FIELDS},
        SPLICEAI_ALTERNATE_FIELDS,
    )
    if delta_scores is None or reference_scores is None or alternate_scores is None:
        return None
    max_field = max(SPLICEAI_DELTA_FIELDS, key=delta_scores.__getitem__)
    return {
        "score": delta_scores[max_field],
        "max_delta_field": max_field,
        "delta_scores": delta_scores,
        "reference_scores": reference_scores,
        "alternate_scores": alternate_scores,
    }


def _row_matches_reference_transcript(gene: str, row: dict) -> bool:
    reference = REFERENCE_TRANSCRIPTS.get(gene)
    if not reference:
        return False

    ensembl = reference["ensembl"]
    ensembl_no_version = ensembl.split(".")[0]
    t_id = str(row.get("t_id") or "")
    if t_id == ensembl or t_id.split(".")[0] == ensembl_no_version:
        return True

    refseq = reference["refseq"]
    refseq_no_version = refseq.split(".")[0]
    for item in row.get("t_refseq_ids") or []:
        item = str(item)
        if item == refseq or item.split(".")[0] == refseq_no_version:
            return True
    return False


def _reference_transcript_match_rank(gene: str, row: dict) -> int:
    reference = REFERENCE_TRANSCRIPTS.get(gene)
    if not reference:
        return 0
    t_id = str(row.get("t_id") or "")
    if t_id == reference["ensembl"]:
        return 4
    if t_id.split(".")[0] == reference["ensembl"].split(".")[0]:
        return 3
    refseq_ids = [str(value) for value in row.get("t_refseq_ids") or []]
    if reference["refseq"] in refseq_ids:
        return 2
    refseq_base = reference["refseq"].split(".")[0]
    if any(value.split(".")[0] == refseq_base for value in refseq_ids):
        return 1
    return 0


def _select_spliceai_score(gene: str, scores: list[dict]) -> dict:
    audited_rows: list[tuple[dict, dict]] = []
    for row in scores:
        audit = _row_score_audit(row)
        if audit is not None:
            audited_rows.append((row, audit))
    if not audited_rows:
        return {
            "score": None,
            "error": "SpliceAI response lacks complete delta, REF, or ALT score fields",
        }

    max_row, max_audit = max(audited_rows, key=lambda item: item[1]["score"])
    ranked_reference_rows = [
        (_reference_transcript_match_rank(gene, row), row, audit)
        for row, audit in audited_rows
        if _reference_transcript_match_rank(gene, row) > 0
    ]
    reference_row = None
    reference_audit = None
    if ranked_reference_rows:
        best_rank = max(item[0] for item in ranked_reference_rows)
        best_matches = [item for item in ranked_reference_rows if item[0] == best_rank]
        signatures = {
            json.dumps(item[2], sort_keys=True, separators=(",", ":"))
            for item in best_matches
        }
        if len(signatures) > 1:
            return {
                "score": None,
                "error": "Ambiguous SpliceAI records for the required reference transcript",
            }
        _, reference_row, reference_audit = best_matches[0]

    max_any_score = max_audit["score"]
    max_any_transcript = str(max_row.get("t_id") or "")
    ref_score = reference_audit["score"] if reference_audit else None
    ref_transcript = str(reference_row.get("t_id") or "") if reference_row else ""

    if SPLICEAI_TRANSCRIPT_POLICY == "max_any_transcript":
        return {
            "score": max_any_score,
            "max_delta_field": max_audit["max_delta_field"],
            "delta_scores": max_audit["delta_scores"],
            "reference_scores": max_audit["reference_scores"],
            "alternate_scores": max_audit["alternate_scores"],
            "selected_transcript": max_any_transcript,
            "selected_transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
            "reference_transcript_score": ref_score,
            "reference_transcript": ref_transcript,
            "max_any_transcript_score": max_any_score,
            "max_any_transcript": max_any_transcript,
        }

    return {
        "score": ref_score,
        "max_delta_field": reference_audit["max_delta_field"] if reference_audit else "",
        "delta_scores": reference_audit["delta_scores"] if reference_audit else {},
        "reference_scores": reference_audit["reference_scores"] if reference_audit else {},
        "alternate_scores": reference_audit["alternate_scores"] if reference_audit else {},
        "selected_transcript": ref_transcript,
        "selected_transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
        "reference_transcript_score": ref_score,
        "reference_transcript": ref_transcript,
        "max_any_transcript_score": max_any_score,
        "max_any_transcript": max_any_transcript,
    }


def _query_spliceai_api(gene: str, chrom: str, pos: int, ref: str, alt: str) -> Optional[dict]:
    """
    Query Broad SpliceAI API for a single variant.
    Returns selected SpliceAI score details, or None on failure.
    """
    chrom_clean = str(chrom).replace("chr", "")
    variant_str = f"chr{chrom_clean}-{pos}-{ref}-{alt}"
    query = {
        "variant": variant_str,
        "hg": 38,
        "distance": SPLICEAI_MAX_DISTANCE,
        "mask": SPLICEAI_MASK,
        "bc": SPLICEAI_ANNOTATION_SUBSET,
    }
    url = f"{SPLICEAI_API_URL}?{urllib.parse.urlencode(query)}"
    try:
        req = urllib.request.Request(
            url,
            headers={"Accept": "application/json", "User-Agent": "BRCA-ACMG-Module1/1.6.4"},
        )
        with urllib.request.urlopen(req, timeout=SPLICEAI_API_TIMEOUT) as resp:
            data = _json.loads(resp.read())

        response_assembly = str(data.get("genomeVersion") or data.get("hg") or "")
        try:
            response_distance = int(data.get("distance"))
            response_mask = int(data.get("mask"))
        except (TypeError, ValueError):
            return {
                "score": None,
                "error": "API response did not report verifiable distance and mask parameters",
            }
        if (
            response_assembly != "38"
            or response_distance != SPLICEAI_MAX_DISTANCE
            or response_mask != SPLICEAI_MASK
        ):
            return {
                "score": None,
                "error": (
                    "API response scoring profile mismatch: "
                    f"assembly={response_assembly!r}, distance={response_distance}, "
                    f"mask={response_mask}"
                ),
            }

        scores = data.get("scores", [])
        if not scores:
            return {"score": None, "error": "API response contained no transcript scores"}

        selected = _select_spliceai_score(gene, scores)
        if selected.get("score") is None:
            return selected
        selected["api_source"] = data.get("source") or SPLICEAI_API_SOURCE
        selected["n_transcript_scores"] = len(scores)
        selected.update({
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
            "genome_assembly": SPLICEAI_GENOME_ASSEMBLY,
            "distance": response_distance,
            "mask": response_mask,
            "annotation_subset": SPLICEAI_ANNOTATION_SUBSET,
            "aggregation": SPLICEAI_AGGREGATION,
        })
        return selected

    except Exception as e:
        return {"score": None, "error": f"{type(e).__name__}: {e}"}


def get_spliceai_score(gene: str, c_notation: str) -> Optional[float]:
    """
    Look up SpliceAI score through the profile-pinned API runtime path.

    Returns a float score or None. None means unavailable, not 0.0.
    Benign criteria must only use confirmed scores <= 0.1.

    Current API-primary lookup order:
      1. In-memory cache (fast, within session)
      2. Precomputed cache only when explicitly enabled after a complete rebuild
      3. Persistent Broad API runtime cache
      4. Broad SpliceAI API

    The precomputed step is disabled by default while Appendix J caches are
    being rebuilt, so legacy cache metadata do not degrade API-primary runs.
    """
    variant_key = f"{gene}:{c_notation}"
    cache_key = _cache_key(gene, c_notation)

    # 1. in-memory cache
    if cache_key in SPLICEAI_CACHE:
        return SPLICEAI_CACHE[cache_key]

    # 2. precomputed BRCA coding SNV cache
    precomputed = _lookup_precomputed_score(gene, c_notation)
    if precomputed is not None:
        score = precomputed["score"]
        SPLICEAI_CACHE[cache_key] = score
        SPLICEAI_STATUS_CACHE[variant_key] = {
            "status": "ok",
            "score": score,
            "reason": "Loaded from local precomputed BRCA SNV SpliceAI cache",
            "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
            "selected_transcript": precomputed.get("selected_transcript"),
            "reference_transcript_score": precomputed.get("reference_transcript_score"),
            "max_any_transcript_score": precomputed.get("max_any_transcript_score"),
            "max_any_transcript": precomputed.get("max_any_transcript", ""),
            "max_delta_field": precomputed.get("max_delta_field", ""),
            "delta_scores": precomputed.get("delta_scores", {}),
            "reference_scores": precomputed.get("reference_scores", {}),
            "alternate_scores": precomputed.get("alternate_scores", {}),
            "cache_key": precomputed.get("cache_key"),
            "source": precomputed.get("source"),
            "grch38": precomputed.get("grch38"),
            "variant": precomputed.get("variant"),
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
            "distance": SPLICEAI_MAX_DISTANCE,
            "mask": SPLICEAI_MASK,
            "annotation_subset": SPLICEAI_ANNOTATION_SUBSET,
            "genome_assembly": SPLICEAI_GENOME_ASSEMBLY,
            "aggregation": SPLICEAI_AGGREGATION,
        }
        return score

    # 3. persistent Broad API cache
    api_cache = _load_api_cache()
    if cache_key in api_cache:
        entry = api_cache[cache_key]
        if _runtime_entry_matches_profile(entry):
            score = float(entry["score"])
            SPLICEAI_CACHE[cache_key] = score
            SPLICEAI_STATUS_CACHE[variant_key] = {
                "status": "ok",
                "score": score,
                "reason": "Loaded from persistent Broad API runtime cache",
                "source": entry.get("source") or entry.get("api_source") or "Broad SpliceAI API runtime cache",
                "transcript_policy": entry.get("transcript_policy"),
                "selected_transcript": entry.get("selected_transcript"),
                "reference_transcript_score": entry.get("reference_transcript_score"),
                "max_any_transcript_score": entry.get("max_any_transcript_score"),
                "max_any_transcript": entry.get("max_any_transcript", ""),
                "max_delta_field": entry.get("max_delta_field", ""),
                "delta_scores": entry.get("delta_scores", {}),
                "reference_scores": entry.get("reference_scores", {}),
                "alternate_scores": entry.get("alternate_scores", {}),
                "grch38": f"{entry.get('chrom')}:{entry.get('pos')}:{entry.get('ref')}>{entry.get('alt')}",
                "cache_key": cache_key,
                "scoring_profile_id": SPLICEAI_PROFILE_ID,
                "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
                "distance": SPLICEAI_MAX_DISTANCE,
                "mask": SPLICEAI_MASK,
                "annotation_subset": SPLICEAI_ANNOTATION_SUBSET,
                "genome_assembly": SPLICEAI_GENOME_ASSEMBLY,
                "aggregation": SPLICEAI_AGGREGATION,
            }
            return score
        register_issue(
            "SpliceAI API cache",
            "ignored a runtime record created with an incompatible or incomplete scoring profile",
        )

    # 4. need GRCh38 coords to call API
    resolved = {}
    resolved_variant = resolve_variant(gene, c_notation)
    if resolved_variant and resolved_variant.status != "failed":
        resolved[variant_key] = resolved_variant

    coords = get_grch38(resolved, gene, c_notation)
    if coords is None:
        SPLICEAI_STATUS_CACHE[variant_key] = {
            "status": "no_grch38_coords",
            "score":  None,
            "reason": "No GRCh38 coordinates available",
            "source": "Broad SpliceAI API",
            "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "distance": SPLICEAI_MAX_DISTANCE,
            "mask": SPLICEAI_MASK,
        }
        return None

    # 5. live API call
    time.sleep(SPLICEAI_API_RATE_SLEEP)
    selected = _query_spliceai_api(
        gene, coords["chrom"], coords["pos"], coords["ref"], coords["alt"]
    )

    if selected is None or selected.get("score") is None:
        SPLICEAI_STATUS_CACHE[variant_key] = {
            "status": "api_error",
            "score":  None,
            "reason": (
                selected.get("error") if isinstance(selected, dict) and selected.get("error")
                else "Broad SpliceAI API returned no score for the required transcript"
            ),
            "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
            "source": SPLICEAI_API_SOURCE,
            "grch38": f"{coords['chrom']}:{coords['pos']}:{coords['ref']}>{coords['alt']}",
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "distance": SPLICEAI_MAX_DISTANCE,
            "mask": SPLICEAI_MASK,
        }
        return None

    score = selected["score"]

    # cache result
    SPLICEAI_CACHE[cache_key] = score
    api_cache[cache_key] = {
        "score":   score,
        "chrom":   str(coords["chrom"]),
        "pos":     coords["pos"],
        "ref":     coords["ref"],
        "alt":     coords["alt"],
        "source":  SPLICEAI_API_SOURCE,
        "api_source": selected.get("api_source"),
        "api_url": SPLICEAI_API_URL,
        "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
        "selected_transcript": selected.get("selected_transcript"),
        "max_delta_field": selected.get("max_delta_field"),
        "delta_scores": selected.get("delta_scores", {}),
        "reference_scores": selected.get("reference_scores", {}),
        "alternate_scores": selected.get("alternate_scores", {}),
        "reference_transcript_score": selected.get("reference_transcript_score"),
        "reference_transcript": selected.get("reference_transcript"),
        "max_any_transcript_score": selected.get("max_any_transcript_score"),
        "max_any_transcript": selected.get("max_any_transcript"),
        "n_transcript_scores": selected.get("n_transcript_scores"),
        "scoring_profile_id": SPLICEAI_PROFILE_ID,
        "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
        "genome_assembly": SPLICEAI_GENOME_ASSEMBLY,
        "distance": SPLICEAI_MAX_DISTANCE,
        "mask": SPLICEAI_MASK,
        "annotation_subset": SPLICEAI_ANNOTATION_SUBSET,
        "aggregation": SPLICEAI_AGGREGATION,
    }
    cache_saved = _save_api_cache(api_cache)

    SPLICEAI_STATUS_CACHE[variant_key] = {
        "status": "ok",
        "score":  score,
        "reason": (
            f"Queried from {SPLICEAI_API_SOURCE} and persisted to the runtime cache"
            if cache_saved
            else f"Queried from {SPLICEAI_API_SOURCE}; available in memory but not persisted"
        ),
        "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
        "source": SPLICEAI_API_SOURCE,
        "selected_transcript": selected.get("selected_transcript"),
        "reference_transcript_score": selected.get("reference_transcript_score"),
        "max_any_transcript_score": selected.get("max_any_transcript_score"),
        "max_any_transcript": selected.get("max_any_transcript", ""),
        "max_delta_field": selected.get("max_delta_field", ""),
        "delta_scores": selected.get("delta_scores", {}),
        "reference_scores": selected.get("reference_scores", {}),
        "alternate_scores": selected.get("alternate_scores", {}),
        "grch38": f"{coords['chrom']}:{coords['pos']}:{coords['ref']}>{coords['alt']}",
        "cache_key": cache_key,
        "scoring_profile_id": SPLICEAI_PROFILE_ID,
        "scoring_profile_sha256": SPLICEAI_PROFILE_SHA256,
        "distance": SPLICEAI_MAX_DISTANCE,
        "mask": SPLICEAI_MASK,
        "annotation_subset": SPLICEAI_ANNOTATION_SUBSET,
        "genome_assembly": SPLICEAI_GENOME_ASSEMBLY,
        "aggregation": SPLICEAI_AGGREGATION,
    }
    return score


# ============================================================
# SpliceAI criterion helper functions
# ============================================================

SPLICEAI_PP3_ALLOWED_TYPES = {
    "synonymous", "silent",
    "missense",
    "inframe_deletion", "inframe_insertion", "inframe_delins", "delins",
    "intronic",
}

SPLICEAI_BP4_ALLOWED_TYPES = {
    "synonymous", "silent",
    "missense",
    "inframe_deletion", "inframe_insertion", "inframe_delins", "delins",
    "intronic",
}


def normalize_variant_type(variant_type: str) -> str:
    return (variant_type or "").strip().lower()


def spliceai_is_confirmed_low(score: Optional[float]) -> bool:
    """True only when SpliceAI is available and <= 0.10."""
    return score is not None and score <= SPLICEAI_LOW_THRESHOLD


def spliceai_predicts_splice_effect(score: Optional[float]) -> bool:
    """True only when SpliceAI is available and >= 0.20."""
    return score is not None and score >= SPLICEAI_HIGH_THRESHOLD


def variant_type_allows_spliceai_pp3(variant_type: str) -> bool:
    """
    Guardrail for PP3 from SpliceAI.

    PP3-SpliceAI is not a generic "any variant" rule. It should not be added
    to nonsense/PTC, frameshift, exon-deletion, or canonical splice-site variants
    where the same loss-of-function/splicing mechanism is evaluated through
    PVS1/RNA logic.
    """
    return normalize_variant_type(variant_type) in SPLICEAI_PP3_ALLOWED_TYPES


def variant_type_allows_spliceai_bp4(variant_type: str) -> bool:
    return normalize_variant_type(variant_type) in SPLICEAI_BP4_ALLOWED_TYPES

def spliceai_lookup_report(gene: str, c_notation: str) -> dict:
    """
    Return a small diagnostic object for one variant.
    This is useful when checking why a score was not used.
    """
    key = f"{gene}:{c_notation}"
    score = get_spliceai_score(gene, c_notation)
    status = SPLICEAI_STATUS_CACHE.get(key, {})
    resolved = {}
    resolved_variant = resolve_variant(gene, c_notation)
    if resolved_variant and resolved_variant.status != "failed":
        resolved[key] = resolved_variant
    coords = get_grch38(resolved, gene, c_notation)
    return {
        "variant": key,
        "coords": coords,
        "score": score,
        "status": status.get("status"),
        "reason": status.get("reason"),
    }


if __name__ == "__main__":
    print(f"SpliceAI API cache: {SPLICEAI_API_CACHE_PATH}")
    existing = _load_api_cache()
    print(f"Cached variants:    {len(existing)}")
