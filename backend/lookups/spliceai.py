# ============================================================
# SpliceAI lookup via an ENIGMA-compatible API and profile-pinned runtime cache
#
# Scores are computed on demand. ARIANE does not use a precomputed gene-wide
# variant space as a classification source. The persistent cache contains only
# results that were returned by the configured service for variants that were
# actually requested.
#
# Default endpoint: http://127.0.0.1:8082/spliceai/
# Variant format: chr{chrom}-{pos}-{ref}-{alt}
# Successful results are kept in the profile-specific runtime cache.
#
# The previous MANE VCF subset approach was removed because the Ensembl MANE v1.0
# file uses an older Gencode version and gives incorrect scores for some variants
# (e.g. BRCA1 c.4185G>A: MANE gives DS_DL=0.01, Broad API gives DS_DL=0.93).
# ============================================================
from typing import Optional, Dict
from pathlib import Path
import json
import socket
import threading
import time
import urllib.error
import urllib.request
import urllib.parse
import json as _json
import os
import tempfile
from backend.infrastructure.health import DataHealthRegistry
from backend.infrastructure.runtime_cache import choose_runtime_cache_dir
from backend.policy.spliceai_profile import (
    SPLICEAI_AGGREGATION,
    SPLICEAI_ALTERNATE_FIELDS,
    SPLICEAI_APPROVED_DOCKER_IMAGE,
    SPLICEAI_ANNOTATION_SUBSET,
    SPLICEAI_DELTA_FIELDS,
    SPLICEAI_GENOME_ASSEMBLY,
    SPLICEAI_MASK,
    SPLICEAI_MAX_DISTANCE,
    SPLICEAI_PROFILE_ID,
    SPLICEAI_PROFILE_SHA256,
    SPLICEAI_REFERENCE_FIELDS,
    SPLICEAI_REFERENCE_TRANSCRIPTS,
    SPLICEAI_TRANSCRIPT_POLICY_REQUIRED,
)

from backend.lookups.coordinates import resolve_variant, get_grch38
from backend.version import ARIANE_VERSION


def choose_project_root() -> Path:
    env_root = os.environ.get("BRCA_ACMG_PROJECT_ROOT")
    if env_root:
        return Path(env_root)
    # The repository root is stable regardless of the process working directory.
    # Deployments with external data must opt in through BRCA_ACMG_PROJECT_ROOT.
    return Path(__file__).resolve().parents[2]

PROJECT_ROOT  = choose_project_root()
SPLICEAI_DIR  = PROJECT_ROOT / "data" / "spliceai"

RUNTIME_CACHE_DIR = choose_runtime_cache_dir()

# Mutable API results are separate from immutable, versioned snapshots.
SPLICEAI_API_CACHE_PATH = RUNTIME_CACHE_DIR / "spliceai_api_cache.json"
# In-memory caches
SPLICEAI_CACHE:        Dict[str, float] = {}   # policy:gene:c_notation -> score
SPLICEAI_STATUS_CACHE: Dict[str, dict]  = {}   # gene:c_notation -> status details

# ARIANE uses its private, digest-pinned service by default. An explicitly
# configured compatible endpoint is accepted, but it is never used as a
# fallback after the configured source fails.
DEFAULT_SPLICEAI_API_URL = "http://127.0.0.1:8082/spliceai/"


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
# Interactive requests are proxied by nginx with a 180-second read timeout.
# Required Figure 1A evidence must complete within the bounded lookup deadline
# or the classifier returns no classification. Offline cache builders do not
# use this web-request deadline.
SPLICEAI_API_TIMEOUT = _env_int("SPLICEAI_API_TIMEOUT", 120)
SPLICEAI_API_RATE_SLEEP = _env_float("SPLICEAI_API_RATE_SLEEP", 0.0)
SPLICEAI_API_ATTEMPTS = max(1, _env_int("SPLICEAI_API_ATTEMPTS", 1))
SPLICEAI_API_RETRY_DELAY = _env_float("SPLICEAI_API_RETRY_DELAY", 2.0)
SPLICEAI_API_MAX_CONCURRENT = max(
    1, _env_int("SPLICEAI_API_MAX_CONCURRENT", 2)
)
SPLICEAI_API_SOURCE = os.environ.get(
    "SPLICEAI_API_SOURCE",
    (
        "ARIANE local SpliceAI service"
        if urllib.parse.urlsplit(SPLICEAI_API_URL).hostname in {"127.0.0.1", "localhost", "::1"}
        else "Configured SpliceAI API"
    ),
)

REFERENCE_TRANSCRIPTS = SPLICEAI_REFERENCE_TRANSCRIPTS

_API_REQUEST_GATE = threading.BoundedSemaphore(SPLICEAI_API_MAX_CONCURRENT)
_API_RATE_LOCK = threading.Lock()
_API_NEXT_REQUEST_AT = 0.0


def spliceai_runtime_health() -> dict:
    """Report local service reachability without running model inference."""
    parsed = urllib.parse.urlsplit(SPLICEAI_API_URL)
    host = parsed.hostname or ""
    local = host in {"127.0.0.1", "localhost", "::1"}
    result = {
        "status": "configured",
        "source": SPLICEAI_API_SOURCE,
        "local": local,
        "profile_id": SPLICEAI_PROFILE_ID,
        "docker_image": SPLICEAI_APPROVED_DOCKER_IMAGE,
    }
    if not local:
        return result
    port = parsed.port or (443 if parsed.scheme == "https" else 80)
    try:
        with socket.create_connection((host, port), timeout=0.25):
            pass
    except OSError as exc:
        result.update({"status": "unavailable", "reason": f"{type(exc).__name__}: {exc}"})
        return result
    result["status"] = "ok"
    return result

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

def _load_api_cache(health: DataHealthRegistry | None = None) -> dict:
    """Load the persistent runtime API cache."""
    RUNTIME_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if SPLICEAI_API_CACHE_PATH.exists():
        try:
            with open(SPLICEAI_API_CACHE_PATH) as f:
                result = _json.load(f)
            if not isinstance(result, dict):
                raise ValueError("runtime cache root must be a JSON object")
            result = _current_profile_cache_entries(result)
            if health is not None:
                health.clear("SpliceAI API cache")
            return result
        except Exception as exc:
            if health is not None:
                health.register(
                    "SpliceAI API cache",
                    f"could not load {SPLICEAI_API_CACHE_PATH}: {type(exc).__name__}: {exc}",
                )
    return {}


def _current_profile_cache_entries(cache: dict) -> dict:
    """Discard entries from retired profiles before the next cache write."""
    prefix = f"{SPLICEAI_PROFILE_ID}:{SPLICEAI_TRANSCRIPT_POLICY}:"
    return {
        key: value
        for key, value in cache.items()
        if isinstance(key, str) and key.startswith(prefix)
    }


def _save_api_cache(
    cache: dict,
    health: DataHealthRegistry | None = None,
) -> bool:
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
        if health is not None:
            health.clear("SpliceAI API cache")
        return True
    except Exception as e:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError:
                pass
        print(f"Warning: could not save SpliceAI cache: {e}")
        if health is not None:
            health.register(
                "SpliceAI API cache",
                f"score was obtained and used, but the runtime cache could not be saved to "
                f"{SPLICEAI_API_CACHE_PATH}; this request is unaffected, but the score may need "
                f"to be fetched again after restart: {type(e).__name__}: {e}",
            )
        return False


def _cache_key(gene: str, c_notation: str) -> str:
    return f"{SPLICEAI_PROFILE_ID}:{SPLICEAI_TRANSCRIPT_POLICY}:{gene}:{c_notation}"


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


def _score_entry_is_complete(entry: object) -> bool:
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


def _runtime_entry_matches_profile(entry: object, gene: str) -> bool:
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
            entry.get("selected_transcript")
            != REFERENCE_TRANSCRIPTS.get(gene, {}).get("ensembl"),
        )
    ):
        return False
    return _score_entry_is_complete({**entry, "status": "ok"})


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
    return (
        str(row.get("t_id") or "") == reference["ensembl"]
        and reference["refseq"]
        in {str(value) for value in row.get("t_refseq_ids") or []}
    )


def _reference_transcript_match_rank(gene: str, row: dict) -> int:
    # The active profile pins transcript versions. A versionless match would
    # make a later annotation release look equivalent to the reviewed source.
    return 1 if _row_matches_reference_transcript(gene, row) else 0


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

    if SPLICEAI_TRANSCRIPT_POLICY == "reference_transcript" and reference_row is None:
        reference = REFERENCE_TRANSCRIPTS.get(gene) or {}
        return {
            "score": None,
            "error": (
                "SpliceAI response has no complete score row for the exact "
                f"reference transcript {reference.get('ensembl') or 'unknown'} / "
                f"{reference.get('refseq') or 'unknown'}"
            ),
        }

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


def _query_spliceai_api_once(
    gene: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[dict]:
    """
    Query the configured SpliceAI API for a single variant.
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
            headers={"Accept": "application/json", "User-Agent": f"ARIANE/{ARIANE_VERSION}"},
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

    except urllib.error.HTTPError as exc:
        return {
            "score": None,
            "error": f"HTTP Error {exc.code}: {exc.reason}",
            "http_status": exc.code,
            "retryable": exc.code in {429, 500, 502, 503, 504},
        }
    except (TimeoutError, socket.timeout) as exc:
        return {
            "score": None,
            "error": f"{type(exc).__name__}: {exc}",
            "retryable": True,
        }
    except urllib.error.URLError as exc:
        return {
            "score": None,
            "error": f"URLError: {exc.reason}",
            "retryable": True,
        }
    except Exception as exc:
        return {
            "score": None,
            "error": f"{type(exc).__name__}: {exc}",
            "retryable": False,
        }


def _wait_for_api_rate_slot() -> None:
    global _API_NEXT_REQUEST_AT
    with _API_RATE_LOCK:
        now = time.monotonic()
        delay = max(0.0, _API_NEXT_REQUEST_AT - now)
        if delay:
            time.sleep(delay)
        _API_NEXT_REQUEST_AT = time.monotonic() + SPLICEAI_API_RATE_SLEEP


def _query_spliceai_api(
    gene: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[dict]:
    """Query the configured source with bounded transient retries."""
    last_result: Optional[dict] = None
    for attempt in range(SPLICEAI_API_ATTEMPTS):
        with _API_REQUEST_GATE:
            _wait_for_api_rate_slot()
            last_result = _query_spliceai_api_once(gene, chrom, pos, ref, alt)
        if last_result is None or last_result.get("score") is not None:
            return last_result
        if not last_result.get("retryable") or attempt + 1 == SPLICEAI_API_ATTEMPTS:
            return last_result
        time.sleep(SPLICEAI_API_RETRY_DELAY * (2 ** attempt))
    return last_result


def get_spliceai_score(
    gene: str,
    c_notation: str,
    *,
    health: DataHealthRegistry | None = None,
) -> Optional[float]:
    """
    Look up SpliceAI score through the profile-pinned API runtime path.

    Returns a float score or None. None means unavailable, not 0.0.
    Benign criteria must only use confirmed scores <= 0.1.

    Lookup order:
      1. In-memory cache (fast, within session)
      2. Persistent profile-pinned runtime cache
      3. Configured ENIGMA-compatible SpliceAI API
    """
    variant_key = f"{gene}:{c_notation}"
    cache_key = _cache_key(gene, c_notation)

    # 1. in-memory cache
    if cache_key in SPLICEAI_CACHE:
        return SPLICEAI_CACHE[cache_key]

    # 2. persistent runtime cache
    api_cache = _load_api_cache(health)
    if cache_key in api_cache:
        entry = api_cache[cache_key]
        if _runtime_entry_matches_profile(entry, gene):
            score = float(entry["score"])
            SPLICEAI_CACHE[cache_key] = score
            SPLICEAI_STATUS_CACHE[variant_key] = {
                "status": "ok",
                "score": score,
                "reason": "Loaded from persistent SpliceAI runtime cache",
                "source": entry.get("source") or entry.get("api_source") or "SpliceAI runtime cache",
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
        if health is not None:
            health.register(
                "SpliceAI API cache",
                "ignored a runtime record created with an incompatible or incomplete scoring profile",
            )

    # 3. need GRCh38 coords to call API
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
            "source": SPLICEAI_API_SOURCE,
            "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "distance": SPLICEAI_MAX_DISTANCE,
            "mask": SPLICEAI_MASK,
        }
        return None

    # 4. live API call
    selected = _query_spliceai_api(
        gene, coords["chrom"], coords["pos"], coords["ref"], coords["alt"]
    )

    if selected is None or selected.get("score") is None:
        SPLICEAI_STATUS_CACHE[variant_key] = {
            "status": "api_error",
            "score":  None,
            "reason": (
                selected.get("error") if isinstance(selected, dict) and selected.get("error")
                else "Configured SpliceAI API returned no score for the required transcript"
            ),
            "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
            "source": SPLICEAI_API_SOURCE,
            "grch38": f"{coords['chrom']}:{coords['pos']}:{coords['ref']}>{coords['alt']}",
            "scoring_profile_id": SPLICEAI_PROFILE_ID,
            "distance": SPLICEAI_MAX_DISTANCE,
            "mask": SPLICEAI_MASK,
            "retryable": bool(
                selected.get("retryable")
                if isinstance(selected, dict)
                else False
            ),
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
    cache_saved = _save_api_cache(api_cache, health)

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


def get_spliceai_status(gene: str, c_notation: str) -> dict:
    """Return recorded provenance/status for the most recent score lookup."""
    return dict(SPLICEAI_STATUS_CACHE.get(f"{gene}:{c_notation}", {}))

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
