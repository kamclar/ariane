# ============================================================
# SpliceAI lookup via precomputed BRCA SNV cache + Broad API fallback
#
# Primary source for BRCA1/2 coding SNVs is the locally precomputed
# reference-transcript cache in data/spliceai. The Broad Institute SpliceAI
# Lookup API remains the fallback for variants outside that cache and for
# periodic validation checks.
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

from backend.lookups.coordinates import resolve_variant, get_grch38


def choose_project_root() -> Path:
    env_root = os.environ.get("BRCA_ACMG_PROJECT_ROOT")
    if env_root:
        return Path(env_root)
    candidates = [
        Path("/content/drive/MyDrive/Enigma"),
        Path("/content/drive/MyDrive/BRCA_ACMG"),
        Path("."),
    ]
    for root in candidates:
        if (root / "data").exists():
            return root
    if Path("/content/drive/MyDrive").exists():
        return Path("/content/drive/MyDrive/Enigma")
    return Path(".")

PROJECT_ROOT  = choose_project_root()
SPLICEAI_DIR  = PROJECT_ROOT / "data" / "spliceai"
SPLICEAI_DIR.mkdir(parents=True, exist_ok=True)

# Local JSON cache - persists across Colab sessions via Drive
SPLICEAI_API_CACHE_PATH = SPLICEAI_DIR / "spliceai_api_cache.json"
SPLICEAI_PRECOMPUTED_CACHE_PATH = Path(os.environ.get(
    "SPLICEAI_PRECOMPUTED_CACHE_PATH",
    SPLICEAI_DIR / "spliceai_brca_snv_reference_cache.json",
))

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
SPLICEAI_API_TIMEOUT = _env_int("SPLICEAI_API_TIMEOUT", 10)
SPLICEAI_API_RATE_SLEEP = _env_float("SPLICEAI_API_RATE_SLEEP", 1.5)
SPLICEAI_API_SOURCE = os.environ.get(
    "SPLICEAI_API_SOURCE",
    "Local Broad SpliceAI API" if "localhost" in SPLICEAI_API_URL else "Broad SpliceAI API",
)

REFERENCE_TRANSCRIPTS = {
    "BRCA1": {
        "refseq": "NM_007294.4",
        "ensembl": "ENST00000357654.9",
    },
    "BRCA2": {
        "refseq": "NM_000059.4",
        "ensembl": "ENST00000380152.8",
    },
}

SPLICEAI_TRANSCRIPT_POLICY = os.environ.get(
    "SPLICEAI_TRANSCRIPT_POLICY",
    "reference_transcript",
).strip().lower()

SPLICEAI_USE_PRECOMPUTED_CACHE = os.environ.get(
    "SPLICEAI_USE_PRECOMPUTED_CACHE",
    "1",
).strip().lower() not in {"0", "false", "no", "off"}

if SPLICEAI_TRANSCRIPT_POLICY not in {"reference_transcript", "max_any_transcript"}:
    SPLICEAI_TRANSCRIPT_POLICY = "reference_transcript"


def _load_api_cache() -> dict:
    """Load persisted API cache from Drive."""
    if SPLICEAI_API_CACHE_PATH.exists():
        try:
            with open(SPLICEAI_API_CACHE_PATH) as f:
                return _json.load(f)
        except Exception:
            pass
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
    if not SPLICEAI_PRECOMPUTED_CACHE_PATH.exists():
        return SPLICEAI_PRECOMPUTED_CACHE

    try:
        with open(SPLICEAI_PRECOMPUTED_CACHE_PATH, encoding="utf-8") as handle:
            raw = _json.load(handle)
    except Exception as exc:
        print(f"Warning: could not load precomputed SpliceAI cache: {exc}")
        return SPLICEAI_PRECOMPUTED_CACHE

    if isinstance(raw, dict):
        SPLICEAI_PRECOMPUTED_CACHE = raw
    return SPLICEAI_PRECOMPUTED_CACHE


def _save_api_cache(cache: dict) -> None:
    """Persist API cache to Drive."""
    try:
        with open(SPLICEAI_API_CACHE_PATH, "w") as f:
            _json.dump(cache, f, indent=2)
    except Exception as e:
        print(f"Warning: could not save SpliceAI cache: {e}")


def _cache_key(gene: str, c_notation: str) -> str:
    return f"{SPLICEAI_TRANSCRIPT_POLICY}:{gene}:{c_notation}"


def _precomputed_cache_keys(gene: str, c_notation: str) -> tuple[str, str]:
    raw_key = f"{gene}:{c_notation}"
    return raw_key, f"reference_transcript:{raw_key}"


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
        if entry.get("status") not in (None, "ok"):
            continue
        score = _entry_score(entry)
        if score is None:
            continue
        return {
            "score": score,
            "max_delta_field": entry.get("max_delta_field", ""),
            "selected_transcript": REFERENCE_TRANSCRIPTS.get(gene, {}).get("ensembl", ""),
            "reference_transcript_score": score,
            "max_any_transcript_score": entry.get("max_any_transcript_score"),
            "source": entry.get("source") or "precomputed BRCA SNV SpliceAI cache",
            "cache_key": key,
            "grch38": entry.get("grch38", ""),
            "variant": entry.get("variant", ""),
        }
    return None


def _score_row(row: dict) -> tuple[float, str]:
    best_score = 0.0
    best_field = ""
    for key in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
        try:
            value = float(row.get(key, 0) or 0)
        except (ValueError, TypeError):
            value = 0.0
        if value > best_score:
            best_score = value
            best_field = key
    return best_score, best_field


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


def _select_spliceai_score(gene: str, scores: list[dict]) -> dict:
    max_any_score = None
    max_any_field = ""
    max_any_transcript = ""
    ref_score = None
    ref_field = ""
    ref_transcript = ""

    for row in scores:
        row_score, row_field = _score_row(row)
        if max_any_score is None or row_score > max_any_score:
            max_any_score = row_score
            max_any_field = row_field
            max_any_transcript = str(row.get("t_id") or "")
        if _row_matches_reference_transcript(gene, row):
            ref_score = row_score
            ref_field = row_field
            ref_transcript = str(row.get("t_id") or "")

    if SPLICEAI_TRANSCRIPT_POLICY == "max_any_transcript":
        return {
            "score": max_any_score,
            "max_delta_field": max_any_field,
            "selected_transcript": max_any_transcript,
            "selected_transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
            "reference_transcript_score": ref_score,
            "reference_transcript": ref_transcript,
            "max_any_transcript_score": max_any_score,
            "max_any_transcript": max_any_transcript,
        }

    return {
        "score": ref_score,
        "max_delta_field": ref_field,
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
    reference = REFERENCE_TRANSCRIPTS.get(gene)
    transcript_arg = ""
    if SPLICEAI_TRANSCRIPT_POLICY == "reference_transcript" and reference:
        transcript_arg = f"&transcript={urllib.parse.quote(reference['ensembl'])}"
    url = (
        f"{SPLICEAI_API_URL}"
        f"?variant={urllib.parse.quote(variant_str)}"
        f"&hg=38"
        f"&distance=50"
        f"&mask=0"
        f"{transcript_arg}"
    )
    try:
        req = urllib.request.Request(
            url,
            headers={"Accept": "application/json", "User-Agent": "BRCA-ACMG-Module1/1.6.4"},
        )
        with urllib.request.urlopen(req, timeout=SPLICEAI_API_TIMEOUT) as resp:
            data = _json.loads(resp.read())

        scores = data.get("scores", [])
        if not scores:
            return None  # no scores from API - do not treat as confirmed low

        selected = _select_spliceai_score(gene, scores)
        selected["api_source"] = data.get("source") or SPLICEAI_API_SOURCE
        selected["n_transcript_scores"] = len(scores)
        return selected

    except Exception as e:
        return None


def get_spliceai_score(gene: str, c_notation: str) -> Optional[float]:
    """
    Look up SpliceAI score for a variant via Broad API with local Drive cache.

    Returns a float score or None. None means unavailable, not 0.0.
    Benign criteria must only use confirmed scores <= 0.1.

    Lookup order:
      1. In-memory cache (fast, within session)
      2. Precomputed BRCA coding SNV cache (reference transcript only)
      3. Persistent Broad API cache
      4. Broad SpliceAI API fallback
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
            "cache_key": precomputed.get("cache_key"),
            "source": precomputed.get("source"),
            "grch38": precomputed.get("grch38"),
            "variant": precomputed.get("variant"),
        }
        return score

    # 3. persistent Broad API cache
    api_cache = _load_api_cache()
    if cache_key in api_cache:
        entry = api_cache[cache_key]
        score = entry.get("score")
        SPLICEAI_CACHE[cache_key] = score
        SPLICEAI_STATUS_CACHE[variant_key] = {
            "status": "ok" if score is not None else "api_no_score",
            "score":  score,
            "reason": "Loaded from Drive cache (Broad API)",
            "transcript_policy": entry.get("transcript_policy"),
            "selected_transcript": entry.get("selected_transcript"),
            "reference_transcript_score": entry.get("reference_transcript_score"),
            "max_any_transcript_score": entry.get("max_any_transcript_score"),
        }
        return score

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
            "reason": "Broad SpliceAI API call failed, returned no scores, or did not include the required transcript",
            "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
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
        "reference_transcript_score": selected.get("reference_transcript_score"),
        "reference_transcript": selected.get("reference_transcript"),
        "max_any_transcript_score": selected.get("max_any_transcript_score"),
        "max_any_transcript": selected.get("max_any_transcript"),
        "n_transcript_scores": selected.get("n_transcript_scores"),
    }
    _save_api_cache(api_cache)

    SPLICEAI_STATUS_CACHE[variant_key] = {
        "status": "ok",
        "score":  score,
        "reason": f"Queried from {SPLICEAI_API_SOURCE} and cached to Drive",
        "transcript_policy": SPLICEAI_TRANSCRIPT_POLICY,
        "selected_transcript": selected.get("selected_transcript"),
        "reference_transcript_score": selected.get("reference_transcript_score"),
        "max_any_transcript_score": selected.get("max_any_transcript_score"),
    }
    return score


# ============================================================
# SpliceAI criterion helper functions
# ============================================================

SPLICEAI_LOW_THRESHOLD = 0.10
SPLICEAI_HIGH_THRESHOLD = 0.20

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
