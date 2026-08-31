# BayesDel_noAF + AlphaMissense lookup via myvariant.info
#
# Single API call returns both scores:
#   dbnsfp.bayesdel.no_af.score  - BayesDel_noAF
#   alphamissense.am_pathogenicity / am_class  - AlphaMissense (missense only)
#
# AlphaMissense is informational only - not used for ENIGMA VCEP scoring.
#
# Persistent cache: ${ARIANE_RUNTIME_CACHE_DIR}/bayesdel_api_cache.json
# Cache entry also preserves the lookup status/reason so absence is explainable.
# Cache entries created before transcript-safe selection are ignored.

from typing import Optional, Dict, Tuple
from pathlib import Path
import json
import os
import tempfile
import threading
import urllib.error
import urllib.parse
import urllib.request

from backend.data_health import clear_issue, register_issue
from backend.gene_policy import reference_transcript
from backend.lookups import coordinates
from backend.runtime_cache import runtime_cache_path
from backend.version import ARIANE_VERSION

MYVARIANT_BASE_URL = "https://myvariant.info/v1/variant"
BAYESDEL_SELECTION_POLICY = "unambiguous_variant_score_v1"

# Cache stores dicts: {"bayesdel": ..., "am_score": ..., "am_class": ...}
BAYESDEL_CACHE: Dict[str, Optional[dict]] = {}
BAYESDEL_STATUS_CACHE: Dict[str, dict] = {}

_CACHE_PATH = runtime_cache_path("bayesdel_api_cache.json")
_FILE_LOCK  = threading.Lock()


def _load_cache() -> None:
    global BAYESDEL_CACHE
    if not _CACHE_PATH.exists():
        clear_issue("BayesDel cache")
        return
    try:
        with open(_CACHE_PATH, encoding="utf-8") as fh:
            raw = json.load(fh)
        for key, val in raw.items():
            if isinstance(val, dict):
                # A successful response without either annotation can be
                # transient (for example while the upstream annotation
                # backend is degraded).  It must not become a permanent
                # negative result that prevents a later retry.
                if (
                    val.get("status") == "no_score"
                    or (
                        val.get("status") in {None, ""}
                        and val.get("bayesdel") is None
                        and val.get("am_score") is None
                    )
                ):
                    continue
                if (
                    val.get("bayesdel") is not None
                    and val.get("selection_policy") != BAYESDEL_SELECTION_POLICY
                ):
                    # Older cache entries may contain a maximum selected from
                    # divergent upstream values. Their provenance is not safe
                    # enough for PP3/BP4 and they must be fetched again.
                    continue
                BAYESDEL_CACHE[key] = val
            # Float-only legacy entries have no selection provenance and are
            # intentionally ignored.
        print(f"Loaded BayesDel/AM cache: {len(BAYESDEL_CACHE)} entries")
        clear_issue("BayesDel cache")
    except Exception as exc:
        print(f"Warning: could not load BayesDel cache: {exc}")
        register_issue("BayesDel cache", f"could not load {_CACHE_PATH}: {type(exc).__name__}: {exc}")


def _save_cache() -> None:
    with _FILE_LOCK:
        temporary_path = None
        try:
            _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=_CACHE_PATH.parent,
                prefix="bayesdel_api_cache.",
                suffix=".tmp",
                delete=False,
            ) as fh:
                temporary_path = Path(fh.name)
                json.dump(BAYESDEL_CACHE, fh)
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(temporary_path, _CACHE_PATH)
            clear_issue("BayesDel cache")
        except Exception as exc:
            if temporary_path is not None:
                try:
                    temporary_path.unlink(missing_ok=True)
                except OSError:
                    pass
            print(f"Warning: could not save BayesDel cache: {exc}")
            register_issue(
                "BayesDel cache",
                "score was obtained and used, but the runtime cache could not "
                f"be saved to {_CACHE_PATH}; this request is unaffected, but the "
                "score may need to be fetched again after restart: "
                f"{type(exc).__name__}: {exc}",
            )


def _numeric_values(value) -> list[float]:
    values = value if isinstance(value, list) else [value]
    parsed: list[float] = []
    for item in values:
        if item is None or isinstance(item, bool):
            continue
        try:
            parsed.append(float(item))
        except (TypeError, ValueError):
            continue
    return parsed


def _select_unambiguous_bayesdel(value) -> tuple[Optional[float], str, Optional[str]]:
    """Accept only a score whose value cannot depend on transcript selection."""
    values = _numeric_values(value)
    if not values:
        return None, "not_available", None
    unique = sorted(set(values))
    if len(unique) == 1:
        basis = "single_value" if len(values) == 1 else "all_returned_values_identical"
        return unique[0], basis, None
    return (
        None,
        "ambiguous_multiple_values",
        "MyVariant returned multiple distinct BayesDel_noAF values without an "
        "explicit value-to-transcript mapping",
    )


def _returned_transcripts(dbnsfp: dict) -> list[str]:
    ensembl = dbnsfp.get("ensembl", {}) if isinstance(dbnsfp, dict) else {}
    value = ensembl.get("transcriptid") if isinstance(ensembl, dict) else None
    values = value if isinstance(value, list) else [value]
    return [str(item) for item in values if item]


def fetch_variant_data_myvariant(gene: str, c_notation: str, hg37_coords: Optional[dict]) -> dict:
    """Fetch BayesDel + AlphaMissense from myvariant.info in one request."""
    result = {
        "bayesdel": None, "am_score": None, "am_class": None,
        "status": "not_queried", "error": None,
        "selection_policy": BAYESDEL_SELECTION_POLICY,
        "selection_basis": "not_available",
        "reference_transcript": reference_transcript(gene),
        "returned_transcripts": [],
    }

    if hg37_coords is None:
        result["status"] = "no_grch37_coords"
        return result

    chrom = hg37_coords["chrom"]
    pos   = hg37_coords["pos"]
    ref   = hg37_coords["ref"]
    alt   = hg37_coords["alt"]

    hgvs = f"chr{chrom}:g.{pos}{ref}>{alt}"
    url = (f"{MYVARIANT_BASE_URL}/{urllib.parse.quote(hgvs)}"
           f"?{urllib.parse.urlencode({'fields': 'dbnsfp.bayesdel,dbnsfp.ensembl.transcriptid,alphamissense'})}")
    req = urllib.request.Request(url, headers={"User-Agent": f"ARIANE/{ARIANE_VERSION}"})

    try:
        with urllib.request.urlopen(req, timeout=20) as response:
            data = json.loads(response.read())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            result["status"] = "not_found"
            return result
        result["status"] = "api_error"
        result["error"] = f"HTTPError {e.code}: {e.reason}"
        return result
    except Exception as e:
        result["status"] = "api_error"
        result["error"] = f"{type(e).__name__}: {e}"
        return result

    if data.get("notfound"):
        result["status"] = "not_found"
        return result

    # BayesDel
    dbnsfp   = data.get("dbnsfp", {})
    result["returned_transcripts"] = _returned_transcripts(dbnsfp)
    bayesdel = dbnsfp.get("bayesdel", {}) if isinstance(dbnsfp, dict) else {}
    bd_score = None
    if isinstance(bayesdel, dict):
        no_af = bayesdel.get("no_af", {})
        if isinstance(no_af, dict):
            bd_score = no_af.get("score")
    selected, basis, selection_error = _select_unambiguous_bayesdel(bd_score)
    result["bayesdel"] = selected
    result["selection_basis"] = basis

    # AlphaMissense (only available for missense variants)
    am = data.get("alphamissense", {})
    if isinstance(am, dict):
        am_score = am.get("am_pathogenicity")
        am_class = am.get("am_class")
        if am_score is not None:
            result["am_score"] = float(am_score)
            result["am_class"] = am_class or ""

    if selection_error:
        result["status"] = "ambiguous_transcript"
        result["error"] = (
            f"{selection_error}; configured reference transcript is "
            f"{result['reference_transcript']}. BayesDel was not used."
        )
    else:
        result["status"] = "ok" if (
            result["bayesdel"] is not None or result["am_score"] is not None
        ) else "no_score"
    return result


def get_bayesdel_and_alphamissense(
    gene: str, c_notation: str
) -> Tuple[Optional[float], Optional[dict]]:
    """
    Return (bayesdel_score, alphamissense_dict) from cache or API.
    alphamissense_dict is {"am_score": float, "am_class": str} or None.
    Single function so both scores come from one API call.
    """
    variant_key = f"{gene}:{c_notation}"

    entry = BAYESDEL_CACHE.get(variant_key)
    if (
        isinstance(entry, dict)
        and entry.get("bayesdel") is not None
        and entry.get("selection_policy") != BAYESDEL_SELECTION_POLICY
    ):
        BAYESDEL_CACHE.pop(variant_key, None)
        entry = None
    if isinstance(entry, dict) and entry.get("status") == "no_score":
        # Defensive migration for a process that obtained a transient empty
        # response before this policy was applied.  Retry instead of treating
        # it as a durable cache hit.
        BAYESDEL_CACHE.pop(variant_key, None)
        entry = None
    if entry is not None:
        bd = entry.get("bayesdel") if isinstance(entry, dict) else entry
        am_score = entry.get("am_score") if isinstance(entry, dict) else None
        am_class = entry.get("am_class") if isinstance(entry, dict) else None
        am = {"am_score": am_score, "am_class": am_class} if am_score is not None else None
        BAYESDEL_STATUS_CACHE[variant_key] = {
            "status": entry.get("status", "ok" if bd is not None or am_score is not None else "no_score"),
            "reason": entry.get("reason", "Loaded from local MyVariant/BayesDel cache"),
        }
        return bd, am

    rv = coordinates.resolve_variant(gene, c_notation)
    hg37 = None
    if rv and rv.has_grch37():
        hg37 = {"chrom": rv.grch37.chrom, "pos": rv.grch37.pos,
                "ref": rv.grch37.ref, "alt": rv.grch37.alt}

    data = fetch_variant_data_myvariant(gene, c_notation, hg37)
    BAYESDEL_STATUS_CACHE[variant_key] = {
        "status": data["status"],
        "reason": data.get("error") or {
            "no_grch37_coords": "No GRCh37 coordinates available for MyVariant",
            "not_found": "Variant was not found by MyVariant",
            "no_score": "MyVariant response contained no BayesDel or AlphaMissense score",
            "ambiguous_transcript": (
                "MyVariant returned multiple distinct BayesDel_noAF values without "
                "an explicit value-to-transcript mapping"
            ),
            "ok": "MyVariant response parsed successfully",
        }.get(data["status"], data["status"]),
    }

    # Only stable responses belong in the persistent cache. Coordinate and API
    # failures remain retryable instead of being converted to a silent null.
    if data["status"] in {"ok", "not_found", "ambiguous_transcript"}:
        BAYESDEL_CACHE[variant_key] = {
            "bayesdel": data["bayesdel"],
            "am_score": data["am_score"],
            "am_class": data["am_class"],
            "status": data["status"],
            "reason": BAYESDEL_STATUS_CACHE[variant_key]["reason"],
            "selection_policy": data.get("selection_policy"),
            "selection_basis": data.get("selection_basis"),
            "reference_transcript": data.get("reference_transcript"),
            "returned_transcripts": data.get("returned_transcripts", []),
        }
        _save_cache()

    am = ({"am_score": data["am_score"], "am_class": data["am_class"]}
          if data["am_score"] is not None else None)
    return data["bayesdel"], am


def get_bayesdel_score(gene: str, c_notation: str) -> Optional[float]:
    """Backward-compatible wrapper returning only the BayesDel score."""
    bd, _ = get_bayesdel_and_alphamissense(gene, c_notation)
    return bd


_load_cache()
