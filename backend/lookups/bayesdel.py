# BayesDel_noAF + AlphaMissense lookup via myvariant.info
#
# Single API call returns both scores:
#   dbnsfp.bayesdel.no_af.score  - BayesDel_noAF
#   alphamissense.am_pathogenicity / am_class  - AlphaMissense (missense only)
#
# AlphaMissense is informational only - not used for ENIGMA VCEP scoring.
#
# Persistent cache: backend/data/bayesdel_cache.json
# Cache entry format: {"bayesdel": float|null, "am_score": float|null, "am_class": str|null}
# Old float-only entries are migrated automatically on load.

from typing import Optional, Dict, Tuple
from pathlib import Path
import json
import threading
import urllib.request
import urllib.parse
import urllib.error

MYVARIANT_BASE_URL = "https://myvariant.info/v1/variant"

# Cache stores dicts: {"bayesdel": ..., "am_score": ..., "am_class": ...}
BAYESDEL_CACHE: Dict[str, Optional[dict]] = {}

_CACHE_PATH = Path(__file__).resolve().parent.parent / "data" / "bayesdel_cache.json"
_FILE_LOCK  = threading.Lock()


def _load_cache() -> None:
    global BAYESDEL_CACHE
    if not _CACHE_PATH.exists():
        return
    try:
        with open(_CACHE_PATH, encoding="utf-8") as fh:
            raw = json.load(fh)
        for key, val in raw.items():
            if isinstance(val, dict):
                BAYESDEL_CACHE[key] = val
            else:
                # Migrate old format (float or null) - AM score not available
                BAYESDEL_CACHE[key] = {
                    "bayesdel": float(val) if isinstance(val, (int, float)) else None,
                    "am_score": None,
                    "am_class": None,
                }
        print(f"Loaded BayesDel/AM cache: {len(BAYESDEL_CACHE)} entries")
    except Exception as exc:
        print(f"Warning: could not load BayesDel cache: {exc}")


def _save_cache() -> None:
    with _FILE_LOCK:
        try:
            with open(_CACHE_PATH, "w", encoding="utf-8") as fh:
                json.dump(BAYESDEL_CACHE, fh)
        except Exception as exc:
            print(f"Warning: could not save BayesDel cache: {exc}")


def fetch_variant_data_myvariant(gene: str, c_notation: str, hg37_coords: Optional[dict]) -> dict:
    """Fetch BayesDel + AlphaMissense from myvariant.info in one request."""
    result = {
        "bayesdel": None, "am_score": None, "am_class": None,
        "status": "not_queried", "error": None,
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
           f"?{urllib.parse.urlencode({'fields': 'dbnsfp.bayesdel,alphamissense'})}")
    req = urllib.request.Request(url, headers={"User-Agent": "BRCA-ACMG/1.8.0"})

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
    bayesdel = dbnsfp.get("bayesdel", {}) if isinstance(dbnsfp, dict) else {}
    bd_score = None
    if isinstance(bayesdel, dict):
        no_af = bayesdel.get("no_af", {})
        if isinstance(no_af, dict):
            bd_score = no_af.get("score")
    if isinstance(bd_score, list):
        valid = [s for s in bd_score if s is not None]
        bd_score = max(valid) if valid else None
    if bd_score is not None:
        result["bayesdel"] = float(bd_score)

    # AlphaMissense (only available for missense variants)
    am = data.get("alphamissense", {})
    if isinstance(am, dict):
        am_score = am.get("am_pathogenicity")
        am_class = am.get("am_class")
        if am_score is not None:
            result["am_score"] = float(am_score)
            result["am_class"] = am_class or ""

    result["status"] = "ok" if (result["bayesdel"] is not None or result["am_score"] is not None) else "no_score"
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
    if entry is not None:
        bd = entry.get("bayesdel") if isinstance(entry, dict) else entry
        am_score = entry.get("am_score") if isinstance(entry, dict) else None
        am_class = entry.get("am_class") if isinstance(entry, dict) else None
        am = {"am_score": am_score, "am_class": am_class} if am_score is not None else None
        return bd, am

    from backend.lookups.coordinates import resolve_variant
    rv  = resolve_variant(gene, c_notation)
    hg37 = None
    if rv and rv.has_grch37():
        hg37 = {"chrom": rv.grch37.chrom, "pos": rv.grch37.pos,
                "ref": rv.grch37.ref, "alt": rv.grch37.alt}

    data = fetch_variant_data_myvariant(gene, c_notation, hg37)

    BAYESDEL_CACHE[variant_key] = {
        "bayesdel": data["bayesdel"],
        "am_score": data["am_score"],
        "am_class": data["am_class"],
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
