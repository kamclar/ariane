# ============================================================
# ENIGMA Table 9 - Load from JSON file
# ============================================================
# v1.6.0: Load Table 9 data from JSON file in data/ directory
#
# Table 9 contains ~4,300 variants with reviewed functional assay evidence

from typing import Optional, Dict
from pathlib import Path
from functools import lru_cache
import hashlib
import json

from backend.domain.lazy_mapping import LazyJsonObject
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Load Table 9 from JSON
TABLE9_JSON_PATH = PROJECT_ROOT / "data" / "enigma_table9.json"

@lru_cache(maxsize=1)
def load_table9_data() -> dict:
    """Load Table 9 explicitly or on first lookup, never during module import."""
    if not TABLE9_JSON_PATH.is_file():
        raise RuntimeError(f"Required ENIGMA Table 9 dataset is missing: {TABLE9_JSON_PATH}")
    try:
        value = json.loads(TABLE9_JSON_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Required ENIGMA Table 9 dataset cannot be loaded: {TABLE9_JSON_PATH}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise RuntimeError("Required ENIGMA Table 9 dataset must be a JSON object")
    return value


@lru_cache(maxsize=1)
def table9_checksum() -> str:
    load_table9_data()
    return hashlib.sha256(TABLE9_JSON_PATH.read_bytes()).hexdigest()


TABLE9_DATA = LazyJsonObject(load_table9_data)


def table9_protein_notation(gene: str, c_notation: str) -> Optional[str]:
    """Return the reviewed protein consequence for an exact Table 9 variant."""
    entry = TABLE9_DATA.get("variants", {}).get(f"{gene}:{c_notation}")
    if not isinstance(entry, dict):
        return None
    value = entry.get("p_notation")
    return str(value) if value else None


def table9_lookup_ps3_bs3(gene: str, c_notation: str) -> Dict:
    """
    Lookup PS3/BS3 functional evidence from Table 9.

    Args:
        gene: "BRCA1" or "BRCA2"
        c_notation: HGVS c. notation (e.g., "c.509G>A")

    Returns:
        Dict with code, strength, points, reason
    """
    result = {
        "code": None,
        "strength": None,
        "points": 0,
        "applies": False,
        "reason": "",
        "reviewed": False,
        "splice_result_published": None,
        "spliceai_prediction": None,
        "predicted_or_observed_splicing": None,
        "source_record": None,
    }

    key = f"{gene}:{c_notation}"

    if key in TABLE9_DATA["variants"]:
        entry = TABLE9_DATA["variants"][key]
        code = entry["code"]
        strength = entry["strength"]
        text = entry.get("text", "")

        result["reviewed"] = True
        result["splice_result_published"] = entry.get("splice_result_published")
        result["spliceai_prediction"] = entry.get("spliceai_prediction")
        result["predicted_or_observed_splicing"] = entry.get(
            "predicted_or_observed_splicing"
        )

        assay_results = [
            value
            for index in range(1, 5)
            if (value := entry.get(f"result_{index}"))
        ]
        result["source_record"] = {
            "gene": entry.get("gene"),
            "c_notation": entry.get("c_notation"),
            "p_notation": entry.get("p_notation"),
            "assigned_code": code,
            "code_weight": strength,
            "standardised_text": text,
            "splice_result_published": entry.get("splice_result_published"),
            "spliceai_prediction": entry.get("spliceai_prediction"),
            "predicted_or_observed_splicing": entry.get(
                "predicted_or_observed_splicing"
            ),
            "publication_count": entry.get("publication_count"),
            "assay_results": assay_results,
        }

        result["code"] = code
        result["strength"] = strength
        result["applies"] = code in {"PS3", "BS3"}
        result["reason"] = f"Table 9: {text}"

        if code == "PS3":
            if strength == "Strong":
                result["points"] = 4
            elif strength == "Moderate":
                result["points"] = 2
            elif strength == "Supporting":
                result["points"] = 1
        elif code == "BS3":
            if strength == "Strong":
                result["points"] = -4
            elif strength == "Moderate":
                result["points"] = -2
            elif strength == "Supporting":
                result["points"] = -1
    else:
        result["reason"] = f"No Table 9 entry for {key}"

    return result
