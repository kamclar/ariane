# ============================================================
# ENIGMA Table 9 - Load from JSON file
# ============================================================
# v1.6.0: Load Table 9 data from JSON file in data/ directory
#
# Table 9 contains ~4,300 variants with reviewed functional assay evidence

from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Load Table 9 from JSON
TABLE9_JSON_PATH = PROJECT_ROOT / "data" / "enigma_table9.json"

if not TABLE9_JSON_PATH.is_file():
    raise RuntimeError(f"Required ENIGMA Table 9 dataset is missing: {TABLE9_JSON_PATH}")
with open(TABLE9_JSON_PATH, encoding="utf-8") as f:
    TABLE9_DATA = json.load(f)


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

        result["code"] = code
        result["strength"] = strength
        result["applies"] = code in {"PS3", "BS3"}
        result["reason"] = f"Table 9: {text[:100]}..." if len(text) > 100 else f"Table 9: {text}"

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


if __name__ == "__main__":
    # Test Table 9 lookup
    print("\nTesting Table 9 PS3/BS3 lookup:")
    print("=" * 70)

    test_cases = [
        ("BRCA1", "c.509G>A"),
        ("BRCA1", "c.1534C>T"),
        ("BRCA1", "c.3891_3893del"),
        ("BRCA1", "c.5217T>A"),
        ("BRCA1", "c.5551_5552insT"),  # not in Table 9
    ]

    for gene, c_notation in test_cases:
        result = table9_lookup_ps3_bs3(gene, c_notation)
        if result["applies"]:
            print(f"  {gene} {c_notation}: {result['code']}_{result['strength']} ({result['points']:+d} points)")
        else:
            print(f"  {gene} {c_notation}: No Table 9 entry")
