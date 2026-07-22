# ============================================================
# PP4/BP5 - Multifactorial likelihood clinical evidence
#
# Source: ENIGMA VCEP v1.2, Supplementary Table 7 (ST7)
# 773 reference set variants with posterior probability from
# multifactorial likelihood analysis.
#
# PP4/BP5 thresholds are based on combined clinical LR,
# converted from posterior probability using prior = 0.10
# (global prior probability of pathogenicity per ENIGMA).
#
# LR = posterior * (1 - prior) / (prior * (1 - posterior))
#
# PP4 (towards pathogenicity):
#   PP4 Supporting:   LR >= 2.08
#   PP4 Moderate:     LR >= 4.3
#   PP4 Strong:       LR >= 18.7
#   PP4 Very Strong:  LR >= 350
#
# BP5 (against pathogenicity):
#   BP5 Supporting:   LR <= 0.48
#   BP5 Moderate:     LR <= 0.23
#   BP5 Strong:       LR <= 0.05
#   BP5 Very Strong:  LR <= 0.00285
# ============================================================
from typing import Optional, Dict
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ST7_PATH = PROJECT_ROOT / "data" / "st7_reference_set.json"

# point values per Tavtigian 2020
PP4_POINTS = {
    "Very Strong": 8,
    "Strong":      4,
    "Moderate":    2,
    "Supporting":  1,
}

BP5_POINTS = {
    "Very Strong": -8,
    "Strong":      -4,
    "Moderate":    -2,
    "Supporting":  -1,
}

PRIOR = 0.10  # global prior probability of pathogenicity

# in-memory cache
_ST7_DATA = {}
_ST7_LOADED = False


def _load_st7():
    global _ST7_DATA, _ST7_LOADED
    if _ST7_LOADED:
        return

    if not ST7_PATH.exists():
        print(f"ST7 reference set not found: {ST7_PATH}")
        _ST7_LOADED = True
        return

    with open(ST7_PATH) as f:
        raw = json.load(f)

    for v in raw.get("variants", []):
        key = f"{v['gene']}:{v['c_notation']}"
        _ST7_DATA[key] = v

    _ST7_LOADED = True
    print(f"ST7 reference set loaded: {len(_ST7_DATA)} variants")


def posterior_to_lr(posterior: float) -> Optional[float]:
    """Convert posterior probability to likelihood ratio using global prior."""
    if posterior is None or posterior < 0 or posterior > 1:
        return None
    if posterior == 0:
        return 0.0
    if posterior == 1:
        return float("inf")
    lr = (posterior * (1 - PRIOR)) / (PRIOR * (1 - posterior))
    return lr


def lr_to_pp4_strength(lr: float) -> Optional[str]:
    """Map LR to PP4 strength per ENIGMA v1.2."""
    if lr >= 350:
        return "Very Strong"
    elif lr >= 18.7:
        return "Strong"
    elif lr >= 4.3:
        return "Moderate"
    elif lr >= 2.08:
        return "Supporting"
    return None


def lr_to_bp5_strength(lr: float) -> Optional[str]:
    """Map LR to BP5 strength per ENIGMA v1.2."""
    if lr <= 0.00285:
        return "Very Strong"
    elif lr <= 0.05:
        return "Strong"
    elif lr <= 0.23:
        return "Moderate"
    elif lr <= 0.48:
        return "Supporting"
    return None


def evaluate_pp4_bp5(gene: str, c_notation: str) -> Dict:
    """
    Look up variant in ST7 reference set and assign PP4 or BP5
    based on posterior probability from multifactorial likelihood analysis.
    """
    _load_st7()

    key = f"{gene}:{c_notation}"
    result = {
        "applies": False,
        "code": None,
        "strength": None,
        "points": 0,
        "reason": "",
        "posterior_probability": None,
        "likelihood_ratio": None,
    }

    entry = _ST7_DATA.get(key)
    if entry is None:
        result["reason"] = f"Not in ST7 reference set"
        return result

    posterior = entry.get("posterior_probability")
    if posterior is None:
        result["reason"] = f"In ST7 but no posterior probability available"
        return result

    result["posterior_probability"] = posterior
    lr = posterior_to_lr(posterior)
    result["likelihood_ratio"] = lr

    if lr is None:
        result["reason"] = f"Could not compute LR from posterior {posterior}"
        return result

    # check PP4 (pathogenic direction)
    pp4_strength = lr_to_pp4_strength(lr)
    if pp4_strength:
        result["applies"] = True
        result["code"] = "PP4"
        result["strength"] = pp4_strength
        result["points"] = PP4_POINTS[pp4_strength]
        result["reason"] = (
            f"ST7 multifactorial likelihood: posterior={posterior:.6f}, "
            f"LR={lr:.1f} - PP4_{pp4_strength} "
            f"(source: {entry.get('source', 'ENIGMA')})"
        )
        return result

    # check BP5 (benign direction)
    bp5_strength = lr_to_bp5_strength(lr)
    if bp5_strength:
        result["applies"] = True
        result["code"] = "BP5"
        result["strength"] = bp5_strength
        result["points"] = BP5_POINTS[bp5_strength]
        result["reason"] = (
            f"ST7 multifactorial likelihood: posterior={posterior:.6f}, "
            f"LR={lr:.4f} - BP5_{bp5_strength} "
            f"(source: {entry.get('source', 'ENIGMA')})"
        )
        return result

    # LR is in the uninformative range (0.48 < LR < 2.08)
    result["reason"] = (
        f"ST7 multifactorial likelihood: posterior={posterior:.6f}, "
        f"LR={lr:.2f} - uninformative range (0.48 < LR < 2.08), "
        f"neither PP4 nor BP5 applicable"
    )
    return result
