# ============================================================
# External comparison - ClinVar + ClinGen ERepo
# Read-only reference, does not feed into classification
# ============================================================
from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse

CLASSIFICATION_MAP = {
    "Pathogenic":               5,
    "Likely pathogenic":        4,
    "Uncertain significance":   3,
    "Likely benign":            2,
    "Benign":                   1,
}


def external_comparison(
    gene: str,
    c_notation: str,
    predicted_class: int,
    clinvar_result: Dict,
    erepo_result: Dict,
) -> Dict:
    """
    Compare predicted class against ENIGMA EP classification
    from ClinVar or ClinGen ERepo.
    """
    cv = clinvar_result
    er = erepo_result

    # ENIGMA EP class - prefer ClinVar (more coverage), fallback to ERepo
    enigma_class_str = ""
    enigma_source    = ""
    if cv.get("enigma_submission"):
        enigma_class_str = cv["enigma_submission"]["class"]
        enigma_source    = "ClinVar EP"
    elif er.get("status") == "ok":
        enigma_class_str = er["classification"]
        enigma_source    = "ClinGen ERepo"

    enigma_class_int = CLASSIFICATION_MAP.get(enigma_class_str)
    match = (
        predicted_class == enigma_class_int
        if enigma_class_int is not None else None
    )

    return {
        "clinvar":          cv,
        "erepo":            er,
        "enigma_class":     enigma_class_str,
        "enigma_source":    enigma_source,
        "enigma_class_int": enigma_class_int,
        "predicted_class":  predicted_class,
        "match":            match,
    }
