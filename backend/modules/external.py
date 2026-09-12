# ============================================================
# Live external comparison - ClinVar + ClinGen ERepo
# Read-only reference. Curated PVS1 RNA uses a separate pinned local registry.
# ============================================================
from typing import Dict

CLASSIFICATION_MAP = {
    "Pathogenic":               5,
    "Likely pathogenic":        4,
    "Likely Pathogenic":        4,
    "Uncertain significance":   3,
    "Uncertain Significance":   3,
    "Likely benign":            2,
    "Likely Benign":            2,
    "Benign":                   1,
}


def external_comparison(
    gene: str,
    c_notation: str,
    predicted_class: int,
    clinvar_result: Dict,
    erepo_result: Dict,
    local_erepo_result: Dict | None = None,
) -> Dict:
    """
    Compare predicted class against ENIGMA EP classification
    from ClinVar or ClinGen ERepo.
    """
    cv = clinvar_result
    er = erepo_result
    local_er = local_erepo_result or {}

    # The checksum-validated local ERepo snapshot is the authoritative current
    # VCEP comparison. Live services remain visible as external context only.
    enigma_class_str = ""
    enigma_source    = ""
    local_record = local_er.get("record") or {}
    if local_er.get("status") == "current_vcep_assertion":
        enigma_class_str = local_record.get("classification", "")
        enigma_source = "Local ClinGen ERepo ENIGMA VCEP v1.2 snapshot"
    elif local_record:
        enigma_class_str = local_record.get("classification", "")
        enigma_source = "Local ClinGen ERepo non-current assertion (context only)"
    elif er.get("status") == "ok":
        enigma_class_str = er["classification"]
        enigma_source = "Live ClinGen ERepo comparison"
    elif cv.get("enigma_submission"):
        enigma_class_str = cv["enigma_submission"]["class"]
        enigma_source = "ClinVar expert-panel comparison"

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
