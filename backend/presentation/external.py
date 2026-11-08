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

ARIANE_CLASS_LABELS = {
    1: "Benign",
    2: "Likely Benign",
    3: "Uncertain Significance",
    4: "Likely Pathogenic",
    5: "Pathogenic",
}


def _discordance_scope(ariane_class: int, erepo_class: int) -> str:
    pathogenic = {4, 5}
    benign = {1, 2}
    if {ariane_class, erepo_class} <= pathogenic or {ariane_class, erepo_class} <= benign:
        return "The results have the same evidence direction but differ in final class."
    if 3 in {ariane_class, erepo_class}:
        return "The difference crosses the VUS boundary."
    return "The results have opposite pathogenic and benign directions."


def erepo_reconciliation(predicted_class: int, local_erepo_result: Dict | None) -> Dict:
    """Compare ARIANE with the exact pinned ERepo assertion without changing points."""
    local = local_erepo_result or {}
    status = str(local.get("status") or "not_found")
    record = dict(local.get("record") or {})
    erepo_class_text = str(record.get("classification") or "")
    erepo_class = CLASSIFICATION_MAP.get(erepo_class_text)
    ariane_label = ARIANE_CLASS_LABELS.get(predicted_class, str(predicted_class))
    method_version = str(record.get("assertion_method_version") or "")
    published_date = str(record.get("published_date") or "")
    criteria = [str(value) for value in record.get("criteria_met", []) if value]

    common = {
        "status": "not_found",
        "title": "",
        "explanation": "",
        "requires_review": False,
        "published_date": published_date,
        "assertion_url": str(record.get("erepo_url") or ""),
        "criteria_met": criteria,
    }
    if status == "current_vcep_assertion" and erepo_class is not None:
        if predicted_class == erepo_class:
            common.update(
                status="current_agreement",
                title="Current ENIGMA VCEP assertion agrees with ARIANE",
                explanation=(
                    f"The exact ENIGMA BRCA1/2 VCEP v{method_version or '1.2.0'} "
                    f"assertion and the ARIANE automated result both classify this "
                    f"variant as {ariane_label}. The ERepo assertion is an external "
                    "confirmation and has not added criteria or points."
                ),
            )
            return common

        common.update(
            status="current_discordance",
            title="Current ENIGMA VCEP and ARIANE results differ",
            explanation=(
                f"The exact current ENIGMA BRCA1/2 VCEP v{method_version or '1.2.0'} "
                f"assertion classifies this variant as {erepo_class_text}, while "
                f"ARIANE's criterion-by-criterion automated result is {ariane_label}. "
                f"{_discordance_scope(predicted_class, erepo_class)} ARIANE has not "
                "overwritten either result. Review the applied criteria, source "
                "versions and evidence dates before using the classification."
            ),
            requires_review=True,
        )
        return common

    if status == "historical_vcep_assertion":
        common.update(
            status="historical_context",
            title="Historical ENIGMA expert-panel assertion available",
            explanation=(
                f"ERepo records an ENIGMA classification of {erepo_class_text or 'unknown'} "
                f"under specification version {method_version or 'unknown'}, not the active "
                "v1.2 method. It is shown as historical context and does not change the "
                "ARIANE result."
            ),
        )
        return common

    if status == "unversioned_vcep_assertion":
        common.update(
            status="unversioned_context",
            title="Unversioned ENIGMA expert-panel assertion available",
            explanation=(
                f"ERepo records an ENIGMA classification of {erepo_class_text or 'unknown'}, "
                "but the assertion method version is not recorded. It is shown as context "
                "and does not change the ARIANE result."
            ),
        )
    return common


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
