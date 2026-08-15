"""Mechanism-aware ENIGMA evidence interaction and deduplication helpers."""

from typing import Any, Dict, Iterable, List


SPECIFICATIONS_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "11e62fec-23b0-4a3e-b2df-751855301746/data"
)
APPENDIX_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "9e6119dc-90b9-42b5-a3b7-1a2eb28b1b12/data"
)

BIOINFORMATIC_CODES = {"PP3", "BP1", "BP4", "BP7"}
PROTEIN_FUNCTION_CODES = {"PS3", "BS3"}


def interaction(
    *,
    status: str,
    mechanism: str,
    criteria: Iterable[str],
    retained: Iterable[str] = (),
    suppressed: Iterable[str] = (),
    reason: str,
    source: str,
    source_url: str,
    review_required: bool = False,
) -> Dict[str, Any]:
    return {
        "status": status,
        "mechanism": mechanism,
        "criteria": list(criteria),
        "retained": list(retained),
        "suppressed": list(suppressed),
        "reason": reason,
        "source": source,
        "source_url": source_url,
        "review_required": review_required,
    }


def automatic_functional_interactions(
    criteria: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Explain why Table 9 evidence does not blanket-suppress predictions."""
    functional = sorted(PROTEIN_FUNCTION_CODES & set(criteria))
    bioinformatic = sorted(BIOINFORMATIC_CODES & set(criteria))
    if not functional or not bioinformatic:
        return []

    functional_points = sum(criteria[code].get("points", 0) for code in functional)
    bioinformatic_points = sum(criteria[code].get("points", 0) for code in bioinformatic)
    contradictory = functional_points * bioinformatic_points < 0
    return [
        interaction(
            status="conflict" if contradictory else "info",
            mechanism="protein_function_and_bioinformatic_context",
            criteria=functional + bioinformatic,
            retained=functional + bioinformatic,
            reason=(
                "Calibrated PS3/BS3 evidence and relevant bioinformatic criteria "
                "are both retained. ENIGMA Figure 1C instructs retaining relevant "
                "bioinformatic codes when PVS1 is not applied. The evidence points "
                "in opposite directions and requires expert review."
                if contradictory
                else
                "Calibrated PS3/BS3 evidence and relevant bioinformatic criteria "
                "are both retained. ENIGMA Figure 1C instructs retaining relevant "
                "bioinformatic codes when PVS1 is not applied."
            ),
            source="ENIGMA v1.2 Figure 1C and Appendix E",
            source_url=SPECIFICATIONS_URL,
            review_required=contradictory,
        )
    ]


def pvs1_prediction_deduplication() -> Dict[str, Any]:
    return interaction(
        status="deduplicated",
        mechanism="splicing_or_loss_of_function",
        criteria=["PVS1", "PP3"],
        retained=["PVS1"],
        suppressed=["PP3"],
        reason=(
            "PP3 was not counted because PVS1 already captures the predicted "
            "loss-of-function mechanism for this variant."
        ),
        source="ENIGMA v1.2 Figure 1A, Figure 1B and Appendix E",
        source_url=SPECIFICATIONS_URL,
    )


def apply_automatic_rna_interactions(
    criteria: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Apply Figure 1B deduplication to automatically accepted PVS1 (RNA)."""
    if "PVS1_RNA" not in criteria:
        return []

    warnings: List[Dict[str, Any]] = []
    replaceable = BIOINFORMATIC_CODES | {"PS1", "PS1_SPLICE"}
    suppressed = sorted(code for code in replaceable if code in criteria)
    for code in suppressed:
        criteria.pop(code, None)
    if suppressed:
        warnings.append(
            interaction(
                status="deduplicated",
                mechanism="experimentally_confirmed_splicing",
                criteria=["PVS1_RNA"] + suppressed,
                retained=["PVS1_RNA"],
                suppressed=suppressed,
                reason=(
                    "PVS1 (RNA) from ENIGMA-curated mRNA evidence replaces "
                    "weaker bioinformatic or predictive evidence for the same "
                    "splicing consequence."
                ),
                source="ENIGMA v1.2 Figure 1B and Appendix E",
                source_url=APPENDIX_URL,
            )
        )

    functional = sorted(PROTEIN_FUNCTION_CODES & set(criteria))
    if functional:
        warnings.append(
            interaction(
                status="review_required",
                mechanism="splicing_and_protein_function",
                criteria=["PVS1_RNA"] + functional,
                retained=["PVS1_RNA"] + functional,
                reason=(
                    "PVS1 (RNA) and protein-functional PS3/BS3 evidence are "
                    "retained as potentially distinct mechanisms. Confirm assay "
                    "scope and independence of the protein result."
                ),
                source="ENIGMA v1.2 Figure 1C and Appendix E",
                source_url=APPENDIX_URL,
                review_required=True,
            )
        )
    return warnings


def apply_manual_rna_interactions(
    combined: Dict[str, Dict[str, Any]],
    applied_manual_codes: set[str],
) -> List[Dict[str, Any]]:
    """Apply the explicit Figure 1B hierarchy after manual criteria are added."""
    warnings: List[Dict[str, Any]] = []

    if "PVS1_RNA" in applied_manual_codes:
        replaceable = BIOINFORMATIC_CODES | {"PS1", "PS1_SPLICE"}
        suppressed = sorted(code for code in replaceable if code in combined)
        for code in suppressed:
            combined.pop(code, None)
        if suppressed:
            warnings.append(
                interaction(
                    status="deduplicated",
                    mechanism="experimentally_confirmed_splicing",
                    criteria=["PVS1_RNA"] + suppressed,
                    retained=["PVS1_RNA"],
                    suppressed=suppressed,
                    reason=(
                        "Accepted damaging mRNA evidence replaces weaker "
                        "bioinformatic or predictive evidence for the same "
                        "splicing consequence."
                    ),
                    source="ENIGMA v1.2 Figure 1B and Appendix E",
                    source_url=SPECIFICATIONS_URL,
                )
            )

        functional = sorted(PROTEIN_FUNCTION_CODES & set(combined))
        if functional:
            warnings.append(
                interaction(
                    status="review_required",
                    mechanism="splicing_and_protein_function",
                    criteria=["PVS1_RNA"] + functional,
                    retained=["PVS1_RNA"] + functional,
                    reason=(
                        "PVS1 (RNA) and protein-functional PS3/BS3 evidence are "
                        "retained as potentially distinct mechanisms. Confirm assay "
                        "scope and whether the protein result is independent of the "
                        "RNA effect."
                    ),
                    source="ENIGMA v1.2 Figure 1C and Appendix E",
                    source_url=APPENDIX_URL,
                    review_required=True,
                )
            )

    if "BP7_RNA" in applied_manual_codes:
        if "BP7" in combined:
            combined.pop("BP7")
            warnings.append(
                interaction(
                    status="deduplicated",
                    mechanism="experimentally_excluded_splicing_effect",
                    criteria=["BP7_RNA", "BP7"],
                    retained=["BP7_RNA"],
                    suppressed=["BP7"],
                    reason=(
                        "BP7 Supporting was replaced by BP7 Strong (RNA). "
                        "The accepted mRNA result is stronger evidence for the "
                        "same absence of a damaging transcript effect."
                    ),
                    source="ENIGMA v1.2 Figure 1B and Appendix E",
                    source_url=SPECIFICATIONS_URL,
                )
            )

        if "PP3" in combined:
            warnings.append(
                interaction(
                    status="conflict",
                    mechanism="splicing",
                    criteria=["BP7_RNA", "PP3"],
                    retained=["BP7_RNA", "PP3"],
                    reason=(
                        "Accepted mRNA evidence shows no damaging transcript effect "
                        "while PP3 predicts a splice effect. Figure 1B retains the "
                        "applicable bioinformatic code, so both are shown and expert "
                        "review is required."
                    ),
                    source="ENIGMA v1.2 Figure 1B and Appendix E",
                    source_url=SPECIFICATIONS_URL,
                    review_required=True,
                )
            )

        functional = sorted(PROTEIN_FUNCTION_CODES & set(combined))
        if functional:
            warnings.append(
                interaction(
                    status="review_required",
                    mechanism="rna_and_protein_function",
                    criteria=["BP7_RNA"] + functional,
                    retained=["BP7_RNA"] + functional,
                    reason=(
                        "BP7 Strong (RNA) co-occurs with protein-functional "
                        "evidence. Confirm that BP7 (RNA) eligibility and the "
                        "protein assay interpretation satisfy Figure 1B, Figure 1C "
                        "and Appendix E."
                    ),
                    source="ENIGMA v1.2 Figure 1B, Figure 1C and Appendix E",
                    source_url=APPENDIX_URL,
                    review_required=True,
                )
            )

    return warnings
