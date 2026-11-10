from typing import Optional, Dict, List

from backend.variant_processing.utils import (
    get_amino_acid_position,
    get_cds_position_from_c_notation,
    get_intron_offset_from_c_notation,
    get_termination_amino_acid_position,
)
from backend.reference_data.table4 import (
    TABLE4_SOURCE_ID,
    parse_exon_from_deletion_notation,
    parse_exon_from_duplication_notation,
    parse_pvs1_code_strength,
    table4_exon_range,
    table4_lookup_deletion,
    table4_lookup_duplication,
    table4_lookup_pvs1_ptc,
    table4_lookup_splice,
    table4_source_url,
)
from backend.policy.gene import (
    decision_asset,
    pvs1_thresholds,
    vcep_specification,
)


APPENDIX_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "9e6119dc-90b9-42b5-a3b7-1a2eb28b1b12/data"
)


def _pvs1_path(
    *,
    gene: str,
    branch_id: str,
    outcome_node: str,
    steps: List[Dict[str, str]],
) -> Dict:
    specification = vcep_specification(gene)
    if branch_id == "canonical-splice-outcome":
        tree_id = "appendix-pvs1-splice"
        asset = decision_asset(gene, "PVS1", "splice")
    else:
        tree_id = "appendix-pvs1-nonsplice"
        asset = decision_asset(gene, "PVS1", "nonsplice")
    figure_number = asset["figure_number"]
    figure_url = asset["figure_url"]
    return {
        "tree_id": tree_id,
        "tree_version": "ENIGMA VCEP 1.2.0",
        "branch_id": branch_id,
        "criterion": "PVS1",
        "outcome": "applied",
        "outcome_node": outcome_node,
        "steps": steps,
        "sources": [
            {
                "source_id": "enigma-v1.2-appendix",
                "label": f"VCEP {specification['id']} Appendix v{specification['version']}",
                "url": APPENDIX_URL,
                "location": f"Figure {figure_number}",
                "figure_url": figure_url,
            },
            {
                "source_id": TABLE4_SOURCE_ID,
                "label": "ENIGMA Specifications Table 4 v1.2",
                "url": table4_source_url(),
                "location": "Exact gene/exon rule",
            },
        ],
    }


LOF_TYPES = frozenset({
    "frameshift",
    "nonsense",
    "ptc_delins",
    "splice_site",
    "initiation_codon",
    "exon_deletion",
    "exon_duplication",
})
PTC_TYPES = frozenset({"frameshift", "nonsense", "ptc_delins"})


def _empty_result(gene: str) -> Dict:
    return {
        "applies": False,
        "strength": None,
        "points": 0,
        "reason": "",
        "requires_rna": False,
        "pm5_code": None,
        "pm5_strength": None,
        "pm5_points": 0,
        "exon": None,
        "pm5_exon": None,
        "pvs1_code": None,
        "review_required": False,
        "source": vcep_specification(gene)["url"],
    }


def _evaluate_splice(
    result: Dict,
    gene: str,
    c_notation: str,
    spliceai_score: Optional[float],
) -> Dict:
    rule = table4_lookup_splice(gene, c_notation)
    if rule["found"]:
        pvs1_code = rule["pvs1_code"]
        strength, points, requires_rna = parse_pvs1_code_strength(pvs1_code)
        result["pvs1_code"] = pvs1_code
        result["exon"] = rule["exon"]
        if strength is None:
            result["reason"] = rule["reason"]
            if rule["notes"]:
                result["reason"] += f" ({rule['notes']})"
        elif requires_rna:
            result["requires_rna"] = True
            result["reason"] = (
                f"{rule['reason']} - the exact RNA-derived rule is evaluated "
                "by the PVS1_RNA branch"
            )
        else:
            result.update({
                "applies": True,
                "strength": strength,
                "points": points,
                "reason": rule["reason"],
            })
            result["decision_path"] = _pvs1_path(
                gene=gene,
                branch_id="canonical-splice-outcome",
                outcome_node="splice-table4",
                steps=[
                    {
                        "node_id": "splice-coding",
                        "question": "Predicted alteration affects coding sequence?",
                        "result": "yes",
                        "observed": rule["reason"],
                    },
                    {
                        "node_id": "splice-consequence",
                        "question": "Predicted splice consequence?",
                        "result": "predicted_alteration",
                        "observed": str(
                            rule.get("notes") or "Table 4 splice consequence"
                        ),
                    },
                    {
                        "node_id": "splice-transcript",
                        "question": (
                            "Biologically relevant transcript and critical "
                            "region assessment"
                        ),
                        "result": "table4",
                        "observed": f"{gene} {rule['exon']}: {pvs1_code}",
                    },
                ],
            )
        return result

    intron_info = get_intron_offset_from_c_notation(c_notation)
    if intron_info is None:
        result["reason"] = "Could not parse intronic offset from c. notation"
        return result

    _, offset = intron_info
    if abs(offset) <= 2:
        result["reason"] = (
            f"Canonical splice site (offset {offset:+d}) - NOT FOUND in Table 4. "
            "Manual review required. Do not auto-apply PVS1 for BRCA splice variants."
        )
    else:
        score = f"{spliceai_score:.3f}" if spliceai_score is not None else "N/A"
        result["reason"] = (
            f"Non-canonical splice (offset {offset:+d}), SpliceAI {score}. "
            "PVS1 requires RNA evidence; use PP3 for predictive splice evidence "
            "when applicable."
        )
    return result


def _ptc_decision_path(
    gene: str,
    termination_aa: int,
    table4_result: Dict,
) -> Dict:
    nmd_boundary = pvs1_thresholds(gene)["nmd_boundary_c_first_not_predicted"]
    predicted_codon_start = termination_aa * 3 - 2
    nmd_predicted = predicted_codon_start < nmd_boundary
    steps: List[Dict[str, str]] = [
        {
            "node_id": "pvs-ptc-nmd",
            "question": "Is nonsense-mediated decay predicted?",
            "result": "yes" if nmd_predicted else "no",
            "observed": (
                f"Predicted termination p.{termination_aa}, near "
                f"c.{predicted_codon_start}; gene-specific NMD boundary "
                f"c.{nmd_boundary - 1}_c.{nmd_boundary}"
            ),
        }
    ]
    if nmd_predicted:
        steps.append({
            "node_id": "pvs-ptc-transcript",
            "question": "Variant present in a biologically relevant transcript?",
            "result": "present",
            "observed": (
                f"Table 4 {gene} {table4_result['exon']} -> "
                f"{table4_result['pvs1_code']}"
            ),
        })
    else:
        steps.append({
            "node_id": "pvs-ptc-critical",
            "question": "Is the truncated or altered region critical to protein function?",
            "result": "yes",
            "observed": table4_result["reason"],
        })
    return _pvs1_path(
        gene=gene,
        branch_id="nonsense-frameshift",
        outcome_node="pvs-ptc-table4",
        steps=steps,
    )


def _evaluate_ptc(
    result: Dict,
    gene: str,
    c_notation: str,
    p_notation: str,
) -> Dict:
    cds_pos = get_cds_position_from_c_notation(c_notation)
    if cds_pos is None:
        result["reason"] = f"Could not parse CDS position from {c_notation}"
        return result

    termination_aa = get_termination_amino_acid_position(p_notation)
    table4_result = table4_lookup_pvs1_ptc(
        gene,
        cds_pos,
        get_amino_acid_position(p_notation),
        termination_aa,
    )
    result.update({
        "exon": table4_result["exon"],
        "pvs1_code": table4_result.get("pvs1_code"),
        "reason": table4_result["reason"],
        "review_required": table4_result.get("review_required", False),
        "application_status": table4_result.get(
            "application_status",
            "not_applied",
        ),
    })

    if table4_result["pvs1_strength"] is not None:
        result.update({
            "applies": True,
            "strength": table4_result["pvs1_strength"],
            "points": table4_result["pvs1_points"],
            "requires_rna": False,
        })

    if table4_result["pm5_code"]:
        result.update({
            "pm5_code": table4_result["pm5_code"],
            "pm5_strength": table4_result["pm5_strength"],
            "pm5_points": table4_result["pm5_points"],
            "pm5_exon": table4_result.get("pm5_exon"),
        })

    if result["applies"] and termination_aa is None:
        result["reason"] += (
            "; NMD decision path is unavailable because the predicted termination "
            "position could not be determined from the normalized protein consequence. "
            "The nucleotide-change position was not used as a substitute."
        )
    elif result["applies"]:
        assert termination_aa is not None
        result["decision_path"] = _ptc_decision_path(
            gene,
            termination_aa,
            table4_result,
        )
    return result


def _apply_exon_rule(result: Dict, exon: str, rule: Dict) -> None:
    result["exon"] = exon
    result["pvs1_code"] = rule.get("pvs1_code")
    result["reason"] = rule["reason"]
    if rule["pvs1_strength"]:
        result["applies"] = True
        result["strength"] = rule["pvs1_strength"]
        result["points"] = rule["pvs1_points"]


def _exon_length(gene: str, exon: str) -> Optional[int]:
    exon_range = table4_exon_range(gene, exon)
    if exon_range is None:
        return None
    return exon_range[1] - exon_range[0] + 1


def _evaluate_exon_deletion(
    result: Dict,
    gene: str,
    c_notation: str,
) -> Dict:
    exon = parse_exon_from_deletion_notation(c_notation, gene)
    rule = table4_lookup_deletion(gene, exon) if exon else None
    if exon and rule and rule["found"]:
        _apply_exon_rule(result, exon, rule)
        if result["applies"]:
            exon_length = _exon_length(gene, exon)
            in_frame = exon_length is not None and exon_length % 3 == 0
            steps: List[Dict[str, str]] = [
                {
                    "node_id": "pvs-del-full",
                    "question": "Full gene deletion?",
                    "result": "no",
                    "observed": f"Single-exon Table 4 lookup: {exon}",
                },
                {
                    "node_id": "pvs-del-domain",
                    "question": "Targets a gene-specific critical coding exon?",
                    "result": "assessed",
                    "observed": rule["reason"],
                },
                {
                    "node_id": "pvs-del-frame",
                    "question": "Predicted consequence: PTC or in-frame?",
                    "result": "in_frame" if in_frame else "ptc",
                    "observed": (
                        f"Coding length {exon_length} nt"
                        if exon_length is not None
                        else "Reading-frame status represented by Table 4"
                    ),
                },
            ]
            if in_frame:
                steps.append({
                    "node_id": "pvs-del-size",
                    "question": "How much coding sequence is removed?",
                    "result": "table4",
                    "observed": f"Exact {gene} {exon} deletion rule",
                })
            result["decision_path"] = _pvs1_path(
                gene=gene,
                branch_id="exon-deletion",
                outcome_node="pvs-del-table4",
                steps=steps,
            )
        return result

    result["reason"] = (
        f"PVS1 was not applied: {c_notation} could not be mapped "
        "unambiguously to an ENIGMA Table 4 deletion row."
    )
    return result


def _evaluate_exon_duplication(
    result: Dict,
    gene: str,
    c_notation: str,
    dup_type: str,
) -> Dict:
    exon = parse_exon_from_duplication_notation(c_notation, gene)
    rule = table4_lookup_duplication(gene, exon, dup_type) if exon else None
    if exon and rule and rule["found"]:
        _apply_exon_rule(result, exon, rule)
        if result["applies"]:
            exon_length = _exon_length(gene, exon)
            in_frame = exon_length is not None and exon_length % 3 == 0
            steps: List[Dict[str, str]] = [
                {
                    "node_id": "pvs-dup-tandem",
                    "question": "Tandem arrangement proven, presumed, or excluded?",
                    "result": "proven_or_presumed",
                    "observed": f"Duplication arrangement: {dup_type}",
                },
                {
                    "node_id": "pvs-dup-frame",
                    "question": "PTC or reading frame preserved?",
                    "result": "in_frame" if in_frame else "ptc",
                    "observed": (
                        f"Duplicated coding length {exon_length} nt"
                        if exon_length is not None
                        else "Reading-frame status represented by Table 4"
                    ),
                },
            ]
            if in_frame:
                steps.append({
                    "node_id": "pvs-dup-domain",
                    "question": (
                        "Is the in-frame duplication contained within a "
                        "critical domain?"
                    ),
                    "result": "table4",
                    "observed": rule["reason"],
                })
            result["decision_path"] = _pvs1_path(
                gene=gene,
                branch_id="duplication",
                outcome_node="pvs-dup-table4",
                steps=steps,
            )
        return result

    result["reason"] = (
        f"PVS1 was not applied: {c_notation} could not be mapped "
        "unambiguously to an ENIGMA Table 4 duplication row."
    )
    return result


def evaluate_pvs1(
    gene: str,
    variant_type: str,
    p_notation: str,
    c_notation: str = "",
    spliceai_score: Optional[float] = None,
    dup_type: str = "Unknown",
) -> Dict:
    """Evaluate the applicable ENIGMA Table 4 PVS1 branch."""
    result = _empty_result(gene)
    normalized_type = variant_type.lower()
    if normalized_type not in LOF_TYPES:
        result["reason"] = f"PVS1 not applicable for {variant_type} variants"
        return result
    if normalized_type == "initiation_codon":
        result["reason"] = (
            "Initiation codon variant recognized. Automated PVS1 is not applied: "
            "the ENIGMA initiation flowchart requires a curated Module 1 data rule."
        )
        return result
    if normalized_type == "splice_site":
        return _evaluate_splice(result, gene, c_notation, spliceai_score)
    if normalized_type in PTC_TYPES:
        return _evaluate_ptc(result, gene, c_notation, p_notation)
    if normalized_type == "exon_deletion":
        return _evaluate_exon_deletion(result, gene, c_notation)
    return _evaluate_exon_duplication(result, gene, c_notation, dup_type)
