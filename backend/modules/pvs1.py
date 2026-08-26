from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse

from backend.modules.utils import (
    get_amino_acid_position,
    get_cds_position_from_c_notation,
    get_intron_offset_from_c_notation,
)
from backend.modules.table4 import (
    TABLE4_DATA,
    table4_lookup_splice,
    parse_pvs1_code_strength,
    table4_lookup_pvs1_ptc,
    table4_lookup_deletion,
    table4_lookup_duplication,
    parse_exon_from_deletion_notation,
    parse_exon_from_duplication_notation,
)
from backend.gene_policy import (
    decision_asset,
    pvs1_thresholds,
    spliceai_thresholds,
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
                "source_id": "enigma-v1.2-table4",
                "label": "ENIGMA Specifications Table 4 v1.2",
                "url": "https://cspec.genome.network/cspec/File/id/ca5cf57b-94df-4ad6-a001-c62ceccb3845/data",
                "location": "Exact gene/exon rule",
            },
        ],
    }


def _termination_aa_position(p_notation: str) -> Optional[int]:
    """Return the predicted stop position for consequence provenance."""
    text = (p_notation or "").replace("*", "Ter")
    frameshift = re.search(r"[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}fsTer(\d+)", text)
    if frameshift:
        return int(frameshift.group(1)) + int(frameshift.group(2)) - 1
    nonsense = re.search(r"[A-Z][a-z]{2}(\d+)Ter", text)
    return int(nonsense.group(1)) if nonsense else None

def evaluate_pvs1(
    gene: str,
    variant_type: str,
    p_notation: str,
    c_notation: str = "",
    spliceai_score: Optional[float] = None,
    dup_type: str = "Unknown",
) -> Dict:
    """
    Evaluate PVS1 using ENIGMA Table 4 decision tree.

    v1.6.0 fixes:
    - Uses parse_pvs1_code_strength() for proper RNA code handling
    - Critical boundary logic: aa <= boundary -> PVS1, aa > boundary -> PVS1_N/A
    """
    result = {
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
        "source": vcep_specification(gene)["url"],
    }

    lof_types = [
        "frameshift", "nonsense", "splice_site", "initiation_codon",
        "exon_deletion", "exon_duplication",
    ]
    if variant_type.lower() not in lof_types:
        result["reason"] = f"PVS1 not applicable for {variant_type} variants"
        return result

    # Get CDS position and AA position
    cds_pos = get_cds_position_from_c_notation(c_notation)
    first_altered_aa = get_amino_acid_position(p_notation)
    termination_aa = _termination_aa_position(p_notation)

    if variant_type.lower() == "initiation_codon":
        result["reason"] = (
            "Initiation codon variant recognized. Automated PVS1 is not applied: "
            "the ENIGMA initiation flowchart requires a curated Module 1 data rule."
        )
        return result

    # ---- Splice site variants ----
    if variant_type.lower() == "splice_site":
        # Try Table 4 lookup for exact variant
        table4_splice = table4_lookup_splice(gene, c_notation)

        if table4_splice["found"]:
            pvs1_code = table4_splice["pvs1_code"]
            result["pvs1_code"] = pvs1_code
            result["exon"] = table4_splice["exon"]

            # Parse the code properly
            strength, points, requires_rna = parse_pvs1_code_strength(pvs1_code)

            if strength is None:  # PVS1_N/A
                result["applies"] = False
                result["reason"] = table4_splice["reason"]
                if table4_splice["notes"]:
                    result["reason"] += f" ({table4_splice['notes']})"
            elif requires_rna:
                result["requires_rna"] = True
                result["reason"] = (
                    f"{table4_splice['reason']} - RNA confirmation required; "
                    "PVS1 (RNA) is outside the automated Module 1 scope"
                )
            else:
                result["applies"] = True
                result["strength"] = strength
                result["points"] = points
                result["requires_rna"] = requires_rna
                result["reason"] = table4_splice["reason"]
                if requires_rna:
                    result["reason"] += " (requires RNA confirmation)"

            if result["applies"]:
                result["decision_path"] = _pvs1_path(
                    gene=gene,
                    branch_id="canonical-splice-outcome",
                    outcome_node="splice-table4",
                    steps=[
                        {
                            "node_id": "splice-coding",
                            "question": "Predicted alteration affects coding sequence?",
                            "result": "yes",
                            "observed": table4_splice["reason"],
                        },
                        {
                            "node_id": "splice-consequence",
                            "question": "Predicted splice consequence?",
                            "result": "predicted_alteration",
                            "observed": str(table4_splice.get("notes") or "Table 4 splice consequence"),
                        },
                        {
                            "node_id": "splice-transcript",
                            "question": "Biologically relevant transcript and critical region assessment",
                            "result": "table4",
                            "observed": f"{gene} {table4_splice['exon']}: {pvs1_code}",
                        },
                    ],
                )

            return result

        # Fallback: parse intron offset and use generic rules
        intron_info = get_intron_offset_from_c_notation(c_notation)
        if intron_info is None:
            result["reason"] = "Could not parse intronic offset from c. notation"
            return result

        _, offset = intron_info
        is_canonical = abs(offset) <= 2

        if is_canonical:
            # SAFETY: Do NOT auto-apply PVS1 for canonical splice not in Table 4
            # BRCA has exceptions (e.g. c.8953+2T>C is PVS1_N/A due to functional GC splice)
            result["applies"] = False
            result["strength"] = None
            result["points"] = 0
            result["reason"] = (
                f"Canonical splice site (offset {offset:+d}) - "
                f"NOT FOUND in Table 4. Manual review required. "
                f"Do not auto-apply PVS1 for BRCA splice variants."
            )
            splice_low = spliceai_thresholds(gene)["bp4"]
            if spliceai_score is not None and spliceai_score < splice_low:
                result["applies"] = False
                result["strength"] = None
                result["points"] = 0
                result["reason"] = f"Canonical splice but SpliceAI {spliceai_score:.3f} < {splice_low} - flag for review"
        else:
            score_str = f"{spliceai_score:.3f}" if spliceai_score is not None else "N/A"
            result["reason"] = (
                f"Non-canonical splice (offset {offset:+d}), SpliceAI {score_str}. "
                "PVS1 requires RNA evidence; use PP3 for predictive splice evidence when applicable."
            )

        return result

    # ---- Frameshift and nonsense variants ----
    if variant_type.lower() in ["frameshift", "nonsense"]:
        if cds_pos is None:
            result["reason"] = f"Could not parse CDS position from {c_notation}"
            return result

        table4_result = table4_lookup_pvs1_ptc(
            gene, cds_pos, first_altered_aa, termination_aa
        )

        result["exon"] = table4_result["exon"]
        result["pvs1_code"] = table4_result.get("pvs1_code")
        result["reason"] = table4_result["reason"]

        if table4_result["pvs1_strength"] is None:  # PVS1_N/A
            result["applies"] = False
        elif table4_result.get("requires_rna"):
            result["applies"] = False
            result["requires_rna"] = True
            result["reason"] += " - PVS1 (RNA) is outside the automated Module 1 scope"
        else:
            result["applies"] = True
            result["strength"] = table4_result["pvs1_strength"]
            result["points"] = table4_result["pvs1_points"]
            result["requires_rna"] = table4_result.get("requires_rna", False)

        # Pass through PM5 info
        if table4_result["pm5_code"]:
            result["pm5_code"] = table4_result["pm5_code"]
            result["pm5_strength"] = table4_result["pm5_strength"]
            result["pm5_points"] = table4_result["pm5_points"]
            result["pm5_exon"] = table4_result.get("pm5_exon")

        if result["applies"]:
            nmd_boundary = pvs1_thresholds(gene)["nmd_boundary_c_first_not_predicted"]
            estimated_ptc_c = (
                (termination_aa * 3 - 2) if termination_aa is not None else cds_pos
            )
            nmd_predicted = estimated_ptc_c < nmd_boundary
            branch_steps: List[Dict[str, str]] = [
                {
                    "node_id": "pvs-ptc-nmd",
                    "question": "Is nonsense-mediated decay predicted?",
                    "result": "yes" if nmd_predicted else "no",
                    "observed": (
                        f"Predicted termination near c.{estimated_ptc_c}; "
                        f"gene-specific NMD boundary c.{nmd_boundary - 1}_c.{nmd_boundary}"
                    ),
                }
            ]
            if nmd_predicted:
                branch_steps.append({
                    "node_id": "pvs-ptc-transcript",
                    "question": "Variant present in a biologically relevant transcript?",
                    "result": "present",
                    "observed": f"Table 4 {gene} {table4_result['exon']} -> {table4_result['pvs1_code']}",
                })
            else:
                branch_steps.append({
                    "node_id": "pvs-ptc-critical",
                    "question": "Is the truncated or altered region critical to protein function?",
                    "result": "yes",
                    "observed": table4_result["reason"],
                })
            result["decision_path"] = _pvs1_path(
                gene=gene,
                branch_id="nonsense-frameshift",
                outcome_node="pvs-ptc-table4",
                steps=branch_steps,
            )

        return result

    # ---- Exon deletion ----
    if variant_type.lower() == "exon_deletion":
        # Try to parse exon from c_notation
        exon = parse_exon_from_deletion_notation(c_notation, gene)

        if exon:
            del_result = table4_lookup_deletion(gene, exon)
            if del_result["found"]:
                result["exon"] = exon
                result["pvs1_code"] = del_result.get("pvs1_code")
                if del_result["pvs1_strength"]:
                    result["applies"] = True
                    result["strength"] = del_result["pvs1_strength"]
                    result["points"] = del_result["pvs1_points"]
                    result["reason"] = del_result["reason"]
                else:
                    result["applies"] = False
                    result["reason"] = del_result["reason"]
                if result["applies"]:
                    exon_range = TABLE4_DATA.get("exon_ranges", {}).get(gene, {}).get(exon)
                    exon_length = (
                        exon_range[1] - exon_range[0] + 1
                        if isinstance(exon_range, list) and len(exon_range) == 2
                        else None
                    )
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
                            "observed": del_result["reason"],
                        },
                        {
                            "node_id": "pvs-del-frame",
                            "question": "Predicted consequence: PTC or in-frame?",
                            "result": "in_frame" if in_frame else "ptc",
                            "observed": (
                                f"Coding length {exon_length} nt" if exon_length is not None
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

        result["applies"] = False
        result["reason"] = (
            f"PVS1 was not applied: {c_notation} could not be mapped "
            "unambiguously to an ENIGMA Table 4 deletion row."
        )
        return result

    # ---- Exon duplication ----
    if variant_type.lower() == "exon_duplication":
        # Try to parse exon from c_notation
        exon = parse_exon_from_duplication_notation(c_notation, gene)

        if exon:
            # Use dup_type parameter - tandem vs unknown affects strength
            dup_result = table4_lookup_duplication(gene, exon, dup_type)
            if dup_result["found"]:
                result["exon"] = exon
                result["pvs1_code"] = dup_result.get("pvs1_code")
                if dup_result["pvs1_strength"]:
                    result["applies"] = True
                    result["strength"] = dup_result["pvs1_strength"]
                    result["points"] = dup_result["pvs1_points"]
                    result["reason"] = dup_result["reason"]
                else:
                    result["applies"] = False
                    result["reason"] = dup_result["reason"]
                if result["applies"]:
                    exon_range = TABLE4_DATA.get("exon_ranges", {}).get(gene, {}).get(exon)
                    exon_length = (
                        exon_range[1] - exon_range[0] + 1
                        if isinstance(exon_range, list) and len(exon_range) == 2
                        else None
                    )
                    in_frame = exon_length is not None and exon_length % 3 == 0
                    steps = [
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
                            "question": "Is the in-frame duplication contained within a critical domain?",
                            "result": "table4",
                            "observed": dup_result["reason"],
                        })
                    result["decision_path"] = _pvs1_path(
                        gene=gene,
                        branch_id="duplication",
                        outcome_node="pvs-dup-table4",
                        steps=steps,
                    )
                return result

        result["applies"] = False
        result["reason"] = (
            f"PVS1 was not applied: {c_notation} could not be mapped "
            "unambiguously to an ENIGMA Table 4 duplication row."
        )
        return result

    return result
