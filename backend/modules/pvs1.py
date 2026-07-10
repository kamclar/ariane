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
    table4_lookup_splice,
    parse_pvs1_code_strength,
    table4_lookup_pvs1_ptc,
    table4_lookup_deletion,
    table4_lookup_duplication,
    parse_exon_from_deletion_notation,
    parse_exon_from_duplication_notation,
)

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
        "exon": None
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
            if spliceai_score is not None and spliceai_score < 0.1:
                result["applies"] = False
                result["strength"] = None
                result["points"] = 0
                result["reason"] = f"Canonical splice but SpliceAI {spliceai_score:.3f} < 0.1 - flag for review"
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

        table4_result = table4_lookup_pvs1_ptc(gene, cds_pos, first_altered_aa)

        result["exon"] = table4_result["exon"]
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

        return result

    # ---- Exon deletion ----
    if variant_type.lower() == "exon_deletion":
        # Try to parse exon from c_notation
        exon = parse_exon_from_deletion_notation(c_notation, gene)

        if exon:
            del_result = table4_lookup_deletion(gene, exon)
            if del_result["found"]:
                result["exon"] = exon
                if del_result["pvs1_strength"]:
                    result["applies"] = True
                    result["strength"] = del_result["pvs1_strength"]
                    result["points"] = del_result["pvs1_points"]
                    result["reason"] = del_result["reason"]
                else:
                    result["applies"] = False
                    result["reason"] = del_result["reason"]
                return result

        result["applies"] = False
        result["reason"] = (
            f"Exon deletion - could not parse exon from {c_notation}. "
            f"Use table4_lookup_deletion(gene, exon) manually."
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
                if dup_result["pvs1_strength"]:
                    result["applies"] = True
                    result["strength"] = dup_result["pvs1_strength"]
                    result["points"] = dup_result["pvs1_points"]
                    result["reason"] = dup_result["reason"]
                else:
                    result["applies"] = False
                    result["reason"] = dup_result["reason"]
                return result

        result["applies"] = False
        result["reason"] = (
            f"Exon duplication - could not parse exon from {c_notation}. "
            f"Use table4_lookup_duplication(gene, exon, dup_type) manually."
        )
        return result

    return result


if __name__ == "__main__":
    # Test evaluate_pvs1
    print("\nTesting evaluate_pvs1:")
    print("=" * 70)

    test_cases = [
        ("BRCA1", "frameshift", "p.(Asp1851ValfsTer29)", "c.5551_5552insT"),  # Should be PVS1 Very Strong
        ("BRCA1", "nonsense", "p.(Gln210Ter)", "c.628C>T"),  # E9(10) -> PVS1_N/A
        ("BRCA1", "frameshift", "p.(Cys1225fs)", "c.3668_3671dup"),  # E10(11) -> PVS1
        ("BRCA2", "splice_site", "p.(?)", "c.8953+2T>C"),  # PVS1_N/A
        ("BRCA2", "splice_site", "p.(?)", "c.8953+2T>A"),  # PVS1
    ]

    for gene, vtype, p, c in test_cases:
        r = evaluate_pvs1(gene, vtype, p, c, spliceai_score=0.9)
        status = f"{r['strength']} ({r['points']} pts)" if r['applies'] else "Not applied"
        print(f"  {gene} {c}: PVS1 {status}")
        print(f"    Exon: {r.get('exon', 'N/A')}")
        print(f"    {r['reason'][:80]}...")
        if r.get('pm5_code'):
            print(f"    PM5: {r['pm5_code']} ({r['pm5_points']} pts)")




