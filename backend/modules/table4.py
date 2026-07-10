# ============================================================
# ENIGMA Table 4 - Load from JSON file
# ============================================================
# v1.6.0: Load Table 4 data from JSON file in data/ directory

from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Load Table 4 from JSON
TABLE4_JSON_PATH = PROJECT_ROOT / "data" / "enigma_table4.json"

# Fallback: try current directory
if not TABLE4_JSON_PATH.exists():
    TABLE4_JSON_PATH = Path("enigma_table4.json")

if TABLE4_JSON_PATH.exists():
    with open(TABLE4_JSON_PATH, 'r') as f:
        TABLE4_DATA = json.load(f)
else:
    TABLE4_DATA = {
        "exon_ranges": {"BRCA1": {}, "BRCA2": {}},
        "ptc_rules": {"BRCA1": {}, "BRCA2": {}},
        "splice_rules": {"BRCA1": {}, "BRCA2": {}},
        "deletion_rules": {"BRCA1": {}, "BRCA2": {}},
        "duplication_rules": {"BRCA1": {}, "BRCA2": {}},
        "critical_boundaries": {
            "BRCA1": {"exon": "E23(24)", "boundary_aa": 1854},
            "BRCA2": {"exon": "E27", "boundary_aa": 3309}
        }
    }


# Validate Table 4 data - check that all exons in rules have ranges
def validate_table4_data():
    """Check that all exons referenced in rules have corresponding ranges."""
    issues = []
    for gene in ["BRCA1", "BRCA2"]:
        exon_ranges = set(TABLE4_DATA.get("exon_ranges", {}).get(gene, {}).keys())

        # Check ptc_rules
        ptc_exons = set(TABLE4_DATA.get("ptc_rules", {}).get(gene, {}).keys())
        missing_ptc = ptc_exons - exon_ranges
        if missing_ptc:
            issues.append(f"{gene} ptc_rules missing ranges: {missing_ptc}")

        # Check deletion_rules
        del_exons = set(TABLE4_DATA.get("deletion_rules", {}).get(gene, {}).keys())
        missing_del = del_exons - exon_ranges
        if missing_del:
            issues.append(f"{gene} deletion_rules missing ranges: {missing_del}")

        # Check duplication_rules
        dup_exons = set(TABLE4_DATA.get("duplication_rules", {}).get(gene, {}).keys())
        missing_dup = dup_exons - exon_ranges
        if missing_dup:
            issues.append(f"{gene} duplication_rules missing ranges: {missing_dup}")

    if issues:
        for issue in issues:
            print(f"WARNING: Table 4 data validation issue: {issue}")

if __name__ == "__main__":
    validate_table4_data()


def get_table4_exon(gene: str, cds_pos: int) -> Optional[str]:
    """Find the Table 4 exon for a given CDS position."""
    if gene not in TABLE4_DATA["exon_ranges"]:
        return None

    for exon_name, bounds in TABLE4_DATA["exon_ranges"][gene].items():
        start, end = bounds[0], bounds[1]
        if start <= cds_pos <= end:
            return exon_name

    return None


def table4_lookup_splice(gene: str, c_notation: str) -> Dict:
    """
    Lookup PVS1 code for splice site variant from Table 4.
    Table 4 has allele-specific rules for canonical splice variants.
    """
    result = {
        "found": False,
        "pvs1_code": None,
        "exon": None,
        "notes": "",
        "reason": ""
    }

    if gene not in TABLE4_DATA["splice_rules"]:
        result["reason"] = f"No Table 4 splice rules for {gene}"
        return result

    # Try different notation formats
    # Table 4 uses format like "c.8953+2TC" (no ">")
    formats_to_try = [
        c_notation,
        c_notation.replace(">", ""),
    ]

    for fmt in formats_to_try:
        if fmt in TABLE4_DATA["splice_rules"][gene]:
            entry = TABLE4_DATA["splice_rules"][gene][fmt]
            result["found"] = True
            result["pvs1_code"] = entry["pvs1_code"]
            result["exon"] = entry["exon"]
            result["notes"] = entry.get("notes", "")
            result["reason"] = f"Table 4 splice lookup: {gene} {c_notation} -> {entry['pvs1_code']}"
            return result

    result["reason"] = f"No Table 4 entry for {gene} {c_notation}"
    return result


def parse_pvs1_code_strength(pvs1_code: str) -> tuple:
    """
    Parse PVS1 code to get strength and points.

    Handles: PVS1, PVS1_N/A, PVS1_Strong, PVS1_Moderate, PVS1_Supporting
    And RNA variants: PVS1 (RNA), PVS1_Strong (RNA), etc.

    Returns: (strength, points, requires_rna)
    """
    if pvs1_code is None:
        return (None, 0, False)

    requires_rna = "(RNA)" in pvs1_code or "RNA" in pvs1_code

    # Remove (RNA) suffix for parsing
    code = pvs1_code.replace(" (RNA)", "").replace("(RNA)", "").strip()

    if code == "PVS1" or code == "PVS1_VeryStrong":
        return ("Very Strong", 8, requires_rna)
    elif code == "PVS1_Strong":
        return ("Strong", 4, requires_rna)
    elif code == "PVS1_Moderate":
        return ("Moderate", 2, requires_rna)
    elif code == "PVS1_Supporting":
        return ("Supporting", 1, requires_rna)
    elif code == "PVS1_N/A" or "N/A" in code:
        return (None, 0, requires_rna)
    else:
        # Unknown code, default to checking for keywords
        if "Strong" in code and "Very" not in code:
            return ("Strong", 4, requires_rna)
        elif "Moderate" in code:
            return ("Moderate", 2, requires_rna)
        elif "Supporting" in code:
            return ("Supporting", 1, requires_rna)
        else:
            return ("Very Strong", 8, requires_rna)


def table4_lookup_pvs1_ptc(
    gene: str,
    cds_pos: int,
    first_altered_aa: Optional[int] = None
) -> Dict:
    """
    Lookup PVS1/PM5 codes for PTC (nonsense/frameshift) variants using Table 4.

    CRITICAL FIX v1.6.0:
    For last exon with critical boundary:
    - aa <= boundary -> PVS1 (truncation removes critical region)
    - aa > boundary -> PVS1_N/A (truncation after critical region)
    """
    result = {
        "pvs1_code": None,
        "pvs1_strength": None,
        "pvs1_points": 0,
        "requires_rna": False,
        "pm5_code": None,
        "pm5_strength": None,
        "pm5_points": 0,
        "exon": None,
        "reason": ""
    }

    # Find exon
    exon = get_table4_exon(gene, cds_pos)
    if exon is None:
        result["reason"] = f"CDS position {cds_pos} not found in Table 4 exon ranges for {gene}"
        return result

    result["exon"] = exon

    # Get PTC rules for this exon
    if gene not in TABLE4_DATA["ptc_rules"] or exon not in TABLE4_DATA["ptc_rules"][gene]:
        result["reason"] = f"No Table 4 PTC rule for {gene} {exon}"
        return result

    rule = TABLE4_DATA["ptc_rules"][gene][exon]
    pvs1_code = rule["pvs1_code"]
    pm5_code = rule.get("pm5_code")

    # Check critical boundary for last exon
    # BRCA1 E23(24): PTC <= p.1854 -> PVS1, PTC > p.1854 -> PVS1_N/A
    # BRCA2 E27: PTC <= p.3309 -> PVS1, PTC > p.3309 -> PVS1_N/A
    if gene in TABLE4_DATA.get("critical_boundaries", {}):
        boundary_info = TABLE4_DATA["critical_boundaries"][gene]
        if exon == boundary_info["exon"]:
            # Safety check: if AA position unknown in critical boundary exon, require manual review
            if first_altered_aa is None:
                result["applies"] = False
                result["pvs1_code"] = None
                result["reason"] = (
                    f"Critical boundary exon {exon} but amino acid position could not be parsed from p_notation. "
                    f"Manual review required to determine if PTC is before or after p.{boundary_info['boundary_aa']}."
                )
                return result
            boundary_aa = boundary_info["boundary_aa"]

            if first_altered_aa <= boundary_aa:
                # Truncation AT OR BEFORE boundary = removes critical region = PVS1
                # PM5 code stays as specified in Table 4 for this exon (PM5_N/A for last exons
                # with critical boundary - E23(24) BRCA1 and E27 BRCA2)
                pvs1_code = "PVS1"
                result["reason"] = (
                    f"PTC at p.{first_altered_aa} <= p.{boundary_aa} (at/before critical boundary) "
                    f"in {exon} -> PVS1 (truncation removes critical C-terminal region)"
                )
            else:
                # Truncation AFTER boundary = critical region preserved = PVS1_N/A
                pvs1_code = "PVS1_N/A"
                pm5_code = "PM5_N/A"
                result["reason"] = (
                    f"PTC at p.{first_altered_aa} > p.{boundary_aa} (after critical boundary) "
                    f"in {exon} -> PVS1_N/A (critical C-terminal region preserved)"
                )
        else:
            result["reason"] = f"Table 4 lookup: {gene} {exon} PTC -> {pvs1_code}"
    else:
        result["reason"] = f"Table 4 lookup: {gene} {exon} PTC -> {pvs1_code}"

    # Parse PVS1 code to get strength and points
    strength, points, requires_rna = parse_pvs1_code_strength(pvs1_code)
    result["pvs1_code"] = pvs1_code
    result["pvs1_strength"] = strength
    result["pvs1_points"] = points
    result["requires_rna"] = requires_rna

    # Parse PM5 code
    # PM5_Strong (PTC), PM5_Supporting (PTC), PM5 (PTC), PM5_N/A
    if pm5_code and "N/A" not in pm5_code:
        result["pm5_code"] = pm5_code
        if "Strong" in pm5_code and "Very" not in pm5_code:
            result["pm5_strength"] = "Strong"
            result["pm5_points"] = 4
        elif "Supporting" in pm5_code:
            result["pm5_strength"] = "Supporting"
            result["pm5_points"] = 1
        elif pm5_code in ["PM5 (PTC)", "PM5", "PM5_Moderate"]:
            # Plain PM5 (PTC) = Moderate
            result["pm5_strength"] = "Moderate"
            result["pm5_points"] = 2

    return result


def table4_lookup_deletion(gene: str, exon: str) -> Dict:
    """Lookup PVS1 code for exon deletion."""
    result = {"found": False, "pvs1_code": None, "pvs1_strength": None, "pvs1_points": 0, "reason": ""}

    if gene in TABLE4_DATA["deletion_rules"] and exon in TABLE4_DATA["deletion_rules"][gene]:
        rule = TABLE4_DATA["deletion_rules"][gene][exon]
        result["found"] = True
        result["pvs1_code"] = rule["pvs1_code"]
        result["notes"] = rule.get("notes", "")
        result["reason"] = f"Table 4 deletion lookup: {gene} DEL {exon} -> {rule['pvs1_code']}"

        # Parse strength
        strength, points, _ = parse_pvs1_code_strength(rule["pvs1_code"])
        result["pvs1_strength"] = strength
        result["pvs1_points"] = points
    else:
        result["reason"] = f"No Table 4 entry for {gene} DEL {exon}"

    return result


def table4_lookup_duplication(gene: str, exon: str, dup_type: str = "Unknown") -> Dict:
    """
    Lookup PVS1 code for exon duplication.

    Args:
        gene: BRCA1 or BRCA2
        exon: Exon name (e.g., "E10(11)")
        dup_type: "Tandem" or "Unknown"

    Table 4 distinguishes:
        DUP Ex, Tandem -> usually PVS1_Moderate
        DUP Ex, Unknown -> usually PVS1_Supporting
    """
    result = {"found": False, "pvs1_code": None, "pvs1_strength": None, "pvs1_points": 0, "reason": ""}

    if gene not in TABLE4_DATA.get("duplication_rules", {}):
        result["reason"] = f"No duplication rules for {gene}"
        return result

    if exon not in TABLE4_DATA["duplication_rules"][gene]:
        result["reason"] = f"No Table 4 entry for {gene} DUP {exon}"
        return result

    exon_rules = TABLE4_DATA["duplication_rules"][gene][exon]

    # Try exact dup_type match, fallback to Unknown
    if dup_type in exon_rules:
        rule = exon_rules[dup_type]
    elif "Unknown" in exon_rules:
        rule = exon_rules["Unknown"]
        dup_type = "Unknown (fallback)"
    else:
        # Take first available
        dup_type = list(exon_rules.keys())[0]
        rule = exon_rules[dup_type]

    result["found"] = True
    result["pvs1_code"] = rule["pvs1_code"]
    result["dup_type"] = dup_type
    result["notes"] = rule.get("notes", "")
    result["reason"] = f"Table 4 duplication lookup: {gene} DUP {exon} ({dup_type}) -> {rule['pvs1_code']}"

    # Parse strength
    strength, points, _ = parse_pvs1_code_strength(rule["pvs1_code"])
    result["pvs1_strength"] = strength
    result["pvs1_points"] = points

    return result



def parse_exon_from_duplication_notation(c_notation: str, gene: str) -> Optional[str]:
    """
    Try to parse exon from complex duplication notation.

    Examples:
        c.(793+1_794-1)_(1909+1_1910-1)dup -> E10 (BRCA2)
        c.(670+1_671-1)_(4096+1_4097-1)dup -> E10(11) (BRCA1)

    Pattern: c.(X+1_Y-1)_(Z+1_W-1)dup
    """
    import re

    match = re.match(r'c\.\((\d+)\+1_(\d+)-1\)_\((\d+)\+1_(\d+)-1\)dup', c_notation)

    if match:
        exon_start = int(match.group(2))
        exon_end = int(match.group(3))

        if gene in TABLE4_DATA["exon_ranges"]:
            for exon_name, (start, end) in TABLE4_DATA["exon_ranges"][gene].items():
                if start == exon_start and end == exon_end:
                    return exon_name
                if abs(start - exon_start) <= 5 and abs(end - exon_end) <= 5:
                    return exon_name

    return None
def parse_exon_from_deletion_notation(c_notation: str, gene: str) -> Optional[str]:
    """
    Try to parse exon from complex deletion notation.

    Examples:
        c.(793+1_794-1)_(1909+1_1910-1)del -> E10 (BRCA2)
        c.(670+1_671-1)_(4096+1_4097-1)del -> E10(11) (BRCA1)

    Pattern: c.(X+1_Y-1)_(Z+1_W-1)del
    Where Y is exon start and Z is exon end
    """
    import re

    # Pattern for exon deletion: c.(X+1_Y-1)_(Z+1_W-1)del
    # Y = first coding position of deleted exon
    # Z = last coding position of deleted exon
    match = re.match(r'c\.\((\d+)\+1_(\d+)-1\)_\((\d+)\+1_(\d+)-1\)del', c_notation)

    if match:
        exon_start = int(match.group(2))  # Y = first coding pos
        exon_end = int(match.group(3))    # Z = last coding pos

        # Find matching exon in Table 4
        if gene in TABLE4_DATA["exon_ranges"]:
            for exon_name, (start, end) in TABLE4_DATA["exon_ranges"][gene].items():
                if start == exon_start and end == exon_end:
                    return exon_name
                # Also try approximate match (within 5bp)
                if abs(start - exon_start) <= 5 and abs(end - exon_end) <= 5:
                    return exon_name

    return None


if __name__ == "__main__":
    # Test Table 4 lookups
    print("\nTesting Table 4 lookups:")
    print("=" * 70)

    # Test critical boundary for BRCA1 c.5551_5552insT
    print("\n--- Critical boundary test (BRCA1 E23(24)) ---")
    test1 = table4_lookup_pvs1_ptc("BRCA1", 5551, first_altered_aa=1851)
    print(f"  c.5551 p.1851: exon={test1['exon']}, PVS1={test1['pvs1_code']} ({test1['pvs1_strength']})")
    print(f"    {test1['reason']}")

    test2 = table4_lookup_pvs1_ptc("BRCA1", 5565, first_altered_aa=1856)
    print(f"  c.5565 p.1856: exon={test2['exon']}, PVS1={test2['pvs1_code']} ({test2['pvs1_strength']})")
    print(f"    {test2['reason']}")

    # Test E9(10) which should be PVS1_N/A
    print("\n--- E9(10) test (should be PVS1_N/A) ---")
    test3 = table4_lookup_pvs1_ptc("BRCA1", 628, first_altered_aa=210)
    print(f"  c.628 p.210: exon={test3['exon']}, PVS1={test3['pvs1_code']}")

    # Test E10(11) which should now work
    print("\n--- E10(11) test (previously missing) ---")
    test4 = table4_lookup_pvs1_ptc("BRCA1", 3668, first_altered_aa=1225)
    print(f"  c.3668 p.1225: exon={test4['exon']}, PVS1={test4['pvs1_code']}")

    # Test splice lookup
    print("\n--- Splice lookup test ---")
    for var in ["c.8953+2T>C", "c.8953+2T>A", "c.8953+2T>G"]:
        r = table4_lookup_splice("BRCA2", var)
        print(f"  BRCA2 {var}: {r['pvs1_code'] if r['found'] else 'NOT FOUND'}")






