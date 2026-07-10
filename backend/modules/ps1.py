# ============================================================
# PS1 - Same amino acid change as known pathogenic variant
#
# Source: ENIGMA VCEP v1.2, Supplementary Table 7 (ST7)
# Reference set of known P/LP variants.
#
# PS1 applies when:
#   - A different nucleotide change leads to the same amino acid change
#     as a known P/LP variant
#   - Neither the variant under assessment nor the reference variant
#     has a predicted splice effect (SpliceAI <= 0.1)
#
# Strength:
#   PS1 Strong:    reference variant is Pathogenic (class 5)
#   PS1 Moderate:  reference variant is Likely Pathogenic (class 4)
#
# Only applies to missense variants.
# ============================================================
from typing import Optional, Dict, List
from pathlib import Path
import json
import re

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ST7_PATH = PROJECT_ROOT / "data" / "st7_reference_set.json"

PS1_POINTS = {
    "Strong":   4,
    "Moderate": 2,
}

# in-memory lookup: gene -> amino_acid_change -> [variants]
_PS1_LOOKUP = {}
_PS1_LOADED = False


def _extract_aa_change(p_notation: str) -> Optional[str]:
    """
    Extract amino acid change from p. notation.
    p.(Cys61Gly) -> Cys61Gly
    p.(His41Arg) -> His41Arg
    Returns None for non-missense (frameshift, synonymous, nonsense).
    """
    if not p_notation:
        return None
    # remove p. prefix and parentheses
    clean = p_notation.replace("p.", "").replace("(", "").replace(")", "").strip()
    if not clean:
        return None
    # skip non-missense
    if any(x in clean.lower() for x in ["fs", "ter", "del", "ins", "dup", "=", "?"]):
        return None
    # basic missense pattern: Xxx999Yyy
    m = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', clean)
    if m:
        return clean
    return None


def _extract_aa_position_and_ref(p_notation: str) -> Optional[tuple]:
    """
    Extract (position, ref_aa, alt_aa) from p. notation.
    p.(Cys61Gly) -> (61, "Cys", "Gly")
    """
    change = _extract_aa_change(p_notation)
    if not change:
        return None
    m = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', change)
    if m:
        return (int(m.group(2)), m.group(1), m.group(3))
    return None


def _load_ps1_reference():
    global _PS1_LOOKUP, _PS1_LOADED
    if _PS1_LOADED:
        return

    if not ST7_PATH.exists():
        print(f"ST7 reference set not found: {ST7_PATH}")
        _PS1_LOADED = True
        return

    with open(ST7_PATH) as f:
        raw = json.load(f)

    count = 0
    for v in raw.get("variants", []):
        # only P/LP variants
        if "Pathogenic" not in (v.get("reference_set") or ""):
            continue
        gene = v["gene"]
        p_not = v.get("p_notation", "")
        parsed = _extract_aa_position_and_ref(p_not)
        if not parsed:
            continue
        pos, ref_aa, alt_aa = parsed
        # key: gene + position + alt_aa (the amino acid it changes TO)
        aa_key = f"{ref_aa}{pos}{alt_aa}"

        if gene not in _PS1_LOOKUP:
            _PS1_LOOKUP[gene] = {}
        if aa_key not in _PS1_LOOKUP[gene]:
            _PS1_LOOKUP[gene][aa_key] = []

        _PS1_LOOKUP[gene][aa_key].append({
            "c_notation": v["c_notation"],
            "p_notation": p_not,
            "iarc_class": v.get("iarc_class"),
            "posterior":  v.get("posterior_probability"),
            "source":     v.get("source", ""),
        })
        count += 1

    _PS1_LOADED = True
    print(f"PS1 reference set loaded: {count} P/LP missense variants")
    for gene in _PS1_LOOKUP:
        print(f"  {gene}: {len(_PS1_LOOKUP[gene])} unique amino acid changes")


def evaluate_ps1(
    gene: str,
    c_notation: str,
    p_notation: str,
    variant_type: str,
    spliceai_score: Optional[float] = None,
) -> Dict:
    """
    Check if a different nucleotide change at the same amino acid position
    is known to be pathogenic (PS1).

    Requirements:
    - Variant must be missense
    - SpliceAI <= 0.1 for the variant under assessment (no splice effect)
    - Reference variant must exist with same amino acid change but different c. notation
    """
    _load_ps1_reference()

    result = {
        "applies": False,
        "code": "PS1",
        "strength": None,
        "points": 0,
        "reason": "",
        "reference_variant": None,
    }

    # only for missense
    if variant_type.lower() not in ("missense",):
        result["reason"] = "PS1 only applies to missense variants"
        return result

    # check splice effect - PS1 requires confirmed absence of a predicted splice effect
    if spliceai_score is None:
        result["reason"] = (
            "PS1 not applicable - SpliceAI score not available, "
            "cannot confirm absence of a predicted splice effect"
        )
        return result
    if spliceai_score > 0.1:
        result["reason"] = (
            f"PS1 not applicable - SpliceAI {spliceai_score:.3f} > 0.1 "
            f"suggests predicted splice effect"
        )
        return result

    # extract amino acid change
    parsed = _extract_aa_position_and_ref(p_notation)
    if not parsed:
        result["reason"] = f"Could not extract amino acid change from {p_notation}"
        return result

    pos, ref_aa, alt_aa = parsed
    aa_key = f"{ref_aa}{pos}{alt_aa}"

    # look up in reference set
    gene_lookup = _PS1_LOOKUP.get(gene, {})
    matches = gene_lookup.get(aa_key, [])

    if not matches:
        result["reason"] = f"No known P/LP variant with {aa_key} in {gene}"
        return result

    # filter out same nucleotide change (PS1 requires DIFFERENT c. notation)
    different = [m for m in matches if m["c_notation"] != c_notation]
    if not different:
        result["reason"] = (
            f"Same amino acid change {aa_key} found but same nucleotide change "
            f"- PS1 requires different nucleotide change"
        )
        return result

    # use the strongest reference (class 5 > class 4)
    best = max(different, key=lambda m: m.get("iarc_class") or 0)

    if best.get("iarc_class") == 5:
        result["strength"] = "Strong"
    elif best.get("iarc_class") == 4:
        result["strength"] = "Moderate"
    else:
        result["strength"] = "Supporting"

    result["applies"] = True
    result["points"] = PS1_POINTS.get(result["strength"], 1)
    result["reference_variant"] = best
    result["reason"] = (
        f"Same amino acid change {aa_key} as known "
        f"{'Pathogenic' if best.get('iarc_class') == 5 else 'Likely Pathogenic'} "
        f"variant {best['c_notation']} {best['p_notation']} "
        f"(class {best.get('iarc_class')}, source: {best.get('source', 'ENIGMA')})"
    )
    return result
