# ============================================================
# Clinically important residues - informational display
#
# Source: ENIGMA VCEP v1.2 Appendix Tables 3 and 4
# Known pathogenic missense residues in BRCA1/2 functional domains.
# This is NOT a classification criterion - it provides context
# for the user to understand whether the variant falls on a
# position with demonstrated clinical importance.
# ============================================================
from typing import Optional, Dict
from pathlib import Path
import json
import re
from backend.data_health import clear_issue, register_issue

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESIDUES_PATH = PROJECT_ROOT / "data" / "clinically_important_residues.json"

_RESIDUES_DATA = {}
_RESIDUES_LOADED = False


def _load_residues():
    global _RESIDUES_DATA, _RESIDUES_LOADED
    if _RESIDUES_LOADED:
        return
    if not RESIDUES_PATH.exists():
        print(f"Clinically important residues not found: {RESIDUES_PATH}")
        register_issue("Clinically important residues", f"required informational dataset is missing: {RESIDUES_PATH}")
        _RESIDUES_LOADED = True
        return
    try:
        with open(RESIDUES_PATH, encoding="utf-8") as f:
            raw = json.load(f)
    except (OSError, json.JSONDecodeError) as exc:
        register_issue(
            "Clinically important residues",
            f"could not load {RESIDUES_PATH}: {type(exc).__name__}: {exc}",
        )
        _RESIDUES_LOADED = True
        return
    _RESIDUES_DATA = raw.get("domains", {})
    _RESIDUES_LOADED = True
    clear_issue("Clinically important residues")
    print(f"Clinically important residues loaded from {RESIDUES_PATH}")


def _extract_aa_position(p_notation: str) -> Optional[int]:
    if not p_notation:
        return None
    clean = p_notation.replace("p.", "").replace("(", "").replace(")", "").strip()
    m = re.search(r'[A-Z][a-z]{2}(\d+)', clean)
    if m:
        return int(m.group(1))
    return None


def check_important_residue(gene: str, p_notation: str) -> Dict:
    """
    Check if variant falls on a known clinically important residue.
    Returns informational dict - not a classification criterion.
    """
    _load_residues()

    result = {
        "is_important_residue": False,
        "domain": None,
        "known_pathogenic_at_position": [],
        "domain_note": None,
        "message": "",
    }

    pos = _extract_aa_position(p_notation)
    if pos is None:
        return result

    gene_domains = _RESIDUES_DATA.get(gene, {})
    for domain_name, domain_data in gene_domains.items():
        aa_range = domain_data.get("range", [0, 0])
        if aa_range[0] <= pos <= aa_range[1]:
            result["domain"] = domain_name

            pathogenic = domain_data.get("pathogenic_residues", [])
            at_position = [r for r in pathogenic if r.get("pos") == pos]

            if at_position:
                result["is_important_residue"] = True
                result["known_pathogenic_at_position"] = at_position
                variants_str = ", ".join(r["aa"] for r in at_position)
                result["message"] = (
                    f"Position {pos} in {domain_name} domain has known pathogenic "
                    f"variant(s): {variants_str}"
                )
            else:
                result["message"] = (
                    f"Position {pos} is in {domain_name} domain "
                    f"(aa {aa_range[0]}-{aa_range[1]}) but no known pathogenic "
                    f"missense at this specific position"
                )

            note = domain_data.get("note")
            if note:
                result["domain_note"] = note

            return result

    result["message"] = f"Position {pos} is outside all clinically important domains"
    return result


def initialize_residue_data() -> None:
    """Load the informational dataset so health status is available at startup."""
    _load_residues()
