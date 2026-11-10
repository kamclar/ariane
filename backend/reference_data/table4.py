"""Read and query the pinned ENIGMA Table 4 snapshot."""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import re
from typing import Any, Dict, NamedTuple, Optional

from backend.reference_data.paths import TABLE4_PATH


TABLE4_SOURCE_ID = "enigma-v1.2-table4"


class CriterionWeight(NamedTuple):
    strength: Optional[str]
    points: int
    requires_rna: bool = False


class PtcSelection(NamedTuple):
    pvs1_code: Optional[str]
    pm5_code: Optional[str]
    reason: str
    review_required: bool


PVS1_WEIGHTS = {
    "PVS1": ("Very Strong", 8),
    "PVS1_VeryStrong": ("Very Strong", 8),
    "PVS1_Strong": ("Strong", 4),
    "PVS1_Moderate": ("Moderate", 2),
    "PVS1_Supporting": ("Supporting", 1),
}

PM5_WEIGHTS = {
    "PM5_Strong (PTC)": ("Strong", 4),
    "PM5 (PTC)": ("Moderate", 2),
    "PM5": ("Moderate", 2),
    "PM5_Moderate": ("Moderate", 2),
    "PM5_Supporting (PTC)": ("Supporting", 1),
}


@lru_cache(maxsize=1)
def load_table4_data() -> dict[str, Any]:
    """Load Table 4 explicitly or on first lookup, never during module import."""
    if not TABLE4_PATH.is_file():
        raise RuntimeError(f"Required ENIGMA Table 4 dataset is missing: {TABLE4_PATH}")
    try:
        value = json.loads(TABLE4_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Required ENIGMA Table 4 dataset cannot be loaded: {TABLE4_PATH}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise RuntimeError("Required ENIGMA Table 4 dataset must be a JSON object")
    return value


def get_table4_exon(gene: str, cds_pos: int) -> Optional[str]:
    """Find the Table 4 exon for a given CDS position."""
    ranges = load_table4_data().get("exon_ranges", {}).get(gene, {})
    for exon_name, bounds in ranges.items():
        start, end = bounds[0], bounds[1]
        if start <= cds_pos <= end:
            return str(exon_name)

    return None


def table4_exon_range(gene: str, exon: str) -> Optional[tuple[int, int]]:
    """Return the exact coding range stored for one Table 4 exon."""
    bounds = load_table4_data().get("exon_ranges", {}).get(gene, {}).get(exon)
    if not isinstance(bounds, list) or len(bounds) != 2:
        return None
    return int(bounds[0]), int(bounds[1])


def table4_deletion_exons(gene: str) -> tuple[str, ...]:
    """Return exon labels that have an exact Table 4 deletion rule."""
    rules = load_table4_data().get("deletion_rules", {}).get(gene, {})
    return tuple(rules)


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

    gene_rules = load_table4_data().get("splice_rules", {}).get(gene)
    if not isinstance(gene_rules, dict):
        result["reason"] = f"No Table 4 splice rules for {gene}"
        return result

    formats_to_try = [
        c_notation,
        c_notation.replace(">", ""),
    ]

    for fmt in formats_to_try:
        if fmt in gene_rules:
            entry = gene_rules[fmt]
            result["found"] = True
            result["pvs1_code"] = entry["pvs1_code"]
            result["exon"] = entry["exon"]
            result["notes"] = entry.get("notes", "")
            result["reason"] = f"Table 4 splice lookup: {gene} {c_notation} -> {entry['pvs1_code']}"
            return result

    result["reason"] = f"No Table 4 entry for {gene} {c_notation}"
    return result


@lru_cache(maxsize=1)
def table4_checksum() -> str:
    """Return the checksum of the exact Table 4 snapshot used at runtime."""
    return hashlib.sha256(TABLE4_PATH.read_bytes()).hexdigest()


def table4_source_url() -> str:
    """Return the source URL recorded in the pinned Table 4 snapshot."""
    return str(load_table4_data().get("source_url") or "")


def lookup_table4_pvs1_rna(gene: str, c_notation: str) -> Dict:
    """Return an exact allele-specific Table 4 RNA decision, if present.

    Table 4 uses the ``(RNA)`` suffix only where the VCEP code is based on
    observed molecular impact at transcript level.  The source code therefore
    already contains the final allele-specific strength; this lookup exposes
    that reviewed fact without re-interpreting individual publications.
    """
    data = load_table4_data()
    source_url = str(data.get("source_url") or "")
    provenance = {
        "source_id": TABLE4_SOURCE_ID,
        "source_version": str(data.get("version") or ""),
        "source_checksum": table4_checksum(),
        "source_url": source_url,
    }
    result = table4_lookup_splice(gene, c_notation)
    if not result.get("found"):
        return {
            "status": "not_in_source",
            "record": None,
            "reason": str(result.get("reason") or "No exact Table 4 splice rule"),
            **provenance,
        }

    source_code = str(result.get("pvs1_code") or "")
    strength, points, requires_rna = parse_pvs1_code_strength(source_code)
    if not requires_rna:
        return {
            "status": "not_applicable",
            "record": None,
            "reason": (
                "The exact Table 4 rule is predictive rather than an observed "
                "RNA decision."
            ),
            **provenance,
        }

    record = {
        "gene": gene,
        "c_notation": c_notation,
        "exon": result.get("exon"),
        "code": "PVS1_RNA",
        "source_code": source_code,
        "strength": strength,
        "points": points,
        "notes": str(result.get("notes") or ""),
        "evidence_mechanism": "rna_splicing",
        "decision_scope": "exact_allele",
    }
    if strength is None:
        return {
            "status": "not_applicable_rna",
            "record": record,
            "reason": (
                f"ENIGMA Table 4 assigns {source_code} to the exact allele "
                f"{gene} {c_notation}."
            ),
            **provenance,
        }
    return {
        "status": "eligible",
        "record": record,
        "reason": (
            f"ENIGMA Table 4 assigns {source_code} to the exact allele "
            f"{gene} {c_notation} based on observed transcript impact."
        ),
        **provenance,
    }


def parse_pvs1_code_strength(pvs1_code: Optional[str]) -> CriterionWeight:
    """Translate one Table 4 PVS1 code into strength and points."""
    if pvs1_code is None:
        return CriterionWeight(None, 0)

    requires_rna = "(RNA)" in pvs1_code
    code = pvs1_code.replace(" (RNA)", "").replace("(RNA)", "").strip()
    strength, points = PVS1_WEIGHTS.get(code, (None, 0))
    return CriterionWeight(strength, points, requires_rna)


def parse_pm5_code_strength(pm5_code: Optional[str]) -> CriterionWeight:
    """Translate one Table 4 PM5 PTC code into strength and points."""
    strength, points = PM5_WEIGHTS.get(str(pm5_code), (None, 0))
    return CriterionWeight(strength, points)


def table4_lookup_initiation(gene: str) -> Dict:
    """Return the gene-level ENIGMA Table 4 initiation-codon rule."""
    entry = load_table4_data().get("initiation_rules", {}).get(gene)
    rule = entry if isinstance(entry, dict) else {}
    result = {
        "found": False,
        "pvs1_code": rule.get("pvs1_code"),
        "pvs1_strength": None,
        "pvs1_points": 0,
        "exon": rule.get("exon"),
        "notes": rule.get("notes", ""),
        "source_row": rule.get("source_row"),
        "alternative_start_assessed": rule.get("alternative_start_assessed"),
        "alternative_start_conclusion": rule.get("alternative_start_conclusion"),
        "nearest_alternative_start": rule.get("nearest_alternative_start"),
        "upstream_pathogenic_evidence": rule.get("upstream_pathogenic_evidence"),
        "functional_domain_impact": rule.get("functional_domain_impact"),
        "appendix_figure": rule.get("appendix_figure"),
        "reason": f"No Table 4 initiation-codon rule for {gene}",
    }
    if not rule:
        return result

    pvs1_code = rule.get("pvs1_code")
    strength, points, requires_rna = parse_pvs1_code_strength(pvs1_code)
    if requires_rna or strength is None:
        result["reason"] = (
            f"Table 4 initiation-codon rule for {gene} has an unsupported "
            f"code: {pvs1_code!r}"
        )
        return result

    result["found"] = True
    result["pvs1_strength"] = strength
    result["pvs1_points"] = points
    result["reason"] = (
        f"ENIGMA Table 4 assigns {pvs1_code} to {gene} initiation-codon variants"
    )
    return result


def _select_ptc_codes(
    data: dict[str, Any],
    gene: str,
    exon: str,
    default_pvs1_code: Optional[str],
    default_pm5_code: Optional[str],
    first_altered_aa: Optional[int],
    termination_aa: Optional[int],
) -> PtcSelection:
    boundary = data.get("critical_boundaries", {}).get(gene)
    if not isinstance(boundary, dict) or exon != boundary.get("exon"):
        return PtcSelection(
            default_pvs1_code,
            default_pm5_code,
            f"Table 4 lookup: {gene} {exon} PTC -> {default_pvs1_code}",
            False,
        )

    boundary_aa = int(boundary["boundary_aa"])
    missing_first_position = first_altered_aa is None
    missing_required_stop = (
        first_altered_aa is not None
        and first_altered_aa <= boundary_aa
        and termination_aa is None
    )
    if missing_first_position or missing_required_stop:
        missing_value = (
            "the first altered amino acid"
            if missing_first_position
            else "the introduced termination position"
        )
        return PtcSelection(
            None,
            None,
            (
                f"Critical boundary exon {exon}, but {missing_value} could "
                f"not be determined. Manual review is required to compare "
                f"the complete protein consequence with p.{boundary_aa}."
            ),
            True,
        )

    assert first_altered_aa is not None
    if (
        first_altered_aa <= boundary_aa
        and termination_aa is not None
        and termination_aa <= boundary_aa
    ):
        branch = boundary["at_or_before"]
        pvs1_code = branch.get("pvs1_code")
        pm5_code = branch.get("pm5_code")
        reason = (
            f"First altered residue p.{first_altered_aa} <= p.{boundary_aa} "
            f"and predicted termination p.{termination_aa} <= p.{boundary_aa} "
            f"in {exon} -> {pvs1_code}, {pm5_code}"
        )
        return PtcSelection(pvs1_code, pm5_code, reason, False)

    branch = boundary["after"]
    pvs1_code = branch.get("pvs1_code")
    pm5_code = branch.get("pm5_code")
    if first_altered_aa > boundary_aa:
        boundary_reason = (
            f"first altered residue p.{first_altered_aa} is after p.{boundary_aa}"
        )
    else:
        boundary_reason = (
            f"predicted termination p.{termination_aa} is after p.{boundary_aa}"
        )
    reason = (
        f"Last-exon PTC in {exon}: {boundary_reason}. ENIGMA requires both "
        "the first altered residue and introduced termination to be at or "
        f"before the critical boundary -> {pvs1_code}, {pm5_code}"
    )
    return PtcSelection(pvs1_code, pm5_code, reason, False)


def table4_lookup_pvs1_ptc(
    gene: str,
    cds_pos: int,
    first_altered_aa: Optional[int] = None,
    termination_aa: Optional[int] = None,
) -> Dict:
    """
    Lookup PVS1/PM5 codes for PTC variants using Table 4.

    PVS1 and PM5 (PTC) are taken from the same Table 4 row. For frameshifts and
    coding delins that introduce a stop codon, Appendix D assigns the PM5
    weight by the exon containing the nucleotide change, even when the
    predicted termination codon lies in a later exon.
    """
    result: Dict = {
        "applies": False,
        "review_required": False,
        "application_status": "not_applied",
        "pvs1_code": None,
        "pvs1_strength": None,
        "pvs1_points": 0,
        "requires_rna": False,
        "pm5_code": None,
        "pm5_strength": None,
        "pm5_points": 0,
        "exon": None,
        "pm5_exon": None,
        "notes": "",
        "reason": ""
    }

    data = load_table4_data()
    exon = get_table4_exon(gene, cds_pos)
    if exon is None:
        result["reason"] = f"CDS position {cds_pos} not found in Table 4 exon ranges for {gene}"
        return result

    result["exon"] = exon

    gene_rules = data.get("ptc_rules", {}).get(gene, {})
    if exon not in gene_rules:
        result["reason"] = f"No Table 4 PTC rule for {gene} {exon}"
        return result

    rule = gene_rules[exon]
    result["pm5_exon"] = exon
    result["notes"] = str(rule.get("notes") or "")
    selection = _select_ptc_codes(
        data,
        gene,
        exon,
        rule.get("pvs1_code"),
        rule.get("pm5_code"),
        first_altered_aa,
        termination_aa,
    )
    result["reason"] = selection.reason
    if selection.review_required:
        result["review_required"] = True
        result["application_status"] = "review_required"
        return result

    pvs1_code = selection.pvs1_code
    pm5_code = selection.pm5_code

    strength, points, requires_rna = parse_pvs1_code_strength(pvs1_code)
    result["pvs1_code"] = pvs1_code
    result["pvs1_strength"] = strength
    result["pvs1_points"] = points
    result["requires_rna"] = requires_rna
    result["applies"] = strength is not None
    result["application_status"] = (
        "applied" if strength is not None else "not_applicable"
    )

    if (
        termination_aa is not None
        and termination_aa != first_altered_aa
        and "N/A" not in str(pvs1_code)
    ):
        result["reason"] += (
            f"; predicted termination p.{termination_aa} recorded for the protein "
            f"consequence, but PM5 follows the selected Table 4 row for the "
            f"nucleotide-change exon {exon}"
        )

    if result["notes"]:
        result["reason"] += f"; Table 4 note: {result['notes']}"

    pm5_weight = parse_pm5_code_strength(pm5_code)
    if pm5_weight.strength is not None:
        result["pm5_code"] = pm5_code
        result["pm5_strength"] = pm5_weight.strength
        result["pm5_points"] = pm5_weight.points

    return result


def table4_lookup_deletion(gene: str, exon: str) -> Dict:
    """Lookup PVS1 code for exon deletion."""
    result = {"found": False, "pvs1_code": None, "pvs1_strength": None, "pvs1_points": 0, "reason": ""}

    gene_rules = load_table4_data().get("deletion_rules", {}).get(gene, {})
    if exon in gene_rules:
        rule = gene_rules[exon]
        result["found"] = True
        result["pvs1_code"] = rule["pvs1_code"]
        result["notes"] = rule.get("notes", "")
        if rule["pvs1_code"] == "PVS1_N/A":
            result["reason"] = (
                f"PVS1 was not applied: ENIGMA Table 4 marks the {gene} "
                f"{exon} deletion as PVS1 N/A."
            )
        else:
            result["reason"] = (
                f"ENIGMA Table 4 assigns {rule['pvs1_code']} to the "
                f"{gene} {exon} deletion."
            )

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

    gene_rules = load_table4_data().get("duplication_rules", {}).get(gene)
    if not isinstance(gene_rules, dict):
        result["reason"] = f"No duplication rules for {gene}"
        return result

    if exon not in gene_rules:
        result["reason"] = f"No Table 4 entry for {gene} DUP {exon}"
        return result

    exon_rules = gene_rules[exon]

    # Arrangement is evidence. Never substitute another arrangement rule.
    if dup_type in exon_rules:
        rule = exon_rules[dup_type]
    else:
        result["reason"] = (
            f"No exact Table 4 duplication rule for {gene} DUP {exon} "
            f"with arrangement {dup_type}; available arrangements: {sorted(exon_rules)}"
        )
        return result

    result["found"] = True
    result["pvs1_code"] = rule["pvs1_code"]
    result["dup_type"] = dup_type
    result["notes"] = rule.get("notes", "")
    result["reason"] = f"Table 4 duplication lookup: {gene} DUP {exon} ({dup_type}) -> {rule['pvs1_code']}"

    strength, points, _ = parse_pvs1_code_strength(rule["pvs1_code"])
    result["pvs1_strength"] = strength
    result["pvs1_points"] = points

    return result


def _parse_exon_from_structural_notation(
    c_notation: str,
    gene: str,
    operation: str,
) -> Optional[str]:
    pattern = rf"c\.\((\d+)\+1_(\d+)-1\)_\((\d+)\+1_(\d+)-1\){operation}$"
    match = re.fullmatch(pattern, c_notation)
    if match is None:
        return None

    exon_start = int(match.group(2))
    exon_end = int(match.group(3))
    ranges = load_table4_data().get("exon_ranges", {}).get(gene, {})
    for exon_name, (start, end) in ranges.items():
        if start == exon_start and end == exon_end:
            return str(exon_name)
    return None


def parse_exon_from_duplication_notation(
    c_notation: str,
    gene: str,
) -> Optional[str]:
    """Map an exact whole-exon duplication HGVS expression to Table 4."""
    return _parse_exon_from_structural_notation(c_notation, gene, "dup")


def parse_exon_from_deletion_notation(
    c_notation: str,
    gene: str,
) -> Optional[str]:
    """Map an exact whole-exon deletion HGVS expression to Table 4."""
    return _parse_exon_from_structural_notation(c_notation, gene, "del")
