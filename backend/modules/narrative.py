# ============================================================
# ARIANE narrative summary generator
#
# Generates a human-readable explanation of a variant classification
# from the structured result produced by evaluate_variant().
# Does NOT introduce new scoring - informational only.
# ============================================================
from typing import Optional, Dict
import re

DOMAIN_DESCRIPTIONS: Dict[str, Dict[str, str]] = {
    "BRCA1": {
        "RING":        "RING domain (aa 2-101), which mediates E3 ubiquitin ligase activity through interaction with BARD1",
        "coiled_coil": "coiled-coil domain (aa 1391-1424), which mediates BRCA1 oligomerisation and interaction with PALB2",
        "BRCT":        "BRCT repeats (aa 1650-1857), which bind phosphopeptides and are essential for DNA damage signalling",
    },
    "BRCA2": {
        "PALB2_binding": "N-terminal PALB2-binding domain (aa 10-40), required for nuclear localisation and RAD51 loading",
        "DBD":           "DNA-binding domain (aa 2481-3186), containing OB folds that directly contact single-stranded DNA",
    },
}

BAYESDEL_THRESHOLDS = {
    "BRCA1": {"pp3": 0.28, "bp4": 0.15},
    "BRCA2": {"pp3": 0.30, "bp4": 0.18},
}


def _parse_intron_offset(c_notation: str) -> Optional[tuple]:
    m = re.match(r"c\.(-?\d+)([+-])(\d+)", c_notation)
    if m:
        sign = 1 if m.group(2) == "+" else -1
        return int(m.group(1)), sign * int(m.group(3))
    return None


def _extract_aa_pos(p_notation: str) -> Optional[int]:
    m = re.search(r"[A-Z][a-z]{2}(\d+)", p_notation or "")
    return int(m.group(1)) if m else None


def _three_to_one(aa3: str) -> str:
    table = {
        "Ala": "Alanine",   "Arg": "Arginine",  "Asn": "Asparagine",
        "Asp": "Aspartate", "Cys": "Cysteine",  "Gln": "Glutamine",
        "Glu": "Glutamate", "Gly": "Glycine",   "His": "Histidine",
        "Ile": "Isoleucine","Leu": "Leucine",   "Lys": "Lysine",
        "Met": "Methionine","Phe": "Phenylalanine","Pro": "Proline",
        "Ser": "Serine",    "Thr": "Threonine", "Trp": "Tryptophan",
        "Tyr": "Tyrosine",  "Val": "Valine",    "Ter": "stop",
    }
    return table.get(aa3.capitalize(), aa3)


def _opening(gene: str, c_notation: str, p_notation: str, variant_type: str) -> str:
    p = p_notation or ""
    vt = variant_type.lower()

    if vt == "frameshift":
        pos = _extract_aa_pos(p)
        tail = f" {p}" if p else ""
        if pos:
            return (f"{gene} {c_notation} introduces a frameshift at codon {pos}, "
                    f"creating a premature termination codon{tail}.")
        return f"{gene} {c_notation} is a frameshift variant{tail}."

    if vt == "nonsense":
        pos = _extract_aa_pos(p)
        tail = f" {p}" if p else ""
        return (f"{gene} {c_notation} introduces a premature termination codon at position {pos}{tail}."
                if pos else f"{gene} {c_notation} is a nonsense variant{tail}.")

    if vt == "splice_site":
        info = _parse_intron_offset(c_notation)
        if info:
            _, offset = info
            site = "donor" if offset > 0 else "acceptor"
            kind = "canonical" if abs(offset) <= 2 else "near-canonical"
            return (f"{gene} {c_notation} affects the {kind} {site} splice site "
                    f"(position {offset:+d} relative to the exon boundary).")
        return f"{gene} {c_notation} is a splice site variant."

    if vt == "missense":
        m = re.match(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?", p)
        if m:
            ref_aa = _three_to_one(m.group(1))
            pos    = m.group(2)
            alt_aa = _three_to_one(m.group(3))
            return (f"{gene} {c_notation} {p} substitutes "
                    f"{ref_aa} with {alt_aa} at position {pos}.")
        return f"{gene} {c_notation} {p} is a missense variant."

    if vt in ("synonymous", "silent"):
        pos = _extract_aa_pos(p)
        if pos:
            return (f"{gene} {c_notation} {p} is a synonymous substitution that does not "
                    f"alter the amino acid at position {pos}.")
        return f"{gene} {c_notation} is a synonymous variant."

    if vt == "intronic":
        info = _parse_intron_offset(c_notation)
        if info:
            _, offset = info
            direction = "downstream of" if offset > 0 else "upstream of"
            site = "donor" if offset > 0 else "acceptor"
            return (f"{gene} {c_notation} is an intronic variant located "
                    f"{abs(offset)} bases {direction} the {site} splice site.")
        return f"{gene} {c_notation} is an intronic variant."

    if vt == "exon_deletion":
        return f"{gene} {c_notation} is an exon-level deletion."

    if vt == "exon_duplication":
        return f"{gene} {c_notation} is an exon-level duplication."

    if vt in ("inframe_deletion",):
        return f"{gene} {c_notation} {p} is an in-frame deletion."

    if vt in ("inframe_insertion", "inframe_delins", "delins"):
        return f"{gene} {c_notation} {p} is an in-frame insertion or deletion variant."

    return f"{gene} {c_notation} {p}.".strip()


def _structural_context(gene: str, variant_type: str, criteria: dict,
                        residue_info: Optional[dict]) -> Optional[str]:
    domain = (residue_info or {}).get("domain")
    if domain:
        desc = DOMAIN_DESCRIPTIONS.get(gene, {}).get(domain, f"{domain} domain")
        known = (residue_info or {}).get("known_pathogenic_at_position", [])
        if known:
            names = ", ".join(r["aa"] for r in known[:3])
            return (f"The position falls within the {desc}. "
                    f"Known pathogenic variant(s) exist at the same position: {names}.")
        return f"The position falls within the {desc}."

    pvs1 = criteria.get("PVS1", {})
    if pvs1.get("applies"):
        reason = pvs1.get("reason", "")
        exon_m = re.search(r"in (E[\w\(\)/]+)", reason)
        exon = exon_m.group(1) if exon_m else None
        if "critical boundary" in reason or "critical C-terminal" in reason:
            if "<=" in reason or "at/before" in reason:
                # Extract the boundary value (after <=, not the PTC position)
                boundary_m = re.search(r"<= p\.(\d+)", reason)
                boundary = f" (critical boundary p.{boundary_m.group(1)})" if boundary_m else ""
                return (f"The truncation falls before the critical domain boundary{boundary}, "
                        f"removing the functionally important C-terminal region."
                        if not exon else
                        f"The truncation in exon {exon} falls before the critical domain boundary{boundary}, "
                        f"removing the functionally important C-terminal region.")
        if exon:
            return f"The variant is located in exon {exon}."

    vt = variant_type.lower()
    if vt in ("missense", "synonymous", "silent",
              "inframe_deletion", "inframe_insertion", "inframe_delins", "delins"):
        return "The position falls outside all characterised functional domains."

    return None


def _insilico_sentence(gene: str, variant_type: str, spliceai_score: Optional[float],
                       bayesdel_score: Optional[float], alphamissense: Optional[dict]) -> Optional[str]:
    parts = []
    vt = variant_type.lower()
    protein_types = {"missense", "inframe_deletion", "inframe_insertion", "inframe_delins", "delins"}
    thr = BAYESDEL_THRESHOLDS.get(gene, {"pp3": 0.28, "bp4": 0.15})

    if bayesdel_score is not None and vt in protein_types:
        if bayesdel_score >= thr["pp3"]:
            parts.append(f"BayesDel_noAF {bayesdel_score:.3f} exceeds the {gene} threshold "
                         f"({thr['pp3']}), predicting a damaging protein effect.")
        elif bayesdel_score <= thr["bp4"]:
            parts.append(f"BayesDel_noAF {bayesdel_score:.3f} is below the {gene} threshold "
                         f"({thr['bp4']}), predicting no significant protein effect.")
        else:
            parts.append(f"BayesDel_noAF {bayesdel_score:.3f} falls in the uninformative range "
                         f"({thr['bp4']}-{thr['pp3']}) for {gene}.")

    if spliceai_score is not None:
        if spliceai_score >= 0.2:
            parts.append(f"SpliceAI {spliceai_score:.3f} predicts a significant splicing effect.")
        elif spliceai_score <= 0.1:
            parts.append(f"SpliceAI {spliceai_score:.3f} predicts no significant splicing effect.")
        else:
            parts.append(f"SpliceAI {spliceai_score:.3f} is in the intermediate range.")

    if alphamissense and vt == "missense":
        am_score = alphamissense.get("am_score")
        am_class = alphamissense.get("am_class", "")
        if am_score is not None:
            labels = {
                "likely_pathogenic": "likely pathogenic",
                "ambiguous": "ambiguous",
                "likely_benign": "likely benign",
            }
            label = labels.get(am_class, am_class.replace("_", " ") if am_class else "unknown")
            parts.append(f"AlphaMissense: {am_score:.3f} ({label}).")

    return " ".join(parts) if parts else None


def _functional_evidence_sentence(criteria: dict) -> Optional[str]:
    ps3 = criteria.get("PS3", {})
    bs3 = criteria.get("BS3", {})
    if ps3.get("applies"):
        return (f"Functional assay data (Table 9) supports pathogenicity "
                f"at {ps3.get('strength', '')} level.")
    if bs3.get("applies"):
        return (f"Functional assay data (Table 9) does not support pathogenicity "
                f"({bs3.get('strength', '')} benign evidence).")
    return None


def _classification_sentence(predicted_class: int, predicted_label: str,
                              total_points: int, criteria: dict) -> str:
    applied = [
        (name, crit.get("strength", ""), crit.get("points", 0))
        for name, crit in criteria.items()
        if crit.get("applies") and crit.get("points", 0) != 0
    ]
    applied.sort(key=lambda x: -abs(x[2]))

    if not applied:
        return (f"No automatable criteria applied. "
                f"Classified as {predicted_label} (class {predicted_class}, {total_points} points).")

    # Show at most 4 most significant criteria
    crit_str = "; ".join(
        f"{name} {strength}" if strength else name
        for name, strength, _ in applied[:4]
    )
    if len(applied) > 4:
        crit_str += f" (+ {len(applied) - 4} more)"

    return (f"Classified as {predicted_label} (class {predicted_class}, {total_points} points): "
            f"{crit_str}.")


def _vus_explanation_sentence(result: dict) -> Optional[str]:
    if result.get("predicted_class") != 3:
        return None
    from backend.modules.vus_explanation import explain_vus

    explanation = explain_vus(result)
    if not explanation:
        return None
    return (
        f"VUS explanation: {explanation['summary']} "
        f"Suggested review: {explanation['what_to_check']}"
    )


def generate_narrative(
    gene: str,
    c_notation: str,
    p_notation: str,
    variant_type: str,
    result: dict,
    spliceai_score: Optional[float] = None,
    bayesdel_score: Optional[float] = None,
    alphamissense: Optional[dict] = None,
) -> str:
    """
    Generate a human-readable explanation of a variant classification.
    Purely informational - does not affect scoring.
    """
    criteria     = result.get("criteria", {})
    residue_info = result.get("residue_info")
    pred_class   = result.get("predicted_class", 3)
    pred_label   = result.get("predicted_label", "VUS")
    total_pts    = result.get("total_points", 0)

    sentences = []

    sentences.append(_opening(gene, c_notation, p_notation, variant_type))

    ctx = _structural_context(gene, variant_type, criteria, residue_info)
    if ctx:
        sentences.append(ctx)

    insilico = _insilico_sentence(gene, variant_type, spliceai_score, bayesdel_score, alphamissense)
    if insilico:
        sentences.append(insilico)

    func_ev = _functional_evidence_sentence(criteria)
    if func_ev:
        sentences.append(func_ev)

    sentences.append(_classification_sentence(pred_class, pred_label, total_pts, criteria))

    vus_sentence = _vus_explanation_sentence(result)
    if vus_sentence:
        sentences.append(vus_sentence)

    return " ".join(sentences)
