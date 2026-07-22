"""Manual ENIGMA evidence review and amended working classification."""

from typing import Any, Dict, List, Optional

from backend.modules.classifier import classify_by_enigma_combination


CSPEC_URL = "https://cspec.genome.network/cspec/ui/svi/doc/GN097"

STRENGTH_POINTS = {
    "Very Strong": 8,
    "Strong": 4,
    "Moderate": 2,
    "Supporting": 1,
}

MANUAL_CRITERIA = {
    "PS4": {
        "direction": "pathogenic",
        "allowed_strengths": ["Strong"],
        "title": "Case-control enrichment",
        "threshold": "Strong when p <= 0.05, OR >= 4, and the lower confidence limit excludes 2.0.",
        "check": "Review the case and control definitions, ancestry matching, independence of observations, odds ratio, confidence interval, and p-value.",
        "literature": "Use peer-reviewed case-control studies and verify that the reported cohort is applicable to BRCA1/2 disease.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PS4 and Appendix F",
    },
    "PM3": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong"],
        "title": "Fanconi anemia and variants in trans",
        "threshold": "Supporting at 1 point, Moderate at 2-3 points, Strong at >= 4 points.",
        "check": "Confirm a BRCA1/2-related Fanconi anemia phenotype, phase in trans, classification of the co-occurring variant, and per-proband scoring.",
        "literature": "Review clinical reports, segregation or phasing evidence, chromosome breakage testing, and Specifications Table 6.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PM3, Specifications Table 6 and Appendix H",
    },
    "PP1": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Quantitative co-segregation",
        "threshold": "Supporting at LR >= 2.08, Moderate at LR >= 4.3, Strong at LR >= 18.7, Very Strong at LR >= 350.",
        "check": "Use a quantitative co-segregation analysis and verify informative meioses, pedigree structure, phenotype definition, and ascertainment assumptions.",
        "literature": "Review family studies and calculate the likelihood ratio using an accepted co-segregation method.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PP1 and Appendix I",
    },
    "PP4": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Combined clinical likelihood ratio",
        "threshold": "Supporting at combined LR >= 2.08, Moderate at LR >= 4.3, Strong at LR >= 18.7, Very Strong at LR >= 350.",
        "check": "Confirm that the value is a variant-specific combined clinical LR, document the included clinical data types, their independence, and the primary publication or curated source.",
        "literature": "Review ENIGMA Appendix B and Specifications Table 7. Eligible inputs may include co-segregation, co-occurrence, family history, tumour pathology, and case-control data.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PP4, Specifications Table 7 and Appendix B",
    },
    "BS2": {
        "direction": "benign",
        "allowed_strengths": ["Supporting", "Moderate", "Strong"],
        "title": "Observation without recessive disease",
        "threshold": "Supporting at 1 point, Moderate at 2-3 points, Strong at >= 4 points.",
        "check": "Confirm absence of a BRCA1/2-related Fanconi anemia phenotype and apply the per-proband stipulations.",
        "literature": "Review clinical records and Specifications Table 8; do not treat general adult non-penetrance as sufficient by itself.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, BS2, Specifications Table 8 and Appendix H",
    },
    "BS4": {
        "direction": "benign",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Quantitative lack of segregation",
        "threshold": "Supporting at LR <= 0.48, Moderate at LR <= 0.23, Strong at LR <= 0.05, Very Strong at LR <= 0.00285.",
        "check": "Use quantitative co-segregation analysis and exclude phenocopies, pedigree errors, and incorrect phenotype assignments.",
        "literature": "Review family studies and calculate the likelihood ratio using an accepted co-segregation method.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, BS4 and Appendix I",
    },
    "PVS1_RNA": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "mRNA-only assay showing damaging transcript effect",
        "threshold": "Use only for well-established mRNA-only assays supportive of a damaging effect; select the ENIGMA PVS1 (RNA) strength justified by the curated RNA interpretation.",
        "check": "Confirm that the assay measures mRNA transcript profile only, documents transcript accession, tissue or cell type, NMD sensitivity, abnormal transcript products, and whether functional transcript remains.",
        "literature": "Review RNA assay reports, Appendix E, and Figure 1B. Protein-only or combined mRNA/protein assays should be evaluated under PS3/BS3 instead.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PVS1 (RNA), Figure 1B and Appendix E",
    },
    "BP7_RNA": {
        "direction": "benign",
        "allowed_strengths": ["Strong"],
        "title": "mRNA-only assay showing no damaging transcript effect",
        "threshold": "Strong only, for well-established mRNA-only assays supportive of no damaging effect on transcript profile and meeting BP7 (RNA) eligibility stipulations.",
        "check": "Confirm assay sensitivity, relevant tissue or cell type, transcript coverage, NMD sensitivity, quantification, and eligibility for BP7_Strong (RNA).",
        "literature": "Review RNA assay reports, Appendix E, and Figure 1B. Missense variants in clinically important domains require appropriate functional context before BP7_Strong (RNA).",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, BP7_Strong (RNA), Figure 1B and Appendix E",
    },
    "PVS1_INIT": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Initiation-codon PVS1 flowchart",
        "threshold": "Use only for curated start-loss variants; select the ENIGMA PVS1 initiation-codon strength justified by Specifications Table 4 and Appendix D.",
        "check": "Confirm Met1/start-loss, whether an in-frame alternative start codon is available, evidence for pathogenic variants upstream of the nearest alternative start, and the expected N-terminal functional impact.",
        "literature": "Review the ENIGMA PVS1 initiation-codon flowchart, Appendix D, gene-specific transcript context, and supporting pathogenic variant evidence.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PVS1 initiation codon flowchart, Specifications Table 4 and Appendix D",
    },
    "PS1_SPLICE": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong"],
        "title": "Same splicing impact as known P/LP variant",
        "threshold": "Use only after curated PS1(splicing) review: the VUA must have the same predicted/proven splice event as a known P/LP reference variant, with similar or stronger prediction evidence; select the ENIGMA Appendix J/Table 17 strength manually.",
        "check": "Confirm the reference variant, its P/LP classification source, the exact shared splice event, prediction strength comparison, and Appendix J/Table 17 weight. For exonic variants, consider any predicted or proven protein/missense effect before applying PS1(splicing).",
        "literature": "Review ENIGMA BRCA1/2 VCEP v1.2 PS1, Appendix J Table 17, and the curated reference source. The pilot splice_ps1_reference_set.json is a seed only and requires expert review before use.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PS1 splicing branch, Specifications PS1/Table 5 and Appendix J Table 17",
    },
}

STRUCTURED_CURATED_CODES = {"PP4", "PVS1_RNA", "BP7_RNA", "PVS1_INIT", "PS1_SPLICE"}

RESOURCE_LINKS = [
    {
        "title": "ENIGMA BRCA1/2 VCEP v1.2 criteria registry",
        "url": CSPEC_URL,
        "description": "Versioned BRCA1/2 criterion specifications and combination rules used by ARIANE.",
    },
    {
        "title": "ACMG/AMP sequence variant interpretation guidelines",
        "url": "https://pubmed.ncbi.nlm.nih.gov/25741868/",
        "description": "Foundational 2015 ACMG/AMP framework for sequence variant interpretation.",
    },
    {
        "title": "Tavtigian et al. point-based classification framework",
        "url": "https://pubmed.ncbi.nlm.nih.gov/32720330/",
        "description": "Naturally scaled point system used for contradictory evidence in ARIANE.",
    },
    {
        "title": "Specifications v1.2",
        "url": "https://cspec.genome.network/cspec/File/id/02537f62-66a3-4e67-8aec-cf44b326534d/data",
        "description": "Full BRCA1/2 criterion specifications, flowcharts, and supporting tables.",
    },
    {
        "title": "Appendix v1.2",
        "url": "https://cspec.genome.network/cspec/File/id/5a75d1a0-1222-46a2-8802-68a4f2251a3a/data",
        "description": "Detailed calibration evidence for the criteria.",
    },
    {
        "title": "Supplementary tables v1.2",
        "url": "https://cspec.genome.network/cspec/File/id/3dadda2f-94a3-497f-aa35-3bb6e828ddd5/data",
        "description": "Supplementary tables including evidence calibration material.",
    },
    {
        "title": "Specifications Table 4",
        "url": "https://cspec.genome.network/cspec/File/id/10301df8-45e0-4309-adba-c121eb057d3e/data",
        "description": "PVS1 and PM5 exon-level lookup used by ARIANE.",
    },
    {
        "title": "Specifications Table 9",
        "url": "https://cspec.genome.network/cspec/File/id/c540f11d-0be2-45d6-a0bf-ae5327a04885/data",
        "description": "Calibrated PS3 and BS3 functional evidence used by ARIANE.",
    },
    {
        "title": "ClinVar review status",
        "url": "https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/",
        "description": "Official explanation of ClinVar review stars.",
    },
]


def _number(evidence: Dict[str, Any], key: str) -> Optional[float]:
    value = evidence.get(key)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def suggest_strength(code: str, evidence: Dict[str, Any]) -> Optional[str]:
    if code == "PP4":
        likelihood_ratio = _number(evidence, "combined_clinical_lr")
        source = (evidence.get("source_citation") or "").strip()
        data_summary = (evidence.get("clinical_data_summary") or "").strip()
        if likelihood_ratio is None or likelihood_ratio < 0 or not source or not data_summary:
            return None
        if likelihood_ratio >= 350:
            return "Very Strong"
        if likelihood_ratio >= 18.7:
            return "Strong"
        if likelihood_ratio >= 4.3:
            return "Moderate"
        if likelihood_ratio >= 2.08:
            return "Supporting"
        return None

    if code in {"PVS1_RNA", "BP7_RNA"}:
        assay_scope = evidence.get("assay_scope")
        transcript_accession = (evidence.get("transcript_accession") or "").strip()
        tissue = (evidence.get("tissue_or_cell_type") or "").strip()
        nmd = evidence.get("nmd_assessed")
        if (
            assay_scope != "mrna_only"
            or not transcript_accession
            or not tissue
            or nmd not in {"yes", "no", "not_applicable"}
        ):
            return None

        if code == "PVS1_RNA":
            if evidence.get("rna_conclusion") != "damaging":
                return None
            if evidence.get("functional_transcript_remaining") not in {
                "absent_or_minimal",
                "reduced",
                "uncertain",
            }:
                return None
            strength = evidence.get("curated_strength")
            return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

        if evidence.get("rna_conclusion") != "no_damaging_effect":
            return None
        if evidence.get("bp7_rna_eligible") is not True:
            return None
        return "Strong"

    if code == "PVS1_INIT":
        if evidence.get("met1_loss_confirmed") is not True:
            return None
        if evidence.get("alternative_start_assessed") not in {"yes", "no"}:
            return None
        if evidence.get("upstream_pathogenic_evidence") not in {
            "yes",
            "no",
            "not_applicable",
        }:
            return None
        if evidence.get("functional_domain_impact") not in {
            "yes",
            "no",
            "uncertain",
        }:
            return None
        nearest_start = (evidence.get("nearest_alternative_start") or "").strip()
        rationale = (evidence.get("initiation_flowchart_rationale") or "").strip()
        if not nearest_start or not rationale:
            return None
        strength = evidence.get("curated_strength")
        return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

    if code == "PS1_SPLICE":
        required_text_fields = [
            "reference_variant",
            "reference_classification_source",
            "vua_splice_event",
            "reference_splice_event",
            "ps1_splice_rationale",
        ]
        if any(not (evidence.get(field) or "").strip() for field in required_text_fields):
            return None
        if evidence.get("reference_classification") not in {
            "Pathogenic",
            "Likely Pathogenic",
        }:
            return None
        if evidence.get("same_splice_event_confirmed") is not True:
            return None
        if evidence.get("prediction_strength_comparison") not in {
            "similar",
            "stronger",
        }:
            return None
        strength = evidence.get("curated_strength")
        return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

    if code == "PS4":
        p_value = _number(evidence, "p_value")
        odds_ratio = _number(evidence, "odds_ratio")
        lower_ci = _number(evidence, "lower_ci")
        if (
            p_value is not None
            and odds_ratio is not None
            and lower_ci is not None
            and p_value <= 0.05
            and odds_ratio >= 4
            and lower_ci > 2
        ):
            return "Strong"
        return None

    if code in {"PM3", "BS2"}:
        points = _number(evidence, "evidence_points")
        if points is None or points < 1:
            return None
        if points >= 4:
            return "Strong"
        if points >= 2:
            return "Moderate"
        return "Supporting"

    likelihood_ratio = _number(evidence, "likelihood_ratio")
    if likelihood_ratio is None or likelihood_ratio < 0:
        return None
    if code == "PP1":
        if likelihood_ratio >= 350:
            return "Very Strong"
        if likelihood_ratio >= 18.7:
            return "Strong"
        if likelihood_ratio >= 4.3:
            return "Moderate"
        if likelihood_ratio >= 2.08:
            return "Supporting"
    elif code == "BS4":
        if likelihood_ratio <= 0.00285:
            return "Very Strong"
        if likelihood_ratio <= 0.05:
            return "Strong"
        if likelihood_ratio <= 0.23:
            return "Moderate"
        if likelihood_ratio <= 0.48:
            return "Supporting"
    return None


def evaluate_manual_evidence(
    base_criteria: List[Dict[str, Any]],
    manual_criteria: List[Dict[str, Any]],
) -> Dict[str, Any]:
    combined = {
        criterion["name"]: {
            "applies": criterion.get("applies", True),
            "strength": criterion.get("strength"),
            "points": criterion.get("points", 0),
            "reason": criterion.get("reason", ""),
        }
        for criterion in base_criteria
        if criterion.get("applies", True)
    }
    results = []

    for item in manual_criteria:
        code = item["code"]
        definition = MANUAL_CRITERIA[code]
        suggested = suggest_strength(code, item.get("evidence", {}))
        override = item.get("override_strength") or None
        evidence = item.get("evidence", {})
        pp4_complete = (
            _number(evidence, "combined_clinical_lr") is not None
            and _number(evidence, "combined_clinical_lr") >= 0
            and bool((evidence.get("source_citation") or "").strip())
            and bool((evidence.get("clinical_data_summary") or "").strip())
        )
        if code == "PP4" and item.get("enabled") and not pp4_complete:
            raise ValueError("PP4 requires combined clinical LR, source citation, and clinical data summary")
        if code in STRUCTURED_CURATED_CODES - {"PP4"} and item.get("enabled") and not suggested:
            raise ValueError(
                f"{code} requires a complete structured curated evidence record"
            )
        if code in STRUCTURED_CURATED_CODES:
            override = None
        if override and override not in definition["allowed_strengths"]:
            raise ValueError(f"{override} is not an ENIGMA v1.2 strength for {code}")
        selected = override or suggested
        applies = bool(item.get("enabled") and selected)
        points = STRENGTH_POINTS.get(selected, 0)
        if definition["direction"] == "benign":
            points *= -1
        reason = (
            f"User-provided evidence; ARIANE suggestion: {suggested or 'threshold not met'}"
        )
        if override:
            reason += f"; reviewer selected: {override}"
        results.append(
            {
                "code": code,
                "applies": applies,
                "suggested_strength": suggested,
                "selected_strength": selected,
                "points": points if applies else 0,
                "reason": reason,
                "threshold_note": definition["threshold"],
                "overridden": bool(override and override != suggested),
                "notes": item.get("notes", ""),
                "references": item.get("references", []),
            }
        )
        if applies:
            if code in {"PVS1_RNA", "PVS1_INIT"}:
                combined.pop("PP3", None)
            combined[code] = {
                "applies": True,
                "strength": selected,
                "points": points,
                "reason": reason,
            }

    total_points = sum(c.get("points", 0) for c in combined.values())
    predicted_class, label, note = classify_by_enigma_combination(
        combined, total_points
    )
    return {
        "predicted_class": predicted_class,
        "predicted_label": label,
        "total_points": total_points,
        "classification_note": note,
        "manual_criteria": results,
    }
