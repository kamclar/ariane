"""Manual ENIGMA evidence review and amended working classification."""

import math
from copy import deepcopy

from typing import Any, Dict, List, Mapping, Optional, Sequence

from backend.classification_dag.policy import classify_by_enigma_combination
from backend.modules.criterion_order import criterion_sort_key
from backend.modules.evidence_interactions import (
    apply_manual_rna_interactions,
    automatic_functional_interactions,
    clinical_functional_risk_interactions,
)
from backend.modules.ps1_splice_evidence import DEFINED_SOURCES as PS1_SPLICE_SOURCES
from backend.modules.bp7_rna import evaluate_bp7_rna_variant_context
from backend.modules.variant_input import normalize_variant_input
from backend.modules.variant_type import infer_variant_type
from backend.gene_policy import (
    active_genes,
    clinical_lr_thresholds,
    implementation_profile,
    policy_name,
    policy_version,
    resolve_policy_gene,
    rule_is_applicable,
    spliceai_thresholds,
    vcep_specification,
)


CSPEC_URL = ""

ENIGMA_PP4_SOURCES = {
    "15290653": "Goldgar et al. 2004",
    "12900794": "Thompson et al. 2003",
    "17924331": "Easton et al. 2007",
    "25857409": "Spurdle et al. 2015",
    "27008870": "de la Hoya et al. 2016",
    "31131967": "Parsons et al. 2019",
    "31853058": "Li et al. 2020",
    "34597585": "Caputo et al. 2021",
    "40413188": "Zanti et al. 2025",
}

STRENGTH_POINTS = {
    "Very Strong": 8,
    "Strong": 4,
    "Moderate": 2,
    "Supporting": 1,
}

MANUAL_CRITERIA = {
    "PS3": {
        "direction": "pathogenic",
        "allowed_strengths": ["Strong"],
        "title": "Calibrated functional evidence showing abnormal function",
        "threshold": "Strong only. Use after expert review confirms that the assay satisfies the ENIGMA PS3 functional-evidence specifications.",
        "check": "Confirm the assay scope, calibration against pathogenic and benign controls, the variant-specific result, the applicable strength, and whether RNA and protein effects are independent of other evidence.",
        "literature": "Use ENIGMA Specifications Table 9 as the accepted v1.2 lookup. Evidence outside Table 9 requires a documented expert calibration review under the same VCEP specifications.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PS3, Specifications Table 9 and Appendix E",
    },
    "PS4": {
        "direction": "pathogenic",
        "allowed_strengths": ["Strong"],
        "title": "Case-control enrichment",
        "threshold": "Strong when p <= 0.05, OR >= 4, the lower confidence limit excludes 2.0, and case and control datasets are matched by country and ethnicity.",
        "check": "Review the case and control definitions, ancestry matching, independence of observations, odds ratio, confidence interval, and p-value.",
        "literature": "Use peer-reviewed case-control studies and verify that the reported cohort is applicable to BRCA1/2 disease.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PS4 and Appendix F",
    },
    "PM3": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong"],
        "title": "Fanconi anemia and variants in trans",
        "threshold": "Supporting at 1 point, Moderate at 2-3 points, Strong at >= 4 points, after confirming a co-occurring P/LP variant classified using VCEP specifications and that the assessed variant does not meet benign population evidence.",
        "check": "Confirm a BRCA1/2-related Fanconi anemia phenotype, phase in trans, classification of the co-occurring variant, and per-proband scoring.",
        "literature": "Review clinical reports, segregation or phasing evidence, chromosome breakage testing, and Specifications Table 6.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PM3, Specifications Table 6 and Appendix H",
    },
    "PP1": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Quantitative co-segregation",
        "threshold": "Derived from the configured pathogenic likelihood-ratio thresholds. Very Strong also requires a predicted or experimentally proven effect on protein or mRNA splicing.",
        "check": "Use a quantitative co-segregation analysis and verify informative meioses, pedigree structure, phenotype definition, and ascertainment assumptions.",
        "literature": "Review family studies and calculate the likelihood ratio using an accepted co-segregation method.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PP1 and Appendix I",
    },
    "PP4": {
        "direction": "pathogenic",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Combined clinical likelihood ratio",
        "threshold": "Derived from the configured pathogenic likelihood-ratio thresholds.",
        "check": "Confirm that the value is a variant-specific combined clinical LR, document the included clinical data types, their independence, and the primary publication or curated source.",
        "literature": "Review ENIGMA Appendix B and Specifications Table 7. Eligible inputs may include co-segregation, co-occurrence, family history, tumour pathology, and case-control data.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PP4, Specifications Table 7 and Appendix B",
        "appendix_b_sources": [
            {"pmid": pmid, "citation": citation}
            for pmid, citation in ENIGMA_PP4_SOURCES.items()
        ],
    },
    "BS2": {
        "direction": "benign",
        "allowed_strengths": ["Supporting", "Moderate", "Strong"],
        "title": "Observation without recessive disease",
        "threshold": "Supporting at 1 point, Moderate at 2-3 points, Strong at >= 4 points, after confirming a co-occurring P/LP variant classified using VCEP specifications.",
        "check": "Confirm absence of a BRCA1/2-related Fanconi anemia phenotype and apply the per-proband stipulations.",
        "literature": "Review clinical records and Specifications Table 8; do not treat general adult non-penetrance as sufficient by itself.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, BS2, Specifications Table 8 and Appendix H",
    },
    "BS4": {
        "direction": "benign",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Quantitative lack of segregation",
        "threshold": "Derived from the configured benign likelihood-ratio thresholds.",
        "check": "Use quantitative co-segregation analysis and exclude phenocopies, pedigree errors, and incorrect phenotype assignments. To use BS4 Strong as the only Strong route to Likely Benign, document at least two independent LR components whose product equals the combined LR.",
        "literature": "Review family studies and calculate the likelihood ratio using an accepted co-segregation method.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, BS4 and Appendix I",
    },
    "BS3": {
        "direction": "benign",
        "allowed_strengths": ["Strong"],
        "title": "Calibrated functional evidence showing normal function",
        "threshold": "Strong only. Use after expert review confirms that the assay satisfies the ENIGMA BS3 functional-evidence specifications.",
        "check": "Confirm the assay scope, calibration against pathogenic and benign controls, the variant-specific result, the applicable strength, and whether RNA and protein effects are independent of other evidence.",
        "literature": "Use ENIGMA Specifications Table 9 as the accepted v1.2 lookup. Evidence outside Table 9 requires a documented expert calibration review under the same VCEP specifications.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, BS3, Specifications Table 9 and Appendix E",
    },
    "BP5": {
        "direction": "benign",
        "allowed_strengths": ["Supporting", "Moderate", "Strong", "Very Strong"],
        "title": "Combined benign clinical likelihood ratio",
        "threshold": "Derived from the configured benign likelihood-ratio thresholds.",
        "check": "Confirm that the value is a variant-specific combined clinical LR, document the included clinical data types, their independence, and the primary publication or curated source.",
        "literature": "Review ENIGMA Appendix B and Specifications Table 7. Eligible inputs may include co-segregation, co-occurrence, family history, tumour pathology, and case-control data.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, BP5, Specifications Table 7 and Appendix B",
        "appendix_b_sources": [
            {"pmid": pmid, "citation": citation}
            for pmid, citation in ENIGMA_PP4_SOURCES.items()
        ],
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
        "threshold": "Strong only, for well-established mRNA-only assays supportive of no damaging effect on transcript profile. Missense variants inside an ENIGMA functional domain must also meet BS3.",
        "check": "Confirm assay sensitivity, relevant tissue or cell type, transcript coverage, NMD sensitivity and quantification. ARIANE checks the variant type, functional-domain location and applied Table 9 BS3 evidence.",
        "literature": "Review RNA assay reports, Appendix E, and Figure 1B. Missense variants in clinically important domains must meet BS3 before BP7 Strong (RNA) can be applied.",
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
        "literature": "Review ENIGMA BRCA1/2 VCEP v1.2 PS1, Appendix J Table 17, and the documented curated reference source. ARIANE does not provide a preapproved splice-PS1 reference registry.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, PS1 splicing branch, Specifications PS1/Table 5 and Appendix J Table 17",
    },
    "PS1_PROTEIN": {
        "direction": "pathogenic",
        "allowed_strengths": ["Moderate", "Strong"],
        "title": "Same missense substitution as a VCEP-classified P/LP reference",
        "threshold": "Strong for a Pathogenic reference and Moderate for a Likely Pathogenic reference, after the complete ENIGMA protein-level PS1 reference and splice review.",
        "check": "Confirm the VCEP classification source, same normalized missense substitution, different nucleotide change, SpliceAI <= 0.1 for both variants, and no damaging splice effect in the defined reviewed sources.",
        "literature": "ST7 supplies a trusted P/LP reference candidate. Automatic protein PS1 requires a separately verified ENIGMA/ClinGen VCEP assertion or documented local VCEP reclassification, plus the complete ENIGMA splice checks. Review ENIGMA BRCA1/2 VCEP v1.2 PS1 and Appendix J.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, protein-level PS1 and Appendix J",
    },
}


def manual_criteria_for_gene(gene: str | None = None) -> Dict[str, Dict[str, Any]]:
    """Return policy-bound form definitions for one active gene."""
    policy_gene = resolve_policy_gene(gene)
    profile = implementation_profile(policy_gene)
    if profile != "enigma_brca_vcep_1_2":
        raise RuntimeError(
            f"No manual-evidence form profile is implemented for {profile!r}"
        )
    values = deepcopy(MANUAL_CRITERIA)
    specification = vcep_specification(policy_gene)
    source_prefix = f"{policy_name(policy_gene)} v{policy_version(policy_gene)}"
    for definition in values.values():
        definition["source_url"] = specification["url"]
        detail = str(definition.get("source_detail") or "")
        _prefix, separator, location = detail.partition(", ")
        definition["source_detail"] = (
            f"{source_prefix}, {location}" if separator else source_prefix
        )
    lr = clinical_lr_thresholds(policy_gene)
    pp4 = lr["pp4"]
    bp5 = lr["bp5"]
    pathogenic_text = (
        f"Supporting at LR >= {pp4['supporting_min_inclusive']:g}, "
        f"Moderate at LR >= {pp4['moderate_min_inclusive']:g}, "
        f"Strong at LR >= {pp4['strong_min_inclusive']:g}, "
        f"Very Strong at LR >= {pp4['very_strong_min_inclusive']:g}."
    )
    benign_text = (
        f"Supporting at LR <= {bp5['supporting_max_inclusive']:g}, "
        f"Moderate at LR <= {bp5['moderate_max_inclusive']:g}, "
        f"Strong at LR <= {bp5['strong_max_inclusive']:g}, "
        f"Very Strong at LR <= {bp5['very_strong_max_inclusive']:g}."
    )
    values["PP4"]["threshold"] = pathogenic_text
    values["BP5"]["threshold"] = benign_text
    values["PP1"]["threshold"] = (
        pathogenic_text[:-1]
        + " Very Strong also requires a predicted or experimentally proven effect on protein or mRNA splicing."
    )
    values["BS4"]["threshold"] = benign_text
    splice_low = spliceai_thresholds(policy_gene)["bp4"]
    values["PS1_PROTEIN"]["check"] = (
        "Confirm the VCEP classification source, same normalized missense "
        "substitution, different nucleotide change, "
        f"SpliceAI <= {splice_low} for both variants, and no damaging splice "
        "effect in the defined reviewed sources."
    )
    return {
        code: definition
        for code, definition in values.items()
        if rule_is_applicable(policy_gene, code)
    }

STRUCTURED_CURATED_CODES = {
    "PS3", "PP4", "BS3", "BP5", "PVS1_RNA", "BP7_RNA", "PVS1_INIT",
    "PS1_SPLICE", "PS1_PROTEIN",
}

_COMMON_RESOURCE_LINKS = [
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
        "url": "https://cspec.genome.network/cspec/File/id/11e62fec-23b0-4a3e-b2df-751855301746/data",
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


def resource_links_for_gene(gene: str | None = None) -> List[Dict[str, str]]:
    """Return VCEP-specific links plus shared interpretation resources."""
    genes = (resolve_policy_gene(gene),) if gene else active_genes()
    policy_links: List[Dict[str, str]] = []
    seen_urls: set[str] = set()
    for symbol in genes:
        specification = vcep_specification(symbol)
        if specification["url"] in seen_urls:
            continue
        seen_urls.add(specification["url"])
        policy_links.append({
            "title": (
                f"{symbol} {policy_name(symbol)} v{policy_version(symbol)} "
                "criteria registry"
            ),
            "url": specification["url"],
            "description": (
                f"Versioned criterion specifications and combination rules for {symbol}."
            ),
        })
    return policy_links + deepcopy(_COMMON_RESOURCE_LINKS)


def _number(evidence: Dict[str, Any], key: str) -> Optional[float]:
    value = evidence.get(key)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def evaluate_bs4_likelihood_ratio(
    evidence: Mapping[str, Any],
) -> tuple[Optional[float], bool, int]:
    """Return the BS4 LR and Table 3 single-Strong eligibility provenance."""
    aggregate = _number(dict(evidence), "likelihood_ratio")
    components = evidence.get("likelihood_ratio_components")
    if not components:
        return aggregate, False, 0
    if not isinstance(components, list):
        raise ValueError("BS4 likelihood-ratio components must be a list")

    component_lrs: list[float] = []
    independence_groups: list[str] = []
    for index, component in enumerate(components, start=1):
        if not isinstance(component, Mapping):
            raise ValueError(f"BS4 LR component {index} must be a structured record")
        try:
            lr = float(component.get("likelihood_ratio"))
        except (TypeError, ValueError):
            raise ValueError(f"BS4 LR component {index} requires a numeric LR") from None
        if lr < 0:
            raise ValueError(f"BS4 LR component {index} cannot be negative")
        source = str(component.get("source") or "").strip()
        group = str(component.get("independence_group") or "").strip()
        if not source or not group:
            raise ValueError(
                f"BS4 LR component {index} requires a source and independence group"
            )
        component_lrs.append(lr)
        independence_groups.append(group)

    if len(independence_groups) != len(set(independence_groups)):
        raise ValueError("BS4 LR components must use distinct independence groups")

    combined = math.prod(component_lrs)
    if aggregate is not None and not math.isclose(
        aggregate, combined, rel_tol=1e-6, abs_tol=1e-12
    ):
        raise ValueError(
            "BS4 reported combined LR does not equal the product of its LR components"
        )
    return combined, len(component_lrs) >= 2, len(component_lrs)


def _pp4_value_and_scale(evidence: Dict[str, Any]) -> tuple[Optional[float], str]:
    """Read the current PP4 input while retaining legacy LR audit records."""
    if evidence.get("clinical_lr_value") not in (None, ""):
        return _number(evidence, "clinical_lr_value"), str(
            evidence.get("clinical_lr_scale") or "lr"
        ).strip().lower()
    return _number(evidence, "combined_clinical_lr"), "lr"


def _pp4_strength(evidence: Dict[str, Any], gene: str) -> Optional[str]:
    value, scale = _pp4_value_and_scale(evidence)
    if value is None or value < 0 or scale not in {"lr", "log10_lr", "acmg_points"}:
        return None
    pp4 = clinical_lr_thresholds(gene)["pp4"]
    lr_thresholds = (
        pp4["supporting_min_inclusive"],
        pp4["moderate_min_inclusive"],
        pp4["strong_min_inclusive"],
        pp4["very_strong_min_inclusive"],
    )
    thresholds = {
        "lr": lr_thresholds,
        "log10_lr": tuple(math.log10(value) for value in lr_thresholds),
        "acmg_points": (1.0, 2.0, 4.0, 8.0),
    }[scale]
    supporting, moderate, strong, very_strong = thresholds
    if value >= very_strong:
        return "Very Strong"
    if value >= strong:
        return "Strong"
    if value >= moderate:
        return "Moderate"
    if value >= supporting:
        return "Supporting"
    return None


def _bp5_strength(evidence: Dict[str, Any], gene: str) -> Optional[str]:
    value, scale = _pp4_value_and_scale(evidence)
    if value is None or scale not in {"lr", "log10_lr", "acmg_points"}:
        return None
    bp5 = clinical_lr_thresholds(gene)["bp5"]
    lr_thresholds = (
        bp5["supporting_max_inclusive"],
        bp5["moderate_max_inclusive"],
        bp5["strong_max_inclusive"],
        bp5["very_strong_max_inclusive"],
    )
    thresholds = {
        "lr": lr_thresholds,
        "log10_lr": tuple(math.log10(item) for item in lr_thresholds),
        "acmg_points": (-1.0, -2.0, -4.0, -8.0),
    }[scale]
    supporting, moderate, strong, very_strong = thresholds
    if value <= very_strong:
        return "Very Strong"
    if value <= strong:
        return "Strong"
    if value <= moderate:
        return "Moderate"
    if value <= supporting:
        return "Supporting"
    return None


def _pp4_source_is_reviewed(evidence: Dict[str, Any]) -> bool:
    status = str(evidence.get("source_review_status") or "unreviewed").strip().lower()
    if status == "appendix_b":
        return str(evidence.get("source_pmid") or "").strip() in ENIGMA_PP4_SOURCES
    if status == "other_reviewed":
        return all(
            bool((evidence.get(field) or "").strip())
            for field in ("source_citation", "source_reviewed_by", "source_review_rationale")
        )
    return False


def _pp4_source_is_recorded(evidence: Dict[str, Any]) -> bool:
    status = str(evidence.get("source_review_status") or "unreviewed").strip().lower()
    if status == "appendix_b":
        return str(evidence.get("source_pmid") or "").strip() in ENIGMA_PP4_SOURCES
    if status == "other_reviewed":
        return _pp4_source_is_reviewed(evidence)
    if status == "unreviewed":
        return bool((evidence.get("source_citation") or "").strip())
    return False


def suggest_strength(
    code: str,
    evidence: Dict[str, Any],
    *,
    variant_context: Mapping[str, Any] | None = None,
    base_criteria: Sequence[Mapping[str, Any]] = (),
) -> Optional[str]:
    gene = resolve_policy_gene(
        str((variant_context or {}).get("gene") or "").strip().upper() or None
    )
    clinical_thresholds = clinical_lr_thresholds(gene)
    pp4_thresholds = clinical_thresholds["pp4"]
    bp5_thresholds = clinical_thresholds["bp5"]
    splice_low = spliceai_thresholds(gene)["bp4"]
    if code in {"PP4", "BP5"}:
        data_summary = (evidence.get("clinical_data_summary") or "").strip()
        if not data_summary or not _pp4_source_is_reviewed(evidence):
            return None
        return (
            _pp4_strength(evidence, gene)
            if code == "PP4"
            else _bp5_strength(evidence, gene)
        )

    if code in {"PS3", "BS3"}:
        expected_conclusion = "abnormal" if code == "PS3" else "normal"
        required_text_fields = (
            "assay_name",
            "source_citation",
            "calibration_summary",
            "variant_result_summary",
            "functional_reviewed_by",
        )
        if (
            evidence.get("assay_scope") not in {
                "protein_only",
                "combined_mrna_protein",
            }
            or evidence.get("functional_conclusion") != expected_conclusion
            or evidence.get("calibration_status") != "reviewed_under_enigma_vcep"
            or evidence.get("pathogenic_and_benign_controls_confirmed") is not True
            or any(
                not str(evidence.get(field) or "").strip()
                for field in required_text_fields
            )
        ):
            return None
        strength = evidence.get("curated_strength")
        return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

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
            }:
                return None
            strength = evidence.get("curated_strength")
            return strength if strength in MANUAL_CRITERIA[code]["allowed_strengths"] else None

        if evidence.get("rna_conclusion") != "no_damaging_effect":
            return None
        context_result = evaluate_bp7_rna_variant_context(
            variant_context, base_criteria
        )
        if not context_result["eligible"]:
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

    if code == "PS1_PROTEIN":
        required_text_fields = [
            "reference_variant",
            "reference_p_notation",
            "classification_source",
            "ps1_protein_rationale",
        ]
        if any(not str(evidence.get(field) or "").strip() for field in required_text_fields):
            return None
        if variant_context is not None:
            try:
                assessed = normalize_variant_input(
                    gene,
                    str(variant_context.get("c_notation") or ""),
                    p_notation=str(variant_context.get("p_notation") or ""),
                )
                reference = normalize_variant_input(
                    gene,
                    str(evidence.get("reference_variant") or ""),
                    p_notation=str(evidence.get("reference_p_notation") or ""),
                )
            except ValueError:
                return None
            if (
                infer_variant_type(assessed.c_notation, assessed.p_notation) != "missense"
                or infer_variant_type(reference.c_notation, reference.p_notation) != "missense"
                or assessed.p_notation != reference.p_notation
                or assessed.c_notation == reference.c_notation
            ):
                return None
        if evidence.get("reference_classification") not in {
            "Pathogenic", "Likely Pathogenic"
        }:
            return None
        if evidence.get("classification_verification") not in {
            "external_vcep_assertion",
            "locally_recurated_under_enigma_vcep",
        }:
            return None
        if (
            evidence.get("same_missense_confirmed") is not True
            or evidence.get("different_nucleotide_change_confirmed") is not True
            or evidence.get("splice_source_check_completed") is not True
        ):
            return None
        sources_checked = evidence.get("splice_sources_checked")
        if (
            not isinstance(sources_checked, list)
            or not set(PS1_SPLICE_SOURCES).issubset(
                {str(source).strip() for source in sources_checked}
            )
        ):
            return None
        if evidence.get("vua_confirmed_splice_status") not in {
            "none_identified", "normal"
        } or evidence.get("reference_confirmed_splice_status") not in {
            "none_identified", "normal"
        }:
            return None
        vua_score = _number(evidence, "vua_spliceai_score")
        reference_score = _number(evidence, "reference_spliceai_score")
        if (
            vua_score is None or reference_score is None
            or vua_score > splice_low or reference_score > splice_low
        ):
            return None
        if evidence.get("reference_classification_used_ps1") == "yes":
            if (
                not str(evidence.get("reference_ps1_dependency_reference") or "").strip()
                or evidence.get("direct_reciprocal_dependency_excluded") is not True
            ):
                return None
        return (
            "Strong"
            if evidence["reference_classification"] == "Pathogenic"
            else "Moderate"
        )

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
            and evidence.get("case_control_country_matched") is True
            and evidence.get("case_control_ethnicity_matched") is True
        ):
            return "Strong"
        return None

    if code in {"PM3", "BS2"}:
        if evidence.get("cooccurring_variant_classification_basis") != "vcep_specifications":
            return None
        if (
            code == "PM3"
            and evidence.get("vua_benign_population_review") != "does_not_meet"
        ):
            return None
        points = _number(evidence, "evidence_points")
        if points is None or points < 1:
            return None
        if points >= 4:
            return "Strong"
        if points >= 2:
            return "Moderate"
        return "Supporting"

    likelihood_ratio = _number(evidence, "likelihood_ratio")
    if code == "BS4":
        likelihood_ratio, _, _ = evaluate_bs4_likelihood_ratio(evidence)
    if likelihood_ratio is None or likelihood_ratio < 0:
        return None
    if code == "PP1":
        if likelihood_ratio >= pp4_thresholds["very_strong_min_inclusive"] and evidence.get("very_strong_effect_basis") in {
            "predicted_protein",
            "predicted_splicing",
            "experimental_protein",
            "experimental_splicing",
        }:
            return "Very Strong"
        if likelihood_ratio >= pp4_thresholds["strong_min_inclusive"]:
            return "Strong"
        if likelihood_ratio >= pp4_thresholds["moderate_min_inclusive"]:
            return "Moderate"
        if likelihood_ratio >= pp4_thresholds["supporting_min_inclusive"]:
            return "Supporting"
    elif code == "BS4":
        if likelihood_ratio <= bp5_thresholds["very_strong_max_inclusive"]:
            return "Very Strong"
        if likelihood_ratio <= bp5_thresholds["strong_max_inclusive"]:
            return "Strong"
        if likelihood_ratio <= bp5_thresholds["moderate_max_inclusive"]:
            return "Moderate"
        if likelihood_ratio <= bp5_thresholds["supporting_max_inclusive"]:
            return "Supporting"
    return None


def evaluate_manual_evidence(
    base_criteria: List[Dict[str, Any]],
    manual_criteria: List[Dict[str, Any]],
    variant_context: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    policy_gene = resolve_policy_gene(
        str((variant_context or {}).get("gene") or "").strip().upper() or None
    )
    combined = {
        criterion["name"]: {
            "applies": criterion.get("applies", True),
            "strength": criterion.get("strength"),
            "points": criterion.get("points", 0),
            "reason": criterion.get("reason", ""),
            "single_strong_likely_benign_eligible": criterion.get(
                "single_strong_likely_benign_eligible", False
            ),
            "single_strong_likely_benign_basis": criterion.get(
                "single_strong_likely_benign_basis", ""
            ),
            "independent_evidence_contribution_count": criterion.get(
                "independent_evidence_contribution_count", 0
            ),
        }
        for criterion in base_criteria
        if criterion.get("applies", True)
    }
    results = []

    enabled_manual = [item for item in manual_criteria if item.get("enabled")]
    bp7_context_criteria = list(base_criteria)
    for item in enabled_manual:
        if item.get("code") != "BS3":
            continue
        manual_bs3_strength = suggest_strength(
            "BS3",
            item.get("evidence", {}),
            variant_context=variant_context,
            base_criteria=base_criteria,
        )
        if manual_bs3_strength:
            bp7_context_criteria.append({
                "name": "BS3",
                "applies": True,
                "strength": manual_bs3_strength,
                "decision_path": {
                    "sources": [{
                        "source_id": "manual-enigma-vcep-functional-review",
                    }],
                },
            })
    clinical_lr_is_used = any(code in combined for code in {"PP4", "BP5"}) or any(
        item.get("code") in {"PP4", "BP5"} for item in enabled_manual
    )
    if clinical_lr_is_used:
        for item in enabled_manual:
            if item.get("code") not in {"PP1", "PS4"}:
                continue
            evidence = item.get("evidence", {})
            rationale = str(evidence.get("independence_rationale") or "").strip()
            if evidence.get("independent_from_pp4_bp5") is not True or not rationale:
                raise ValueError(
                    f"{item['code']} cannot be combined with PP4/BP5 until the reviewer "
                    "confirms independent observations and records an independence rationale"
                )

    for item in manual_criteria:
        code = item["code"]
        if item.get("enabled") and not rule_is_applicable(policy_gene, code):
            raise ValueError(
                f"{code} is not applicable under the configured VCEP policy for {policy_gene}"
            )
        definition = MANUAL_CRITERIA[code]
        if item.get("override_strength") not in {None, ""}:
            raise ValueError(
                "Manual strength overrides are not permitted. "
                f"ARIANE derives criterion strength from the configured VCEP "
                f"policy for {policy_gene}."
            )
        suggested = suggest_strength(
            code,
            item.get("evidence", {}),
            variant_context=variant_context,
            base_criteria=(
                bp7_context_criteria if code == "BP7_RNA" else base_criteria
            ),
        )
        evidence = item.get("evidence", {})
        pp4_value, pp4_scale = _pp4_value_and_scale(evidence)
        pp4_source_status = str(
            evidence.get("source_review_status") or "unreviewed"
        ).strip().lower()
        pp4_source_recorded = _pp4_source_is_recorded(evidence)
        clinical_lr_complete = (
            pp4_value is not None
            and (pp4_value >= 0 if pp4_scale == "lr" else True)
            and pp4_scale in {"lr", "log10_lr", "acmg_points"}
            and pp4_source_status in {"appendix_b", "other_reviewed", "unreviewed"}
            and pp4_source_recorded
            and bool((evidence.get("clinical_data_summary") or "").strip())
        )
        if code in {"PP4", "BP5"} and item.get("enabled") and not clinical_lr_complete:
            raise ValueError(
                f"{code} requires a clinical LR value and scale, recorded source, "
                "and clinical data summary"
            )
        if code == "BP7_RNA" and item.get("enabled") and not suggested:
            context_result = evaluate_bp7_rna_variant_context(
                variant_context, bp7_context_criteria
            )
            if not context_result["eligible"]:
                raise ValueError(context_result["reason"])
        if code in STRUCTURED_CURATED_CODES - {"PP4", "BP5"} and item.get("enabled") and not suggested:
            raise ValueError(
                f"{code} requires a complete structured curated evidence record"
            )
        selected = suggested
        applies = bool(item.get("enabled") and selected)
        points = STRENGTH_POINTS.get(selected, 0)
        if definition["direction"] == "benign":
            points *= -1
        reason = (
            f"Reviewer-curated evidence meets {code} {suggested} requirements"
            if suggested
            else f"Reviewer-curated evidence does not meet the structured {code} requirements"
        )
        single_strong_eligible = False
        single_strong_basis = ""
        contribution_count = 0
        if code == "BS4" and applies:
            _, has_multiple_lrs, contribution_count = evaluate_bs4_likelihood_ratio(
                evidence
            )
            single_strong_eligible = bool(
                selected == "Strong" and has_multiple_lrs
            )
            if single_strong_eligible:
                single_strong_basis = (
                    "Multiple independently identified segregation likelihood ratios "
                    "contribute to BS4 Strong"
                )
                reason += "; multiple independent LR components satisfy the ENIGMA Table 3 single-Strong condition"
            elif selected == "Strong":
                reason += "; BS4 Strong is valid, but the ENIGMA Table 3 single-Strong condition is not documented"
        if code == "BP5" and applies:
            evidence_types = {
                str(item).strip().lower()
                for item in evidence.get("clinical_evidence_types", [])
                if str(item).strip()
            }
            contribution_count = len(evidence_types)
            single_strong_eligible = bool(
                selected == "Strong"
                and contribution_count >= 2
                and evidence.get("independence_review_confirmed") is True
            )
            if single_strong_eligible:
                single_strong_basis = (
                    "At least two independently reviewed clinical evidence types "
                    "contribute to BP5 Strong"
                )
                reason += "; multiple clinical evidence types satisfy the ENIGMA Table 3 single-Strong condition"
            elif selected == "Strong":
                reason += "; BP5 Strong is valid, but the ENIGMA Table 3 single-Strong condition is not documented"
        results.append(
            {
                "code": code,
                "applies": applies,
                "suggested_strength": suggested,
                "selected_strength": selected,
                "points": points if applies else 0,
                "reason": reason,
                "threshold_note": definition["threshold"],
                "overridden": False,
                "notes": item.get("notes", ""),
                "references": item.get("references", []),
                "single_strong_likely_benign_eligible": single_strong_eligible,
                "single_strong_likely_benign_basis": single_strong_basis,
                "independent_evidence_contribution_count": contribution_count,
            }
        )
        if applies:
            if code == "PVS1_INIT":
                combined.pop("PP3", None)
            output_code = "PS1" if code == "PS1_PROTEIN" else code
            if output_code in combined:
                raise ValueError(
                    f"{output_code} is already present in the automated result and cannot be counted twice"
                )
            combined[output_code] = {
                "applies": True,
                "strength": selected,
                "points": points,
                "reason": reason,
                "single_strong_likely_benign_eligible": single_strong_eligible,
                "single_strong_likely_benign_basis": single_strong_basis,
                "independent_evidence_contribution_count": contribution_count,
            }

    applied_manual_codes = {
        result["code"] for result in results if result["applies"]
    }
    evidence_interactions = apply_manual_rna_interactions(
        combined, applied_manual_codes
    )
    evidence_interactions.extend(automatic_functional_interactions(combined))
    evidence_interactions.extend(clinical_functional_risk_interactions(combined))
    total_points = sum(c.get("points", 0) for c in combined.values())
    predicted_class, label, note = classify_by_enigma_combination(
        combined, total_points, gene=policy_gene
    )
    return {
        "predicted_class": predicted_class,
        "predicted_label": label,
        "total_points": total_points,
        "classification_note": note,
        "manual_criteria": sorted(
            results, key=lambda item: criterion_sort_key(item["code"])
        ),
        "evidence_interactions": evidence_interactions,
    }
