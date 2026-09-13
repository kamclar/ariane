"""Policy-bound definitions and reference links for manual evidence forms."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, List

from backend.policy.gene import (
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

ENIGMA_RECOGNISED_PP4_SOURCES = {
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
        "recognised_sources": [
            {"pmid": pmid, "citation": citation}
            for pmid, citation in ENIGMA_RECOGNISED_PP4_SOURCES.items()
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
        "recognised_sources": [
            {"pmid": pmid, "citation": citation}
            for pmid, citation in ENIGMA_RECOGNISED_PP4_SOURCES.items()
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
        "literature": "An ENIGMA ST7 v1.2 P/LP record is accepted as the reference classification basis following expert methodological review on 2026-09-07. Protein PS1 still requires the complete identity, protein-mechanism, splice and dependency checks. Review ENIGMA BRCA1/2 VCEP v1.2 PS1 and Appendix J.",
        "source_url": CSPEC_URL,
        "source_detail": "ENIGMA BRCA1/2 VCEP v1.2, protein-level PS1 and Appendix J",
    },
}

MANUAL_CRITERION_PRESENTATION = {
    "PVS1_RNA": ("rna_splicing", "RNA and splicing", 10, 10),
    "BP7_RNA": ("rna_splicing", "RNA and splicing", 10, 20),
    "PVS1_INIT": ("rna_splicing", "RNA and splicing", 10, 30),
    "PS1_SPLICE": ("rna_splicing", "RNA and splicing", 10, 40),
    "PS1_PROTEIN": ("prior_variant", "Prior variant evidence", 20, 10),
    "PS3": ("functional", "Functional evidence", 30, 10),
    "BS3": ("functional", "Functional evidence", 30, 20),
    "PS4": ("clinical", "Clinical and case-control evidence", 40, 10),
    "PP4": ("clinical", "Clinical and case-control evidence", 40, 20),
    "BP5": ("clinical", "Clinical and case-control evidence", 40, 30),
    "PP1": ("segregation", "Family and segregation evidence", 50, 10),
    "BS4": ("segregation", "Family and segregation evidence", 50, 20),
    "PM3": ("cooccurrence", "Co-occurrence evidence", 60, 10),
    "BS2": ("cooccurrence", "Co-occurrence evidence", 60, 20),
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
    for code, definition in values.items():
        group_id, group_title, group_order, criterion_order = (
            MANUAL_CRITERION_PRESENTATION[code]
        )
        definition.update({
            "group_id": group_id,
            "group_title": group_title,
            "group_order": group_order,
            "criterion_order": criterion_order,
        })
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
        "effect in the defined reviewed sources. Confirm whether the reference "
        "classification used PS1; an unknown dependency status cannot be scored."
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

