import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "notebooks" / "BRCA_ACMG_Criteria_Module1_v1_7_0.ipynb"
TARGET = ROOT / "notebooks" / "BRCA_ACMG_Criteria_Module1_v1_8_0.ipynb"


def read(path):
    return path.read_text(encoding="utf-8")


def source_lines(text):
    return text.splitlines(keepends=True)


def extract(text, start, end=None):
    start_index = text.index(start)
    if end is None:
        return text[start_index:]
    return text[start_index:text.index(end, start_index)]


def code_cell(text):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source_lines(text),
    }


def markdown_cell(text):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": source_lines(text),
    }


def replace_cell(cells, marker, replacement):
    for index, cell in enumerate(cells):
        if marker in "".join(cell.get("source", [])):
            cells[index] = replacement
            return index
    raise ValueError(f"Cell marker not found: {marker}")


def find_cell(cells, marker):
    for index, cell in enumerate(cells):
        if marker in "".join(cell.get("source", [])):
            return index
    raise ValueError(f"Cell marker not found: {marker}")


def strip_module_imports(text, first_symbol, end_marker=None):
    body = extract(text, first_symbol, end_marker)
    return body.rstrip() + "\n"


def build_parity_helpers():
    hgvs = read(ROOT / "backend" / "modules" / "hgvs.py")
    variant_type = read(ROOT / "backend" / "modules" / "variant_type.py")
    ps1 = strip_module_imports(
        read(ROOT / "backend" / "modules" / "ps1.py"),
        "ST7_PATH =",
    )
    pp4_bp5 = strip_module_imports(
        read(ROOT / "backend" / "modules" / "pp4_bp5.py"),
        "PP4_POINTS =",
    )
    classifier = strip_module_imports(
        read(ROOT / "backend" / "modules" / "classifier.py"),
        "def classify_by_points",
        "def verify_acmg_combination",
    )
    return f"""# ============================================================
# v1.8.0 application parity helpers
# HGVS normalization, normalized variant type, ST7 PS1, ST7 PP4/BP5,
# and ENIGMA VCEP v1.2 Table 3 combination classification.
# ============================================================

import json
import re
from pathlib import Path
from typing import Optional, Dict, List

{hgvs}

{variant_type}

# ST7 reference file used by PS1 and PP4/BP5.
ST7_PATH = PROJECT_ROOT / "data" / "st7_reference_set.json"
if not ST7_PATH.exists():
    ST7_PATH = Path("st7_reference_set.json")

{ps1}

{pp4_bp5}

{classifier}
"""


def build_evaluate_variant():
    return """def evaluate_variant(variant: Dict, use_gnomad_cache: bool = True) -> Dict:
    # Evaluate the same automatable Module 1 criteria as the ARIANE application.
    gene = variant["gene"].strip().upper()
    c_notation, p_notation = split_combined_hgvs(
        variant["c_notation"], variant.get("p_notation", "")
    )
    variant_type = infer_variant_type(c_notation, p_notation)

    results = {
        "variant": f"{gene} {c_notation} {p_notation}",
        "expert_class": variant.get("expert_class"),
        "criteria": {},
        "total_points": 0,
        "warnings": [],
        "has_functional_evidence": False,
        "classification_note": "",
    }
    results["warnings"].append(
        "FIRST PASS - automatable ENIGMA VCEP v1.2 rules only. "
        "Manual review is still required for PS4, PM3, PP1, BS2 and BS4."
    )

    # Resolve frequency data first because BA1 is stand-alone benign.
    if use_gnomad_cache:
        grch37 = get_grch37(RESOLVED_VARIANTS, gene, c_notation)
        grch38 = get_grch38(RESOLVED_VARIANTS, gene, c_notation)
        gnomad_data = get_gnomad_frequencies(grch37=grch37, grch38=grch38)
    else:
        gnomad_data = {
            "status": "not_queried",
            "found": None,
            "datasets": {},
            "coverage": {"status": "not_evaluated", "passes_pm2": False, "datasets": {}},
            "max_af": None,
            "frequency_metric": None,
            "pm2_absence_established": False,
            "pm2_coverage_ok": False,
            "source": None,
            "cache_mode": "not_used",
        }
    results["gnomad"] = gnomad_data
    results["gnomad_summary"] = gnomad_status_summary(gnomad_data)

    freq_criteria = evaluate_frequency_criteria(gnomad_data, variant_type)
    for crit_name, crit_data in freq_criteria.items():
        if crit_data.get("applies"):
            results["criteria"][crit_name] = crit_data
            results["total_points"] += crit_data["points"]
        elif crit_name == "PM2":
            results["warnings"].append(crit_data["reason"])
    if use_gnomad_cache:
        results["warnings"].append("gnomAD: " + results["gnomad_summary"])

    if "BA1" in results["criteria"]:
        results["predicted_class"] = 1
        results["classification_note"] = "BA1 stand-alone benign"
        return results

    spliceai_score = get_spliceai_score(gene, c_notation)
    bayesdel_score = get_bayesdel_score(gene, c_notation)

    # Table 9 functional evidence.
    table9 = table9_lookup_ps3_bs3(gene, c_notation)
    if table9["applies"]:
        results["criteria"][table9["code"]] = {
            "applies": True,
            "strength": table9["strength"],
            "points": table9["points"],
            "reason": table9["reason"],
        }
        results["total_points"] += table9["points"]
        results["has_functional_evidence"] = True

    # Table 4 PVS1 and PM5.
    pvs1 = evaluate_pvs1(
        gene, variant_type, p_notation,
        c_notation=c_notation,
        spliceai_score=spliceai_score,
        dup_type=variant.get("dup_type", "Unknown"),
    )
    if pvs1["applies"]:
        results["criteria"]["PVS1"] = pvs1
        results["total_points"] += pvs1["points"]
    elif pvs1.get("requires_rna") or variant_type in [
        "initiation_codon", "exon_deletion", "exon_duplication"
    ]:
        results["warnings"].append(pvs1["reason"])

    if pvs1.get("pm5_code") and pvs1.get("pm5_strength"):
        results["criteria"]["PM5_PTC"] = {
            "applies": True,
            "strength": pvs1["pm5_strength"],
            "points": pvs1["pm5_points"],
            "reason": f"Table 4: {pvs1['pm5_code']} for PTC in {pvs1.get('exon', 'unknown exon')}",
        }
        results["total_points"] += pvs1["pm5_points"]

    # ST7 multifactorial likelihood and same-amino-acid reference lookup.
    pp4_bp5 = evaluate_pp4_bp5(gene, c_notation)
    if pp4_bp5.get("applies"):
        results["criteria"][pp4_bp5["code"]] = pp4_bp5
        results["total_points"] += pp4_bp5["points"]

    ps1 = evaluate_ps1(
        gene, c_notation, p_notation,
        variant_type=variant_type,
        spliceai_score=spliceai_score,
    )
    if ps1.get("applies"):
        results["criteria"]["PS1"] = ps1
        results["total_points"] += ps1["points"]

    # Computational evidence. Table 9 suppresses benign computational evidence.
    suppress_benign_comp = results["has_functional_evidence"]
    pp3_bp4 = evaluate_pp3_bp4(
        gene, variant_type, p_notation,
        bayesdel_score=bayesdel_score,
        spliceai_score=spliceai_score,
    )
    for crit_name, crit_data in pp3_bp4.items():
        if not crit_data.get("applies"):
            continue
        if crit_name == "PP3" or not suppress_benign_comp:
            results["criteria"][crit_name] = crit_data
            results["total_points"] += crit_data["points"]
        else:
            results["warnings"].append(
                f"{crit_name} not applied - functional evidence in Table 9 "
                "overrides benign computational prediction."
            )

    # BP7 is added to BP4 for eligible synonymous and intronic variants.
    if variant_type in ["synonymous", "silent", "intronic"] and not suppress_benign_comp:
        aa_pos = get_amino_acid_position(p_notation)
        in_domain = is_in_functional_domain(gene, aa_pos)[0] if aa_pos else False
        bp4_met = "BP4" in results["criteria"] and results["criteria"]["BP4"].get("applies", False)
        bp7 = evaluate_bp7(
            variant_type,
            spliceai_score=spliceai_score,
            in_domain=in_domain,
            bp4_met=bp4_met,
            c_notation=c_notation,
        )
        if bp7["applies"]:
            results["criteria"]["BP7"] = bp7
            results["total_points"] += bp7["points"]

    bp1 = evaluate_bp1(gene, variant_type, p_notation, spliceai_score=spliceai_score)
    if bp1["applies"]:
        results["criteria"]["BP1"] = bp1
        results["total_points"] += bp1["points"]

    if spliceai_score is None:
        results["warnings"].append(
            f"SpliceAI not available for {gene} {c_notation} - "
            "benign criteria BP1/BP4/BP7 and PS1 require confirmed low score"
        )
    if bayesdel_score is None:
        results["warnings"].append(f"BayesDel_noAF not available for {gene} {c_notation}")

    cls, label, note = classify_by_enigma_combination(
        results["criteria"], results["total_points"]
    )
    results["predicted_class"] = cls
    results["predicted_label"] = label
    results["classification_note"] = note or (
        "First pass - automatable ENIGMA VCEP v1.2 rules only. "
        "Non-automated criteria may affect final classification."
    )
    return results


# ============================================================
# Run evaluation on all test variants
# ============================================================

print("Evaluating all test variants...")
print("=" * 80)
print()

all_results = []
USE_GNOMAD_CACHE = True

for i, variant in enumerate(test_variants):
    print(f"[{i+1}/{len(test_variants)}] {variant['gene']} {variant['c_notation']}")
    print("-" * 40)
    result = evaluate_variant(variant, use_gnomad_cache=USE_GNOMAD_CACHE)
    all_results.append(result)
    print(f"  Expert class:    {result['expert_class']}")
    print(f"  Predicted class: {result['predicted_class']}")
    print(f"  Total points:    {result['total_points']}")
    print(f"  Criteria:        {list(result['criteria'].keys())}")
    for warning in result["warnings"]:
        print(f"  Warning:         {warning}")
    print()

print("=" * 80)
print("Evaluation complete!")
"""


SUMMARY = """## Summary - Version 1.8.0

### What's New in v1.8.0

This notebook now mirrors the active ARIANE application Module 1 scoring logic.

1. Added normalized HGVS input handling, including combined input such as
   `c.6147_6149del (p.Val2050del)`.
2. Added normalized variant type inference, including initiation codon,
   intronic BP7, 5' UTR and 3' UTR guards.
3. Added ST7 `PS1` lookup for a different nucleotide change producing the same
   missense amino-acid change as a known P/LP variant. `PS1` requires confirmed
   `SpliceAI <= 0.1`.
4. Added ST7 multifactorial likelihood criteria `PP4` and `BP5`.
5. Added ENIGMA VCEP v1.2 Table 3 combination classification as the default.
   Tavtigian 2020 points are used only when benign and pathogenic evidence
   conflict.
6. Corrected BRCA1 RING domain to aa 2-101.
7. Corrected PVS1 handling: RNA-dependent PVS1 is not automatically scored,
   non-canonical splice prediction does not create PVS1, and initiation codon
   variants are recognized but flagged for manual review.
8. Added intronic `BP7_Supporting` outside conserved donor and acceptor motifs.
9. Benign computational `BP4` and `BP7` are suppressed when curated Table 9
   functional evidence is present.

### Implemented Module 1 Criteria

| Criterion | Status | Source |
|-----------|--------|--------|
| BA1 | Complete | local gnomAD cache, depth >= 20 |
| BS1_Strong / BS1_Supporting | Complete | local gnomAD cache, depth >= 20 |
| PM2_Supporting | Complete | absent in gnomAD v2.1.1 and v3.1.2, depth >= 25, non-indels |
| PVS1 and PM5_PTC | Table 4 | PTC, curated splice variants, exon deletion and duplication |
| PS3 / BS3 | Table 9 | curated functional evidence |
| PS1 | ST7 | same missense amino-acid change, different nucleotide change |
| PP4 / BP5 | ST7 | multifactorial likelihood |
| PP3 / BP4 | Complete | SpliceAI and gene-specific BayesDel thresholds |
| BP1_Strong | Complete | outside functional domain with confirmed low SpliceAI |
| BP7_Supporting | Complete | synonymous and intronic branches, requires BP4 |

### Intentionally Not Automatically Scored

1. `PS4`, `PM3`, `PP1`, `BS2`, `BS4`: require structured clinical or family data.
2. RNA-dependent `PVS1 (RNA)` and `BP7_Strong (RNA)`: require curated RNA input.
3. Initiation codon flowchart: recognized and flagged, pending curated data rule.
4. Splice `PS1`: the implemented `PS1` branch covers protein missense changes.
5. Exon CNV scoring requires a parseable HGVS notation. For duplications,
   `dup_type="Tandem"` must be supplied only when supported by laboratory data.

### Data Files Required

Place these files in `PROJECT_ROOT/data/`:

- `enigma_table4.json`
- `enigma_table9.json`
- `st7_reference_set.json`
- local gnomAD regional cache under `data/gnomad/`

### External Read-Only Comparison

ClinVar and ClinGen Evidence Repository lookups remain informational. They do
not add automatic ACMG/ENIGMA points.

### Important Classification Note

The notebook reports a first-pass Module 1 result. It does not replace expert
review or a complete clinical classification.
"""


def main():
    notebook = json.loads(SOURCE.read_text(encoding="utf-8"))
    cells = notebook["cells"]

    # Clear stale v1.7.0 outputs before writing the new notebook.
    for cell in cells:
        if cell.get("cell_type") == "code":
            cell["execution_count"] = None
            cell["outputs"] = []

    # Intro and corrected domain boundary.
    intro = "".join(cells[0]["source"])
    intro = intro.replace("v1.7.0", "v1.8.0")
    intro = re.sub(
        r"Version 1\.7\.0 adds.*?\n\n",
        "Version 1.8.0 aligns the standalone notebook with the active ARIANE "
        "application Module 1 rules, including ST7 evidence and ENIGMA Table 3 "
        "combination classification.\\n\\n",
        intro,
        flags=re.DOTALL,
    )
    cells[0]["source"] = source_lines(intro)

    domain_index = find_cell(cells, "FUNCTIONAL_DOMAINS =")
    domain_source = "".join(cells[domain_index]["source"])
    cells[domain_index]["source"] = source_lines(domain_source.replace('"RING": (1, 101)', '"RING": (2, 101)'))

    # Replace rule implementations with their active backend counterparts.
    pvs1 = strip_module_imports(
        read(ROOT / "backend" / "modules" / "pvs1.py"),
        "def evaluate_pvs1",
        'if __name__ == "__main__":',
    )
    pp3_bp4_module = read(ROOT / "backend" / "modules" / "pp3_bp4.py")
    pp3_bp4 = extract(pp3_bp4_module, "BAYESDEL_THRESHOLDS =", 'if __name__ == "__main__":')
    bp7 = strip_module_imports(
        read(ROOT / "backend" / "modules" / "bp7.py"),
        "def evaluate_bp7",
        'if __name__ == "__main__":',
    )
    replace_cell(cells, "def evaluate_pvs1(", code_cell(pvs1))
    replace_cell(cells, "BAYESDEL_THRESHOLDS =", code_cell(pp3_bp4))
    replace_cell(cells, "def evaluate_bp7(", code_cell(bp7))

    # Add new parity helpers immediately before the main orchestrator.
    orchestrator_index = find_cell(cells, "def evaluate_variant(variant:")
    cells.insert(orchestrator_index, markdown_cell("## Application parity helpers added in v1.8.0\n"))
    cells.insert(orchestrator_index + 1, code_cell(build_parity_helpers()))
    replace_cell(cells, "def evaluate_variant(variant:", code_cell(build_evaluate_variant()))

    replace_cell(cells, "## Summary - Version 1.7.0", markdown_cell(SUMMARY))

    TARGET.write_text(
        json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {TARGET}")


if __name__ == "__main__":
    main()
