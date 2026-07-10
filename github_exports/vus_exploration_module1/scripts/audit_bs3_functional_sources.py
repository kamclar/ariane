"""Audit BRCA1 BS3 evidence against downloaded Findlay 2018 SGE source data."""

from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
SOURCE_BASE_DIR = ANALYSIS_DIR / "external_sources" / "brca1_bs3_functional_evidence_audit_sources"
SOURCE_DIR = SOURCE_BASE_DIR / "findlay_2018_brca1_sge"
RAW_FINDLAY_CSV = SOURCE_DIR / "Findlay_2018_Nature_Supplementary_Table_1_SGE_scores.csv"
NORMALIZED_FINDLAY_CSV = SOURCE_DIR / "Findlay_2018_Nature_Supplementary_Table_1_SGE_scores.normalized.csv"
ENIGMA_TABLE9 = ROOT / "backend" / "data" / "enigma_table9.json"
TABLE_DIR = ANALYSIS_DIR / "tables" / "bs3_functional_source_audit"
REPORT = ANALYSIS_DIR / "bs3_functional_source_audit_report.md"


def read_findlay_normalized() -> list[dict[str, str]]:
    with RAW_FINDLAY_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    header = rows[2]
    output = []
    for row in rows[3:]:
        if not row or not any(row):
            continue
        padded = row + [""] * (len(header) - len(row))
        output.append({header[i]: padded[i] for i in range(len(header))})
    with NORMALIZED_FINDLAY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(output)
    return output


def load_enigma_brca1_bs3() -> list[dict[str, str]]:
    data = json.load(ENIGMA_TABLE9.open(encoding="utf-8"))["variants"]
    return [
        dict(value)
        for value in data.values()
        if value.get("gene") == "BRCA1" and value.get("code") == "BS3"
    ]


def cds_pos(c_notation: str) -> int | None:
    if not c_notation.startswith("c."):
        return None
    rest = c_notation[2:]
    digits = ""
    for ch in rest:
        if ch.isdigit():
            digits += ch
        elif digits:
            break
    return int(digits) if digits else None


def region(pos: int | None) -> str:
    if pos is None:
        return "unknown"
    aa = (pos + 2) // 3
    if 1 <= aa <= 109:
        return "BRCA1_N_terminal_RING_region"
    if 1642 <= aa <= 1855:
        return "BRCA1_BRCT_region"
    if 24 <= aa <= 65:
        return "BRCA1_RING_zinc_finger_core"
    return "other_BRCA1_region"


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    if not rows:
        return "| none |\n"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def main() -> None:
    findlay = read_findlay_normalized()
    table9_bs3 = load_enigma_brca1_bs3()
    findlay_by_c = {row["transcript_variant"]: row for row in findlay}

    audit_rows = []
    for row in table9_bs3:
        c = row["c_notation"]
        f = findlay_by_c.get(c)
        pos = cds_pos(c)
        audit_rows.append(
            {
                "gene": "BRCA1",
                "c_notation": c,
                "p_notation": row.get("p_notation", ""),
                "region": region(pos),
                "table9_strength": row.get("strength", ""),
                "table9_text": row.get("text", ""),
                "in_findlay_sge": "yes" if f else "no",
                "findlay_func_class": "" if not f else f.get("func.class", ""),
                "findlay_function_score_mean": "" if not f else f.get("function.score.mean", ""),
                "findlay_mean_rna_score": "" if not f else f.get("mean.rna.score", ""),
                "findlay_consequence": "" if not f else f.get("consequence", ""),
                "findlay_experiment": "" if not f else f.get("experiment", ""),
            }
        )

    region_counts = Counter(row["region"] for row in audit_rows)
    region_in_findlay = Counter(row["region"] for row in audit_rows if row["in_findlay_sge"] == "yes")
    func_counts = Counter(
        (row["region"], row["findlay_func_class"] or "not_in_findlay")
        for row in audit_rows
    )

    summary_rows = []
    for key, total in region_counts.most_common():
        matched = region_in_findlay[key]
        summary_rows.append(
            {
                "region": key,
                "brca1_bs3_table9_count": total,
                "matched_findlay_sge_count": matched,
                "matched_percent": f"{matched / total * 100:.2f}" if total else "0.00",
            }
        )

    func_rows = [
        {"region": key[0], "findlay_func_class": key[1], "count": value}
        for key, value in func_counts.most_common()
    ]

    write_csv(
        TABLE_DIR / "brca1_bs3_table9_vs_findlay_sge.csv",
        audit_rows,
        [
            "gene",
            "c_notation",
            "p_notation",
            "region",
            "table9_strength",
            "table9_text",
            "in_findlay_sge",
            "findlay_func_class",
            "findlay_function_score_mean",
            "findlay_mean_rna_score",
            "findlay_consequence",
            "findlay_experiment",
        ],
    )
    write_csv(
        TABLE_DIR / "brca1_bs3_findlay_match_summary.csv",
        summary_rows,
        ["region", "brca1_bs3_table9_count", "matched_findlay_sge_count", "matched_percent"],
    )
    write_csv(
        TABLE_DIR / "brca1_bs3_findlay_function_class_summary.csv",
        func_rows,
        ["region", "findlay_func_class", "count"],
    )

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# BS3 Functional Source Audit

Generated: {generated}

## Sources

- Local ENIGMA Table 9 JSON: `{ENIGMA_TABLE9.relative_to(ROOT)}`
- Downloaded source: Findlay et al. 2018 Nature Supplementary Table 1,
  `external_sources/brca1_bs3_functional_evidence_audit_sources/findlay_2018_brca1_sge/Findlay_2018_Nature_Supplementary_Table_1_SGE_scores.xlsx`
- Downloaded source guide: Findlay et al. 2018 Nature Supplementary Information,
  `external_sources/brca1_bs3_functional_evidence_audit_sources/findlay_2018_brca1_sge/Findlay_2018_Nature_Supplementary_Information.pdf`
- Downloaded update source: Dace & Findlay 2023 interim report cited by ENIGMA Table 9,
  `external_sources/brca1_bs3_functional_evidence_audit_sources/dace_findlay_2023_interim/Dace_Findlay_2023_BRCA1_interim_report.pdf`

## Purpose

This audit checks whether BRCA1 `BS3` entries in the local ENIGMA Table 9
resource are supported variant-by-variant by Findlay et al. 2018 saturation
genome editing source rows, rather than being inferred only from broad domain
membership.

## Match Summary

{markdown_table(summary_rows, ["region", "brca1_bs3_table9_count", "matched_findlay_sge_count", "matched_percent"])}

## Function Class Summary

{markdown_table(func_rows[:30], ["region", "findlay_func_class", "count"])}

## Preliminary Interpretation

Most BRCA1 `BS3` records in Table 9 appear to cite PMID:30209399, which is the
Findlay et al. BRCA1 saturation genome editing study. The audit therefore
supports treating these as variant-level functional evidence from a calibrated
assay source, not merely as a domain-level rule. However, this does not by
itself prove that the local application strength is always correct. The next
step is to inspect discordant or unmatched records, and to verify how ENIGMA
Table 9 converts Findlay functional classes and RNA scores into `BS3`.

## Outputs

- `tables/bs3_functional_source_audit/brca1_bs3_table9_vs_findlay_sge.csv`
- `tables/bs3_functional_source_audit/brca1_bs3_findlay_match_summary.csv`
- `tables/bs3_functional_source_audit/brca1_bs3_findlay_function_class_summary.csv`
- `external_sources/brca1_bs3_functional_evidence_audit_sources/findlay_2018_brca1_sge/Findlay_2018_Nature_Supplementary_Table_1_SGE_scores.normalized.csv`
"""
    REPORT.write_text(text, encoding="utf-8")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
