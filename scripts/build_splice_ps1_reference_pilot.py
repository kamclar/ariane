"""Build a pilot splice PS1 reference set from ENIGMA Supplementary Table 2.

This is a curation seed, not an automatically scored evidence source.
"""

from __future__ import annotations

import json
import re
import zipfile
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path


NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
}

PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PROJECT_ROOT.parents[1]
SUPPLEMENTARY_TABLES = WORKSPACE_ROOT / "SupplementaryTables_V1.2_2024-11-18.xlsx"
OUTPUT_PATH = PROJECT_ROOT / "backend" / "data" / "splice_ps1_reference_set.json"
SHEET_NAME = "ST2 splicing dataset codes"
CSPEC_URL = "https://cspec.genome.network/cspec/ui/svi/doc/GN097"
SPECIFICATIONS_SOURCE = "Specifications_V1.2_2024-11-18.docx"
APPENDICES_SOURCE = "Appendices_V1.2_2024-11-18.docx"


def _xml(zipf: zipfile.ZipFile, name: str) -> ET.Element:
    return ET.fromstring(zipf.read(name))


def _shared_strings(zipf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zipf.namelist():
        return []
    root = _xml(zipf, "xl/sharedStrings.xml")
    strings = []
    for si in root.findall("main:si", NS):
        strings.append("".join(node.text or "" for node in si.findall(".//main:t", NS)))
    return strings


def _sheets(zipf: zipfile.ZipFile) -> list[dict]:
    wb = _xml(zipf, "xl/workbook.xml")
    rels = _xml(zipf, "xl/_rels/workbook.xml.rels")
    rid_to_target = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall("pkgrel:Relationship", NS)
    }
    sheets = []
    for sheet in wb.findall("main:sheets/main:sheet", NS):
        rid = sheet.attrib[f"{{{NS['rel']}}}id"]
        sheets.append({
            "name": sheet.attrib["name"],
            "path": "xl/" + rid_to_target[rid].lstrip("/"),
        })
    return sheets


def _column_index(cell_ref: str) -> int:
    letters = re.match(r"([A-Z]+)", cell_ref).group(1)
    value = 0
    for char in letters:
        value = value * 26 + (ord(char) - ord("A") + 1)
    return value - 1


def _cell_text(cell: ET.Element, shared: list[str]) -> str:
    value = cell.find("main:v", NS)
    if value is None or value.text is None:
        inline = cell.find("main:is", NS)
        if inline is None:
            return ""
        return "".join(node.text or "" for node in inline.findall(".//main:t", NS)).strip()
    if cell.attrib.get("t") == "s":
        try:
            return shared[int(value.text)].strip()
        except (IndexError, ValueError):
            return ""
    return value.text.strip()


def _sheet_rows(path: Path, sheet_name: str) -> list[list[str]]:
    with zipfile.ZipFile(path) as zipf:
        shared = _shared_strings(zipf)
        sheet = next(s for s in _sheets(zipf) if s["name"] == sheet_name)
        root = _xml(zipf, sheet["path"])
        out = []
        for row in root.findall(".//main:sheetData/main:row", NS):
            values = []
            for cell in row.findall("main:c", NS):
                idx = _column_index(cell.attrib["r"])
                while len(values) <= idx:
                    values.append("")
                values[idx] = _cell_text(cell, shared)
            out.append(values)
        return out


def _classification_label(value: str) -> str:
    return {
        "5": "Pathogenic",
        "4": "Likely Pathogenic",
    }.get(value.strip(), "")


def _event_type(result: str) -> str:
    text = result.lower()
    if "no aberration" in text:
        return "no_aberration"
    if "deletion" in text or "skipping" in text:
        return "exon_skipping_or_deletion"
    if "insertion" in text:
        return "insertion"
    if "cryptic" in text:
        return "cryptic_splice_site"
    return "aberrant_splicing"


def build() -> dict:
    rows = _sheet_rows(SUPPLEMENTARY_TABLES, SHEET_NAME)
    header = rows[4]
    variants = []
    seen = set()
    for excel_row, row in enumerate(rows[5:], start=6):
        if len(row) < len(header):
            row += [""] * (len(header) - len(row))
        record = dict(zip(header, row))
        final_class = record.get(
            "Final Class (Based on multifactorial likelihood analysis). See Supp T7 for details.",
            "",
        ).strip()
        classification = _classification_label(final_class)
        result = record.get("Result", "").strip()
        if not classification or not result or "no aberration" in result.lower():
            continue
        key = (record.get("Gene", ""), record.get("HGVS c. nomenclature", ""))
        if key in seen:
            continue
        seen.add(key)
        variants.append({
            "gene": record.get("Gene", "").strip(),
            "reference_variant": record.get("HGVS c. nomenclature", "").strip(),
            "p_notation": record.get("HGVS p. nomenclature", "").strip(),
            "classification": classification,
            "classification_numeric": int(final_class),
            "classification_source": "ENIGMA Supplementary Table 2 / ST7 multifactorial likelihood class",
            "splice_event_label": result,
            "event_type": _event_type(result),
            "assay_result_category": record.get("Splicing Assay Result Category", "").strip(),
            "variant_context": record.get(
                "Summary of variant type, location, assay overview", ""
            ).strip(),
            "prior_probability": record.get(
                "Prior probability (from combined missense/ splicing predictions)", ""
            ).strip(),
            "included_in_analysis": record.get("Included in analysis", "").strip(),
            "transcript": "",
            "source_table": "SupplementaryTables_V1.2_2024-11-18.xlsx / ST2 splicing dataset codes",
            "source_file": "SupplementaryTables_V1.2_2024-11-18.xlsx",
            "source_sheet": SHEET_NAME,
            "source_row": excel_row,
            "rule_source": (
                "ENIGMA BRCA1/2 VCEP v1.2 Specifications PS1; "
                "Appendix J Table 17 PS1(splicing) weighting"
            ),
            "source_url": CSPEC_URL,
            "curation_status": "pilot_from_ST2_requires_review",
            "curation_note": (
                "Seed record only. Before PS1(splicing) use, confirm that the VUA "
                "has the same predicted event and similar/higher prediction strength "
                "than this reference variant, and apply Appendix J/Table 17 weighting."
            ),
        })

    return {
        "schema_version": "0.1",
        "generated_at": date.today().isoformat(),
        "curation_status": "pilot_unreviewed_not_for_automatic_scoring",
        "description": (
            "Pilot reference candidates for ENIGMA PS1(splicing). These records are "
            "not automatically scored and require expert review before use."
        ),
        "sources": [
            {
                "name": "ENIGMA BRCA1/2 VCEP Specifications v1.2",
                "version_date": "2024-11-18",
                "detail": "PS1 variable-weight splicing rule and Table 5/Table 17 weighting",
                "file": SPECIFICATIONS_SOURCE,
                "appendix_file": APPENDICES_SOURCE,
                "url": CSPEC_URL,
            },
            {
                "name": "ENIGMA Supplementary Table 2",
                "version_date": "2024-11-18",
                "detail": "mRNA assay result dataset with final multifactorial class",
                "file": "SupplementaryTables_V1.2_2024-11-18.xlsx",
                "sheet": SHEET_NAME,
                "path": str(SUPPLEMENTARY_TABLES),
            },
        ],
        "usage_requirements": [
            "The VUA predicted event must precisely match the reference splice event.",
            "The VUA prediction strength must be similar to or higher than the reference prediction strength.",
            "Reference variant classification should be assigned using VCEP specifications.",
            "For exonic VUAs, consider any predicted or proven missense/protein functional effect before applying PS1(splicing).",
            "Use Appendix J Table 17 to determine PS1(splicing) strength from VUA baseline code and relative position to the reference variant.",
        ],
        "variants": variants,
    }


def main() -> None:
    data = build()
    OUTPUT_PATH.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    counts = {}
    for variant in data["variants"]:
        counts[variant["classification"]] = counts.get(variant["classification"], 0) + 1
    print(json.dumps({
        "output": str(OUTPUT_PATH),
        "variant_count": len(data["variants"]),
        "counts": counts,
    }, indent=2))


if __name__ == "__main__":
    main()
