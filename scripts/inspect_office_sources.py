"""Inspect local Office sources without third-party dependencies."""

from __future__ import annotations

import argparse
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path


NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
}


def _xml(zipf: zipfile.ZipFile, name: str) -> ET.Element:
    return ET.fromstring(zipf.read(name))


def _xlsx_shared_strings(zipf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zipf.namelist():
        return []
    root = _xml(zipf, "xl/sharedStrings.xml")
    strings = []
    for si in root.findall("main:si", NS):
        parts = [node.text or "" for node in si.findall(".//main:t", NS)]
        strings.append("".join(parts))
    return strings


def _xlsx_sheets(zipf: zipfile.ZipFile) -> list[dict]:
    wb = _xml(zipf, "xl/workbook.xml")
    rels = _xml(zipf, "xl/_rels/workbook.xml.rels")
    rid_to_target = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall("pkgrel:Relationship", NS)
    }
    out = []
    for sheet in wb.findall("main:sheets/main:sheet", NS):
        rid = sheet.attrib[f"{{{NS['rel']}}}id"]
        target = rid_to_target[rid]
        path = "xl/" + target.lstrip("/")
        out.append({"name": sheet.attrib["name"], "path": path})
    return out


def _cell_text(cell: ET.Element, shared: list[str]) -> str:
    value = cell.find("main:v", NS)
    if value is None or value.text is None:
        inline = cell.find("main:is", NS)
        if inline is None:
            return ""
        return "".join(node.text or "" for node in inline.findall(".//main:t", NS))
    if cell.attrib.get("t") == "s":
        try:
            return shared[int(value.text)]
        except (IndexError, ValueError):
            return ""
    return value.text


def _xlsx_rows(path: Path, sheet_name: str, start: int, end: int) -> dict:
    with zipfile.ZipFile(path) as zipf:
        shared = _xlsx_shared_strings(zipf)
        sheets = _xlsx_sheets(zipf)
        matches = [sheet for sheet in sheets if sheet["name"] == sheet_name]
        if not matches:
            return {"file": str(path), "sheet": sheet_name, "rows": []}
        root = _xml(zipf, matches[0]["path"])
        rows = []
        for row in root.findall(".//main:sheetData/main:row", NS):
            row_num = int(row.attrib.get("r", "0"))
            if row_num < start or row_num > end:
                continue
            rows.append({
                "row": row_num,
                "cells": [_cell_text(cell, shared) for cell in row.findall("main:c", NS)],
            })
        return {"file": str(path), "sheet": sheet_name, "rows": rows}


def inspect_xlsx(path: Path, pattern: str) -> dict:
    regex = re.compile(pattern, re.I)
    matches = []
    with zipfile.ZipFile(path) as zipf:
        shared = _xlsx_shared_strings(zipf)
        sheets = _xlsx_sheets(zipf)
        for sheet in sheets:
            root = _xml(zipf, sheet["path"])
            rows = root.findall(".//main:sheetData/main:row", NS)
            for row in rows:
                cells = [
                    _cell_text(cell, shared)
                    for cell in row.findall("main:c", NS)
                ]
                joined = " | ".join(c for c in cells if c)
                if joined and regex.search(joined):
                    matches.append({
                        "sheet": sheet["name"],
                        "row": row.attrib.get("r", ""),
                        "text": joined[:1200],
                    })
    return {"file": str(path), "sheets": sheets, "matches": matches}


def inspect_docx(path: Path, pattern: str) -> dict:
    regex = re.compile(pattern, re.I)
    with zipfile.ZipFile(path) as zipf:
        root = _xml(zipf, "word/document.xml")
    paragraphs = []
    for para in root.findall(".//w:p", NS):
        text = "".join(node.text or "" for node in para.findall(".//w:t", NS))
        if text.strip():
            paragraphs.append(text.strip())
    matches = []
    for idx, text in enumerate(paragraphs):
        if regex.search(text):
            context = paragraphs[max(0, idx - 2): idx + 3]
            matches.append({"paragraph": idx + 1, "text": "\n".join(context)[:1600]})
    return {"file": str(path), "paragraph_count": len(paragraphs), "matches": matches}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--pattern", default="PS1|Table 5|splic")
    parser.add_argument("--sheet")
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int, default=20)
    args = parser.parse_args()
    results = []
    for raw in args.paths:
        path = Path(raw)
        if path.suffix.lower() == ".xlsx" and args.sheet:
            results.append(_xlsx_rows(path, args.sheet, args.start, args.end))
        elif path.suffix.lower() == ".xlsx":
            results.append(inspect_xlsx(path, args.pattern))
        elif path.suffix.lower() == ".docx":
            results.append(inspect_docx(path, args.pattern))
    print(json.dumps(results, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
