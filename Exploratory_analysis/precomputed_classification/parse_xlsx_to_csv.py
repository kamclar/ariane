"""Minimal XLSX-to-CSV extractor for downloaded source workbooks."""

from __future__ import annotations

import csv
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def column_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    index = 0
    for ch in letters:
        index = index * 26 + (ord(ch.upper()) - ord("A") + 1)
    return index - 1


def load_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    values: list[str] = []
    for si in root.findall("a:si", NS):
        texts = [node.text or "" for node in si.findall(".//a:t", NS)]
        values.append("".join(texts))
    return values


def cell_value(cell: ET.Element, shared: list[str]) -> str:
    value_node = cell.find("a:v", NS)
    if value_node is None:
        inline = cell.find("a:is", NS)
        if inline is None:
            return ""
        return "".join(node.text or "" for node in inline.findall(".//a:t", NS))
    value = value_node.text or ""
    if cell.attrib.get("t") == "s":
        try:
            return shared[int(value)]
        except (ValueError, IndexError):
            return value
    return value


def parse_sheet(path: Path) -> list[list[str]]:
    with zipfile.ZipFile(path) as zf:
        shared = load_shared_strings(zf)
        root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        rows: list[list[str]] = []
        for row in root.findall(".//a:sheetData/a:row", NS):
            values: dict[int, str] = {}
            for cell in row.findall("a:c", NS):
                ref = cell.attrib.get("r", "")
                if not ref:
                    continue
                values[column_index(ref)] = cell_value(cell, shared)
            if values:
                max_col = max(values)
                rows.append([values.get(i, "") for i in range(max_col + 1)])
        return rows


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: parse_xlsx_to_csv.py input.xlsx output.csv")
    source = Path(sys.argv[1])
    target = Path(sys.argv[2])
    rows = parse_sheet(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)
    print(f"Wrote {target} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
