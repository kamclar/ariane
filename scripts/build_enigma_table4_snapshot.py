"""Build a lossless, executable ENIGMA BRCA1/2 VCEP Table 4 snapshot."""
from __future__ import annotations

import argparse
import json
import re
from enum import IntEnum
from pathlib import Path
from typing import Any, TypeAlias

import openpyxl


SOURCE_URL = (
    "https://cspec.genome.network/cspec/File/id/"
    "ca5cf57b-94df-4ad6-a001-c62ceccb3845/data"
)
SNAPSHOT_DATE = "2026-11-10"


class Column(IntEnum):
    GENE = 0
    EXON = 1
    ACCEPTOR_NOTES = 2
    ACCEPTOR_CODE = 3
    ACCEPTOR_PATTERN = 4
    ACCEPTOR_ALLELES = 5
    REGION = 7
    EXON_START_FALLBACK = 8
    PRIMARY_VALUE = 9
    MARKER = 10
    SECONDARY_VALUE = 11
    EXON_END = 12
    DOMAIN_VALUE = 13
    DONOR_PATTERN = 15
    DONOR_ALLELES = 16
    DONOR_CODE = 17
    DONOR_NOTES = 18


EXPECTED_HEADERS = {
    Column.GENE: "Gene",
    Column.EXON: "Exon",
    Column.ACCEPTOR_NOTES: "Notes & Warnings",
    Column.ACCEPTOR_CODE: "Code for IVS variant",
    Column.ACCEPTOR_PATTERN: "IVS 1,2 Variants",
    Column.ACCEPTOR_ALLELES: "Vt Allele",
    Column.DONOR_PATTERN: "IVS 1,2 Variants",
    Column.DONOR_ALLELES: "Vt Allele",
    Column.DONOR_CODE: "Code for IVS variant",
    Column.DONOR_NOTES: "Notes & Warnings",
}

SourceRow: TypeAlias = list[Any]
Rule: TypeAlias = dict[str, Any]
GeneRules: TypeAlias = dict[str, dict[str, Rule]]
ContextNotes: TypeAlias = dict[str, dict[str, dict[str, list[str]]]]


# Appendix D records the gene-level initiation-codon branch that was used to
# derive the final Table 4 code.  Keep these source facts next to the generated
# initiation rule so the manual-review form can reproduce the curated branch
# without asking the reviewer to re-enter facts already fixed by ENIGMA.
APPENDIX_D_INITIATION_FACTS = {
    "BRCA1": {
        "alternative_start_assessed": "yes",
        "alternative_start_conclusion": "no_known_alternative_start",
        "nearest_alternative_start": "p.Met18",
        "upstream_pathogenic_evidence": "yes",
        "functional_domain_impact": "not_applicable",
        "appendix_figure": "Appendix D Figure 3",
    },
    "BRCA2": {
        "alternative_start_assessed": "yes",
        "alternative_start_conclusion": "no_known_alternative_start",
        "nearest_alternative_start": "p.Met124",
        "upstream_pathogenic_evidence": "yes",
        "functional_domain_impact": "not_applicable",
        "appendix_figure": "Appendix D Figure 5",
    },
}


def _clean(value: Any) -> Any:
    if isinstance(value, str):
        return value.strip()
    return value


def _validate_header(header: SourceRow) -> None:
    required_width = max(int(column) for column in EXPECTED_HEADERS) + 1
    if len(header) < required_width:
        raise RuntimeError(
            "Unexpected Table 4 workbook width: "
            f"expected at least {required_width} columns, found {len(header)}"
        )
    mismatches = {
        column.name: (expected, str(header[column]).strip())
        for column, expected in EXPECTED_HEADERS.items()
        if str(header[column]).strip() != expected
    }
    if mismatches:
        raise RuntimeError(f"Unexpected Table 4 workbook columns: {mismatches}")


def _position(value: Any) -> int | None:
    if value is None:
        return None
    match = re.fullmatch(r"c\.(-?\d+)\*?", str(value).strip())
    return int(match.group(1)) if match else None


def _extend_unique_notes(target: list[str], *values: Any) -> None:
    for value in values:
        if value in (None, ""):
            continue
        note = str(value).strip()
        if note and note not in target:
            target.append(note)


def _notes(*values: Any) -> str:
    notes: list[str] = []
    _extend_unique_notes(notes, *values)
    return " | ".join(notes)


def _with_notes(rule: Rule, notes: str) -> Rule:
    if notes:
        rule["notes"] = notes
    return rule


def _alleles(value: Any) -> list[str]:
    if value is None:
        return []
    return [part.strip() for part in str(value).split(",") if part.strip()]


def _context_notes(rows: list[SourceRow]) -> ContextNotes:
    notes: ContextNotes = {"BRCA1": {}, "BRCA2": {}}
    for row in rows[1:]:
        gene, exon = row[Column.GENE], row[Column.EXON]
        if gene not in notes or not exon:
            continue
        entry = notes[gene].setdefault(
            exon,
            {"acceptor": [], "donor": [], "exon": []},
        )
        if not (row[Column.ACCEPTOR_CODE] and row[Column.ACCEPTOR_PATTERN]):
            _extend_unique_notes(entry["acceptor"], row[Column.ACCEPTOR_NOTES])
        if not (row[Column.DONOR_PATTERN] and row[Column.DONOR_CODE]):
            _extend_unique_notes(entry["donor"], row[Column.DONOR_NOTES])

    for gene_notes in notes.values():
        for entry in gene_notes.values():
            shared = [
                note for note in entry["acceptor"] if note in entry["donor"]
            ]
            entry["exon"] = shared
            entry["acceptor"] = [
                note for note in entry["acceptor"] if note not in shared
            ]
            entry["donor"] = [
                note for note in entry["donor"] if note not in shared
            ]
    return notes


def _add_splice_rules(
    target: GeneRules,
    gene: str,
    exon: str,
    pattern: Any,
    alleles: Any,
    code: Any,
    notes: str,
) -> None:
    if not pattern or not code:
        return
    pattern = str(pattern).strip()
    if not pattern.startswith("c."):
        raise RuntimeError(f"Unexpected Table 4 splice pattern: {pattern!r}")
    base = pattern[:-1] if pattern.endswith(">") else pattern
    for allele in _alleles(alleles):
        key = base + allele
        if key in target[gene]:
            previous = target[gene][key]
            new_code = str(code).strip()
            # Table 4 deliberately lists a generic prediction and then a
            # variant-specific RNA result for a few alleles.  The latter is
            # the executable rule while both rows remain in source_rows.
            if "(RNA)" not in new_code and previous["pvs1_code"] != new_code:
                raise RuntimeError(f"Conflicting Table 4 splice rule: {gene}:{key}")
        target[gene][key] = {
            "exon": exon,
            "pvs1_code": str(code).strip(),
            "notes": _notes(notes),
        }


def build_snapshot(source: Path) -> dict[str, Any]:
    workbook = openpyxl.load_workbook(source, read_only=True, data_only=True)
    readme_sheet = workbook["Specifications - Table 4 Readme"]
    sheet = workbook["Table 4 - Annotated Exons"]
    rows = [[_clean(value) for value in row] for row in sheet.iter_rows(values_only=True)]
    _validate_header(rows[0])
    context_notes = _context_notes(rows)

    exon_ranges: dict[str, dict[str, list[int]]] = {"BRCA1": {}, "BRCA2": {}}
    initiation_rules: dict[str, Rule] = {"BRCA1": {}, "BRCA2": {}}
    ptc_rules: GeneRules = {"BRCA1": {}, "BRCA2": {}}
    splice_rules: GeneRules = {"BRCA1": {}, "BRCA2": {}}
    deletion_rules: GeneRules = {"BRCA1": {}, "BRCA2": {}}
    duplication_rules: GeneRules = {"BRCA1": {}, "BRCA2": {}}
    critical_boundaries: dict[str, Rule] = {}

    for row_number, row in enumerate(rows[1:], 2):
        gene, exon = row[Column.GENE], row[Column.EXON]
        if gene not in {"BRCA1", "BRCA2"} or not exon:
            continue

        start = next(
            (
                position
                for column in (
                    Column.PRIMARY_VALUE,
                    Column.EXON_START_FALLBACK,
                    Column.REGION,
                )
                if (position := _position(row[column])) is not None
            ),
            None,
        )
        end_positions = [
            position
            for column in (
                Column.SECONDARY_VALUE,
                Column.EXON_END,
                Column.DOMAIN_VALUE,
            )
            if (position := _position(row[column])) is not None
        ]
        end = max(end_positions) if end_positions else None
        if start is not None and end is not None:
            exon_ranges[gene][exon] = [start, end]

        marker_value = row[Column.MARKER]
        marker = str(marker_value).strip() if marker_value is not None else ""
        code = row[Column.PRIMARY_VALUE]
        note_context = context_notes[gene].get(
            exon,
            {"acceptor": [], "donor": [], "exon": []},
        )
        exon_notes = note_context["exon"]
        structural_notes = _notes(
            row[Column.ACCEPTOR_NOTES],
            row[Column.DONOR_NOTES],
            *exon_notes,
        )
        if marker == "START codon":
            if initiation_rules[gene]:
                raise RuntimeError(
                    f"Duplicate Table 4 initiation-codon rule for {gene} at row {row_number}"
                )
            initiation_rules[gene] = {
                "exon": exon,
                "pvs1_code": str(code).strip(),
                "notes": structural_notes,
                "source_row": row_number,
                **APPENDIX_D_INITIATION_FACTS[gene],
            }
        elif marker == "PTC":
            ptc_rules[gene][exon] = _with_notes(
                {
                    "pvs1_code": code,
                    "pm5_code": row[Column.SECONDARY_VALUE],
                    "condition": None,
                },
                structural_notes,
            )
        elif marker.startswith("PTC<"):
            boundary_match = re.search(r"p\.[A-Z][a-z]{0,2}(\d+)", marker)
            if not boundary_match:
                raise RuntimeError(f"Cannot parse PTC boundary at row {row_number}: {marker}")
            first_na_aa = int(boundary_match.group(1))
            critical_boundaries[gene] = _with_notes(
                {
                    "exon": exon,
                    "boundary_aa": first_na_aa - 1,
                    "rule": marker,
                    "at_or_before": {
                        "pvs1_code": code,
                        "pm5_code": row[Column.SECONDARY_VALUE],
                    },
                },
                structural_notes,
            )
        elif marker.startswith("PTC>"):
            ptc_rules[gene][exon] = _with_notes(
                {
                    "pvs1_code": code,
                    "pm5_code": row[Column.SECONDARY_VALUE],
                    "condition": marker,
                },
                structural_notes,
            )
            if gene in critical_boundaries:
                critical_boundaries[gene]["rule"] += f"; {marker}"
                boundary_notes = _notes(
                    critical_boundaries[gene].get("notes"),
                    structural_notes,
                )
                if boundary_notes:
                    critical_boundaries[gene]["notes"] = boundary_notes
                    ptc_rules[gene][exon]["notes"] = boundary_notes
                critical_boundaries[gene]["after"] = {
                    "pvs1_code": code,
                    "pm5_code": row[Column.SECONDARY_VALUE],
                }
        elif marker.startswith("DEL "):
            deletion_rules[gene][exon] = {
                "pvs1_code": code,
                "notes": structural_notes,
            }
        elif marker.startswith("DUP "):
            arrangement = str(row[Column.SECONDARY_VALUE]).strip()
            duplication_rules[gene].setdefault(exon, {})[arrangement] = {
                "pvs1_code": code,
                "notes": structural_notes,
            }

        _add_splice_rules(
            splice_rules, gene, exon,
            row[Column.ACCEPTOR_PATTERN],
            row[Column.ACCEPTOR_ALLELES],
            row[Column.ACCEPTOR_CODE],
            _notes(
                row[Column.ACCEPTOR_NOTES],
                *exon_notes,
                *note_context["acceptor"],
            ),
        )
        _add_splice_rules(
            splice_rules, gene, exon,
            row[Column.DONOR_PATTERN],
            row[Column.DONOR_ALLELES],
            row[Column.DONOR_CODE],
            _notes(
                row[Column.DONOR_NOTES],
                *exon_notes,
                *note_context["donor"],
            ),
        )

    for gene in ("BRCA1", "BRCA2"):
        expected = set(exon_ranges[gene])
        for name, section in (
            ("deletion", deletion_rules),
            ("duplication", duplication_rules),
        ):
            missing = expected - set(section[gene])
            if missing:
                raise RuntimeError(f"Missing {name} rules for {gene}: {sorted(missing)}")

    readme = [_clean(row[0]) for row in readme_sheet.iter_rows(values_only=True)]
    return {
        "schema_version": 2,
        "version": "1.2.0",
        "released": "2025-01-09",
        "source": "ClinGen ENIGMA BRCA1/2 VCEP Specifications Table 4",
        "source_url": SOURCE_URL,
        "generated": SNAPSHOT_DATE,
        "readme": readme,
        "source_columns": 20,
        "source_rows": rows,
        "exon_ranges": exon_ranges,
        "initiation_rules": initiation_rules,
        "ptc_rules": ptc_rules,
        "splice_rules": splice_rules,
        "deletion_rules": deletion_rules,
        "duplication_rules": duplication_rules,
        "critical_boundaries": critical_boundaries,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path("backend/data/enigma_table4.json")
    )
    args = parser.parse_args()
    snapshot = build_snapshot(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"Wrote lossless Table 4 snapshot with {len(snapshot['source_rows'])} source rows "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
