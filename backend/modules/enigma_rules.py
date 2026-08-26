"""Validated, public ENIGMA rules catalogue and table access.

The browser receives only curated metadata and paginated runtime records. Local
filesystem paths and the source workbooks themselves are never exposed.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import json
import re
from typing import Any

from backend.config import (
    ENIGMA_REFERENCE_TABLES_PATH,
    ENIGMA_RULE_CATALOG_PATH,
    ENIGMA_RULE_DIAGRAMS_PATH,
    TABLE9_PATH,
)
from backend.gene_policy import active_genes


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_PUBLIC_TABLE9_FIELDS = {
    "gene",
    "c_notation",
    "p_notation",
    "code",
    "strength",
    "text",
    "splice_result_published",
    "spliceai_prediction",
    "predicted_or_observed_splicing",
    "publication_count",
}

_TABLE_USAGE_ROLES = {
    "used_by_ariane": {
        "title": "Used by ARIANE",
        "description": "ENIGMA rule definitions, runtime lookups and approved candidate registries used by ARIANE.",
        "usages": {"rule_definition", "runtime_rule", "runtime_lookup", "candidate_registry"},
    },
    "expert_review": {
        "title": "Expert review",
        "description": "ENIGMA tables that support criteria requiring expert assessment.",
        "usages": {"expert_review"},
    },
    "supporting_reference": {
        "title": "Supporting and calibration",
        "description": "Official supporting evidence and calibration datasets that are not queried as standalone criteria.",
        "usages": {"rule_support", "calibration_reference"},
    },
}

_ARIANE_TABLE_CONSUMERS = {
    "specification-table-1": "Criterion definitions and gene-specific applicability",
    "specification-table-2": "Evidence-strength calibration",
    "specification-table-3": "ENIGMA criterion-combination policy",
    "specification-table-4": "PVS1 and PM5 PTC lookup",
    "specification-table-7": "PP4 and BP5 clinical-LR interpretation",
    "specification-table-9": "PS3, BS3 and curated splice lookup",
    "appendix-table-3": "BRCA1 functional-domain rules",
    "appendix-table-4": "BRCA2 functional-domain rules",
    "appendix-table-9": "RNA evidence weighting",
    "appendix-table-11": "Variant-specific moderate-risk and reduced-penetrance annotations",
    "appendix-table-14": "Missense prior-probability rules",
    "appendix-table-15": "Missense predictor calibration",
    "appendix-table-16": "Bioinformatic and domain decision rules",
    "supplementary-table-2": "Curated mRNA assay lookup",
    "supplementary-table-7": "Protein PS1 candidate-reference discovery",
}


def _table_role(usage: str) -> str:
    matches = [role for role, definition in _TABLE_USAGE_ROLES.items() if usage in definition["usages"]]
    if len(matches) != 1:
        raise RuntimeError(f"Unknown or ambiguous ENIGMA table usage: {usage}")
    return matches[0]


def _load_json(path) -> dict[str, Any]:
    if not path.is_file():
        raise RuntimeError(f"Required rules dataset is missing: {path.name}")
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise RuntimeError(f"Rules dataset must contain a JSON object: {path.name}")
    return payload


@lru_cache(maxsize=1)
def _load_reference_tables() -> dict[str, Any]:
    payload = _load_json(ENIGMA_REFERENCE_TABLES_PATH)
    if payload.get("schema_version") != 1:
        raise RuntimeError("Unsupported ENIGMA reference-table schema")
    if payload.get("criteria_specification_version") != "1.2.0":
        raise RuntimeError("Reference-table version does not match ENIGMA v1.2")
    tables = payload.get("tables")
    if not isinstance(tables, list) or payload.get("table_count") != 42:
        raise RuntimeError("ENIGMA reference-table inventory must contain 42 tables")
    table_ids = [table.get("id") for table in tables]
    if None in table_ids or len(table_ids) != len(set(table_ids)):
        raise RuntimeError("ENIGMA reference-table inventory contains duplicate ids")
    for table in tables:
        if table.get("group") not in {"specification", "appendix", "supplementary"}:
            raise RuntimeError(f"Invalid table group for {table.get('id')}")
        _table_role(str(table.get("usage") or ""))
        sections = table.get("sections")
        if not isinstance(sections, list) or not sections:
            raise RuntimeError(f"Reference table has no sections: {table.get('id')}")
        section_ids = [section.get("id") for section in sections]
        if None in section_ids or len(section_ids) != len(set(section_ids)):
            raise RuntimeError(f"Duplicate section id in {table.get('id')}")
        for section in sections:
            rows = section.get("rows")
            if not isinstance(rows, list) or section.get("row_count") != len(rows):
                raise RuntimeError(f"Invalid row count in {table.get('id')}/{section.get('id')}")
            width = max((len(row) for row in rows), default=0)
            if section.get("column_count") != width:
                raise RuntimeError(f"Invalid column count in {table.get('id')}/{section.get('id')}")
    used_ids = {
        table["id"] for table in tables
        if _table_role(table["usage"]) == "used_by_ariane"
    }
    if used_ids != set(_ARIANE_TABLE_CONSUMERS):
        missing = sorted(used_ids - set(_ARIANE_TABLE_CONSUMERS))
        stale = sorted(set(_ARIANE_TABLE_CONSUMERS) - used_ids)
        raise RuntimeError(
            f"ENIGMA table consumer registry mismatch: missing={missing}, stale={stale}"
        )
    return payload


def _public_reference_table_catalog() -> dict[str, Any]:
    payload = _load_reference_tables()
    roles = [
        {
            "id": role,
            "title": definition["title"],
            "description": definition["description"],
            "table_count": sum(
                1 for table in payload["tables"] if _table_role(table["usage"]) == role
            ),
        }
        for role, definition in _TABLE_USAGE_ROLES.items()
    ]
    return {
        "description": payload["view_description"],
        "groups": deepcopy(payload["groups"]),
        "roles": roles,
        "table_count": payload["table_count"],
        "items": [
            {
                key: deepcopy(table.get(key))
                for key in (
                    "id", "group", "number", "title", "source_id",
                    "source_sha256", "usage", "usage_text",
                )
            } | {
                "role": _table_role(table["usage"]),
                "consumer": _ARIANE_TABLE_CONSUMERS.get(table["id"]),
                "sections": [
                    {
                        key: section.get(key)
                        for key in ("id", "title", "row_count", "column_count")
                    }
                    for section in table["sections"]
                ],
            }
            for table in payload["tables"]
        ],
    }


def _validate_tree(tree: dict[str, Any], source_ids: set[str]) -> None:
    if tree.get("source_id") not in source_ids:
        raise RuntimeError(f"Unknown source for decision tree {tree.get('id')}")
    branch_ids: set[str] = set()
    for branch in tree.get("branches", []):
        branch_id = branch.get("id")
        if not branch_id or branch_id in branch_ids:
            raise RuntimeError(f"Duplicate or missing branch id in {tree.get('id')}")
        branch_ids.add(branch_id)
        nodes = branch.get("nodes", [])
        node_ids = {node.get("id") for node in nodes}
        if None in node_ids or len(node_ids) != len(nodes):
            raise RuntimeError(f"Duplicate or missing node id in {tree.get('id')}/{branch_id}")
        if branch.get("entry_node") not in node_ids:
            raise RuntimeError(f"Invalid entry node in {tree.get('id')}/{branch_id}")
        for edge in branch.get("edges", []):
            if edge.get("from") not in node_ids or edge.get("to") not in node_ids:
                raise RuntimeError(f"Invalid edge in {tree.get('id')}/{branch_id}")
        layout = branch.get("layout", {})
        positions = layout.get("positions", {})
        if set(positions) != node_ids:
            raise RuntimeError(f"Incomplete graph layout in {tree.get('id')}/{branch_id}")
        if not all(
            isinstance(position, list) and len(position) in {2, 4}
            and all(isinstance(value, (int, float)) for value in position)
            for position in positions.values()
        ):
            raise RuntimeError(f"Invalid graph coordinates in {tree.get('id')}/{branch_id}")
        for edge in branch.get("edges", []):
            points = edge.get("points")
            if points is not None and not (
                isinstance(points, list)
                and len(points) >= 2
                and all(
                    isinstance(point, list)
                    and len(point) == 2
                    and all(isinstance(value, (int, float)) for value in point)
                    for point in points
                )
            ):
                raise RuntimeError(f"Invalid edge geometry in {tree.get('id')}/{branch_id}")
            label_at = edge.get("label_at")
            if label_at is not None and not (
                isinstance(label_at, list)
                and len(label_at) == 2
                and all(isinstance(value, (int, float)) for value in label_at)
            ):
                raise RuntimeError(f"Invalid edge label coordinates in {tree.get('id')}/{branch_id}")

    provenance = tree.get("diagram_provenance")
    if provenance not in {"official_redraw", "ariane_derived"}:
        raise RuntimeError(f"Invalid diagram provenance for {tree.get('id')}")
    if provenance == "official_redraw" and not tree.get("original_figure_ids"):
        raise RuntimeError(f"Official redraw {tree.get('id')} has no original figure reference")


@lru_cache(maxsize=1)
def load_rule_catalog() -> dict[str, Any]:
    catalog = _load_json(ENIGMA_RULE_CATALOG_PATH)
    if catalog.get("schema_version") != 1:
        raise RuntimeError("Unsupported ENIGMA rule catalogue schema")
    sources = catalog.get("sources")
    if not isinstance(sources, list) or not sources:
        raise RuntimeError("ENIGMA rule catalogue has no sources")
    source_ids: set[str] = set()
    for source in sources:
        source_id = source.get("id")
        if not source_id or source_id in source_ids:
            raise RuntimeError("ENIGMA rule catalogue contains duplicate source ids")
        source_ids.add(source_id)
        if not str(source.get("official_url", "")).startswith("https://"):
            raise RuntimeError(f"Source {source_id} has no official HTTPS URL")
        if not _SHA256_RE.fullmatch(str(source.get("sha256", ""))):
            raise RuntimeError(f"Source {source_id} has an invalid SHA256")
        if any(key in source for key in ("path", "local_path", "file_path")):
            raise RuntimeError(f"Source {source_id} exposes a local path")
    source_by_id = {source["id"]: source for source in sources}
    for table in _load_reference_tables()["tables"]:
        source = source_by_id.get(table.get("source_id"))
        if source is None:
            raise RuntimeError(f"Unknown source for reference table {table.get('id')}")
        if table.get("source_sha256") != source.get("sha256"):
            raise RuntimeError(f"Source checksum mismatch for reference table {table.get('id')}")
    diagram_payload = _load_json(ENIGMA_RULE_DIAGRAMS_PATH)
    if diagram_payload.get("schema_version") != 1:
        raise RuntimeError("Unsupported ENIGMA rule diagram schema")
    merged_trees = [*catalog.get("decision_trees", []), *diagram_payload.get("decision_trees", [])]
    tree_ids = [tree.get("id") for tree in merged_trees]
    if None in tree_ids or len(tree_ids) != len(set(tree_ids)):
        raise RuntimeError("ENIGMA rule catalogue contains duplicate decision tree ids")
    catalog["decision_trees"] = merged_trees
    for tree in merged_trees:
        _validate_tree(tree, source_ids)
    for figure in catalog.get("figures", []):
        if figure.get("source_id") not in source_ids:
            raise RuntimeError(f"Unknown source for figure {figure.get('id')}")
        asset_url = str(figure.get("asset_url", ""))
        if not asset_url.startswith("/static/enigma/") or ".." in asset_url:
            raise RuntimeError(f"Invalid public asset URL for figure {figure.get('id')}")
    figure_ids = {figure.get("id") for figure in catalog.get("figures", [])}
    if None in figure_ids or len(figure_ids) != len(catalog.get("figures", [])):
        raise RuntimeError("ENIGMA rule catalogue contains duplicate or missing figure ids")
    for tree in merged_trees:
        unknown_figures = set(tree.get("original_figure_ids", [])) - figure_ids
        if unknown_figures:
            raise RuntimeError(
                f"Decision tree {tree.get('id')} references unknown figures: "
                f"{sorted(unknown_figures)}"
            )
    return catalog


@lru_cache(maxsize=1)
def _load_table9() -> dict[str, Any]:
    payload = _load_json(TABLE9_PATH)
    if payload.get("version") != load_rule_catalog().get("criteria_specification_version"):
        raise RuntimeError("Table 9 version does not match the ENIGMA rule catalogue")
    variants = payload.get("variants")
    if not isinstance(variants, dict) or payload.get("row_count") != len(variants):
        raise RuntimeError("Table 9 row count does not match its runtime records")
    expected = next(
        source for source in load_rule_catalog()["sources"]
        if source["id"] == "enigma-v1.2-table9"
    ).get("runtime_rows")
    if expected != len(variants):
        raise RuntimeError("Table 9 row count does not match the source registry")
    return payload


def public_catalog() -> dict[str, Any]:
    catalog = load_rule_catalog()
    return {
        "schema_version": catalog["schema_version"],
        "catalog_version": catalog["catalog_version"],
        "criteria_specification_id": catalog["criteria_specification_id"],
        "criteria_specification_version": catalog["criteria_specification_version"],
        "release_date": catalog["release_date"],
        "registry_url": catalog["registry_url"],
        "sources": deepcopy(catalog["sources"]),
        "tables": _public_reference_table_catalog(),
        "figures": deepcopy(catalog.get("figures", [])),
        "decision_trees": [
            {
                "id": tree["id"],
                "title": tree["title"],
                "description": tree["description"],
                "source_id": tree["source_id"],
                "source_location": tree["source_location"],
                "diagram_provenance": tree["diagram_provenance"],
                "original_figure_ids": deepcopy(tree.get("original_figure_ids", [])),
                "branches": [
                    {"id": branch["id"], "title": branch["title"]}
                    for branch in tree["branches"]
                ],
            }
            for tree in catalog["decision_trees"]
        ],
    }


def _public_cell_text(cell: Any) -> str:
    if isinstance(cell, dict):
        value = cell.get("value")
        return str(value if value not in {None, ""} else cell.get("formula") or "")
    return "" if cell is None else str(cell)


def _column_label(index: int) -> str:
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(65 + remainder) + result
    return result


def search_reference_table(
    table_id: str,
    *,
    section_id: str | None = None,
    query: str = "",
    page: int = 1,
    page_size: int = 25,
) -> dict[str, Any] | None:
    table = next(
        (item for item in _load_reference_tables()["tables"] if item["id"] == table_id),
        None,
    )
    if table is None:
        return None
    section = next(
        (item for item in table["sections"] if item["id"] == section_id),
        None,
    ) if section_id else table["sections"][0]
    if section is None:
        return None

    needle = query.strip().casefold()
    indexed_rows = list(enumerate(section["rows"], 1))
    if needle:
        indexed_rows = [
            (source_row, row)
            for source_row, row in indexed_rows
            if needle in " ".join(_public_cell_text(cell) for cell in row).casefold()
        ]
    total = len(indexed_rows)
    start = (page - 1) * page_size
    selected = indexed_rows[start:start + page_size]
    return {
        "table_id": table["id"],
        "table_title": table["title"],
        "section_id": section["id"],
        "section_title": section["title"],
        "source_id": table["source_id"],
        "page": page,
        "page_size": page_size,
        "total": total,
        "column_count": section["column_count"],
        "columns": [_column_label(index) for index in range(section["column_count"])],
        "items": [
            {"source_row": source_row, "cells": deepcopy(row)}
            for source_row, row in selected
        ],
    }


@lru_cache(maxsize=1)
def _reduced_penetrance_annotations() -> dict[tuple[str, str], dict[str, Any]]:
    """Extract explicit reduced-penetrance assertions from Appendix Table 11.

    This is deliberately narrower than a generic moderate-risk inference. A
    record is exposed only when the official ENIGMA row itself calls the allele
    a proven reduced-penetrance allele.
    """
    table = next(
        (
            item
            for item in _load_reference_tables()["tables"]
            if item["id"] == "appendix-table-11"
        ),
        None,
    )
    if table is None:
        raise RuntimeError("ENIGMA Appendix Table 11 is missing")
    rows = table["sections"][0]["rows"]
    if not rows:
        raise RuntimeError("ENIGMA Appendix Table 11 is empty")
    header = [_public_cell_text(cell).strip() for cell in rows[0]]
    required = {
        "Gene",
        "Variant",
        "Other information relating to risk",
    }
    if not required.issubset(header):
        raise RuntimeError("ENIGMA Appendix Table 11 columns are incomplete")
    gene_index = header.index("Gene")
    variant_index = header.index("Variant")
    risk_index = header.index("Other information relating to risk")
    source = next(
        item
        for item in load_rule_catalog()["sources"]
        if item["id"] == table["source_id"]
    )

    result: dict[tuple[str, str], dict[str, Any]] = {}
    for source_row, row in enumerate(rows[1:], 2):
        padded = list(row) + [None] * max(0, len(header) - len(row))
        risk_text = _public_cell_text(padded[risk_index]).strip()
        if "proven reduced penetrance allele" not in risk_text.casefold():
            continue
        gene = _public_cell_text(padded[gene_index]).strip().upper()
        variant_text = re.sub(
            r"\s+", "", _public_cell_text(padded[variant_index]).strip()
        )
        match = re.search(r"c\.\d+(?:[+-]\d+)?[ACGT]>[ACGT]", variant_text, re.I)
        if gene not in set(active_genes()) or match is None:
            raise RuntimeError(
                f"Could not normalize reduced-penetrance row {source_row} in Appendix Table 11"
            )
        c_notation = match.group(0)
        pmids = sorted(set(re.findall(r"PMID\s*:\s*(\d+)", risk_text, re.I)))
        key = (gene, c_notation)
        if key in result:
            raise RuntimeError(
                f"Duplicate reduced-penetrance annotation in Appendix Table 11: {gene} {c_notation}"
            )
        result[key] = {
            "category": "reduced_penetrance",
            "label": "Variant with reduced penetrance",
            "summary": (
                "ENIGMA identifies this allele as having proven reduced penetrance. "
                "This is a variant-specific clinical-risk annotation and does not "
                "change the criteria, points or class calculated by ARIANE."
            ),
            "evidence": risk_text,
            "source": "ENIGMA BRCA1/2 VCEP v1.2 Appendix Table 11",
            "source_url": source["official_url"],
            "source_row": source_row,
            "publications": [
                {
                    "pmid": pmid,
                    "label": f"PMID {pmid}",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                }
                for pmid in pmids
            ],
            "affects_classification": False,
        }
    if not result:
        raise RuntimeError(
            "ENIGMA Appendix Table 11 contains no explicit reduced-penetrance annotation"
        )
    return result


def clinical_annotations_for_variant(
    gene: str,
    c_notation: str,
) -> list[dict[str, Any]]:
    annotation = _reduced_penetrance_annotations().get(
        (gene.strip().upper(), c_notation.strip())
    )
    return [deepcopy(annotation)] if annotation else []


def get_decision_tree(tree_id: str) -> dict[str, Any] | None:
    for tree in load_rule_catalog().get("decision_trees", []):
        if tree.get("id") == tree_id:
            return deepcopy(tree)
    return None


def _public_table9_record(record: dict[str, Any]) -> dict[str, Any]:
    public = {key: record.get(key) for key in _PUBLIC_TABLE9_FIELDS}
    public["results"] = [
        record.get(f"result_{index}") for index in range(1, 5)
        if record.get(f"result_{index}")
    ]
    return public


def search_table9(
    *,
    gene: str | None = None,
    query: str = "",
    code: str | None = None,
    page: int = 1,
    page_size: int = 25,
) -> dict[str, Any]:
    records = _load_table9()["variants"].values()
    normalized_gene = gene.upper() if gene else None
    normalized_code = code.upper() if code else None
    needle = query.strip().casefold()
    filtered: list[dict[str, Any]] = []
    for record in records:
        if normalized_gene and record.get("gene") != normalized_gene:
            continue
        if normalized_code and str(record.get("code", "")).upper() != normalized_code:
            continue
        if needle:
            searchable = " ".join(str(record.get(field) or "") for field in (
                "gene", "c_notation", "p_notation", "code", "strength", "text"
            )).casefold()
            if needle not in searchable:
                continue
        filtered.append(record)
    filtered.sort(key=lambda item: (item.get("gene", ""), item.get("c_notation", "")))
    total = len(filtered)
    start = (page - 1) * page_size
    return {
        "table_id": "table9",
        "version": _load_table9()["version"],
        "source_id": "enigma-v1.2-table9",
        "page": page,
        "page_size": page_size,
        "total": total,
        "items": [_public_table9_record(item) for item in filtered[start:start + page_size]],
    }


def validate_rule_catalog() -> None:
    load_rule_catalog()
    _load_table9()
    _load_reference_tables()
    _reduced_penetrance_annotations()
