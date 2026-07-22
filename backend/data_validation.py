"""Fail-fast validation for datasets required by automated classification."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping


def _load_required_json(label: str, path: Path):
    if not path.is_file():
        raise RuntimeError(f"Required {label} dataset is missing: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Required {label} dataset cannot be loaded: {path}: {exc}") from exc


def validate_required_datasets(paths: Mapping[str, Path]) -> None:
    """Raise RuntimeError before the API starts if a core dataset is unusable."""
    table4 = _load_required_json("ENIGMA Table 4", paths["table4"])
    if (
        table4.get("schema_version") != 2
        or table4.get("source_columns") != 20
        or len(table4.get("source_rows", [])) != 493
    ):
        raise RuntimeError(
            "Required ENIGMA Table 4 dataset is not the complete lossless v1.2 snapshot"
        )
    required_table4_sections = {
        "exon_ranges", "ptc_rules", "splice_rules", "deletion_rules", "duplication_rules"
    }
    missing_sections = sorted(required_table4_sections - set(table4))
    if missing_sections:
        raise RuntimeError(f"Required ENIGMA Table 4 dataset is incomplete: missing {missing_sections}")
    for section in required_table4_sections:
        for gene in ("BRCA1", "BRCA2"):
            if not table4.get(section, {}).get(gene):
                raise RuntimeError(f"Required ENIGMA Table 4 dataset has no {section}/{gene} records")
    for gene in ("BRCA1", "BRCA2"):
        exon_ranges = set(table4["exon_ranges"][gene])
        for section in ("ptc_rules", "deletion_rules", "duplication_rules"):
            missing_ranges = sorted(set(table4[section][gene]) - exon_ranges)
            if missing_ranges:
                raise RuntimeError(
                    f"Required ENIGMA Table 4 dataset is inconsistent: "
                    f"{section}/{gene} exons lack ranges: {missing_ranges}"
                )
        for key, entry in table4["splice_rules"][gene].items():
            if entry.get("exon") not in exon_ranges:
                raise RuntimeError(
                    f"Required ENIGMA Table 4 dataset is inconsistent: "
                    f"splice_rules/{gene}/{key} references unknown exon {entry.get('exon')!r}"
                )

    allowed_pvs1_codes = {
        "PVS1", "PVS1 (RNA)", "PVS1_Strong", "PVS1_Strong (RNA)",
        "PVS1_Moderate", "PVS1_Moderate (RNA)", "PVS1_Supporting",
        "PVS1_Supporting (RNA)", "PVS1_N/A", "PVS1_N/A (RNA)",
    }
    allowed_pm5_codes = {
        None, "PM5_Strong (PTC)", "PM5 (PTC)",
        "PM5_Supporting (PTC)", "PM5_N/A",
    }
    for section in ("ptc_rules", "splice_rules", "deletion_rules"):
        for gene in ("BRCA1", "BRCA2"):
            for key, entry in table4[section][gene].items():
                if entry.get("pvs1_code") not in allowed_pvs1_codes:
                    raise RuntimeError(
                        f"Required ENIGMA Table 4 dataset has unknown PVS1 code at "
                        f"{section}/{gene}/{key}: {entry.get('pvs1_code')!r}"
                    )
                if section == "ptc_rules" and entry.get("pm5_code") not in allowed_pm5_codes:
                    raise RuntimeError(
                        f"Required ENIGMA Table 4 dataset has unknown PM5 code at "
                        f"{section}/{gene}/{key}: {entry.get('pm5_code')!r}"
                    )
    for gene in ("BRCA1", "BRCA2"):
        for exon, arrangements in table4["duplication_rules"][gene].items():
            for arrangement, entry in arrangements.items():
                if arrangement not in {"Unknown", "Tandem"}:
                    raise RuntimeError(
                        f"Required ENIGMA Table 4 dataset has unknown duplication arrangement "
                        f"at {gene}/{exon}: {arrangement!r}"
                    )
                if entry.get("pvs1_code") not in allowed_pvs1_codes:
                    raise RuntimeError(
                        f"Required ENIGMA Table 4 dataset has unknown PVS1 code at "
                        f"duplication_rules/{gene}/{exon}/{arrangement}: {entry.get('pvs1_code')!r}"
                    )

    table9 = _load_required_json("ENIGMA Table 9", paths["table9"])
    if not isinstance(table9.get("variants"), dict) or not table9["variants"]:
        raise RuntimeError("Required ENIGMA Table 9 dataset has no variant records")
    if table9.get("schema_version") != 2 or table9.get("row_count") != 4731:
        raise RuntimeError(
            "Required ENIGMA Table 9 dataset is not the complete lossless v1.2 snapshot"
        )
    allowed_functional = {
        ("PS3", "Strong"), ("PS3", "Moderate"), ("PS3", "Supporting"),
        ("BS3", "Strong"), ("BS3", "Moderate"), ("BS3", "Supporting"),
        ("None", "N/A"),
    }
    required_table9_fields = {
        "gene", "c_notation", "p_notation", "code", "strength", "text",
        "splice_result_published", "spliceai_prediction",
        "predicted_or_observed_splicing", "publication_count",
        "result_1", "result_2", "result_3", "result_4",
    }
    for key, entry in table9["variants"].items():
        missing_fields = required_table9_fields - set(entry)
        if missing_fields:
            raise RuntimeError(
                f"Required ENIGMA Table 9 dataset is lossy at {key}: "
                f"missing {sorted(missing_fields)}"
            )
        if (entry.get("code"), entry.get("strength")) not in allowed_functional:
            raise RuntimeError(
                f"Required ENIGMA Table 9 dataset has unsupported code/strength at "
                f"{key}: {entry.get('code')!r}/{entry.get('strength')!r}"
            )
        if not key.startswith(("BRCA1:c.", "BRCA2:c.")):
            raise RuntimeError(f"Required ENIGMA Table 9 dataset has invalid variant key: {key!r}")

    st7 = _load_required_json("ENIGMA Supplementary Table 7", paths["st7"])
    if (
        not isinstance(st7, dict)
        or st7.get("schema_version") != 2
        or st7.get("source_columns") != 28
        or st7.get("total_variants") != 773
    ):
        raise RuntimeError(
            "Required ENIGMA Supplementary Table 7 is not the complete lossless v1.2 snapshot"
        )
    records = st7.get("variants") if isinstance(st7, dict) else st7
    if not records:
        raise RuntimeError("Required ENIGMA Supplementary Table 7 dataset has no variant records")
    seen = set()
    for index, record in enumerate(records):
        if len(record) != 28:
            raise RuntimeError(
                f"Required ENIGMA Supplementary Table 7 is lossy at record #{index}: "
                f"expected 28 fields, found {len(record)}"
            )
        gene = record.get("gene")
        c_notation = record.get("c_notation")
        key = f"{gene}:{c_notation}"
        if gene not in {"BRCA1", "BRCA2"} or not str(c_notation).startswith("c."):
            raise RuntimeError(f"Required ENIGMA Supplementary Table 7 has invalid record #{index}: {key}")
        if key in seen:
            raise RuntimeError(f"Required ENIGMA Supplementary Table 7 has duplicate variant: {key}")
        seen.add(key)
        posterior = record.get("posterior_probability")
        if posterior is not None and (
            not isinstance(posterior, (int, float)) or not 0 <= posterior <= 1
        ):
            raise RuntimeError(
                f"Required ENIGMA Supplementary Table 7 has invalid posterior for {key}: {posterior!r}"
            )
        if record.get("iarc_class") not in {1, 2, 4, 5}:
            raise RuntimeError(
                f"Required ENIGMA Supplementary Table 7 has invalid IARC class for {key}: "
                f"{record.get('iarc_class')!r}"
            )
