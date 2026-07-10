from __future__ import annotations

import csv
import json
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
DATA_DIR = ANALYSIS_DIR / "data" / "structure_function_sources"
TABLE_DIR = ANALYSIS_DIR / "tables" / "structure_function_mapping"


UNIPROT = {
    "BRCA1": "P38398",
    "BRCA2": "P51587",
}


FIELDNAMES = [
    "gene",
    "protein_position",
    "annotation_type",
    "annotation_label",
    "evidence_source",
    "source_id",
    "structure_id",
    "curation_status",
    "notes",
]

REGION_FIELDNAMES = [
    "gene",
    "protein_start",
    "protein_end",
    "annotation_type",
    "annotation_label",
    "evidence_source",
    "source_id",
    "structure_id",
    "curation_status",
    "notes",
]

INTERACTION_TERMS = [
    "interaction",
    "affinity",
    "binding",
    "bind",
    "localization",
]


def fetch_json(url: str) -> dict[str, Any]:
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def feature_range(feature: dict[str, Any]) -> tuple[int | None, int | None]:
    location = feature.get("location", {})
    start = location.get("start", {}).get("value")
    end = location.get("end", {}).get("value")
    return start, end


def feature_description(feature: dict[str, Any]) -> str:
    return feature.get("description") or feature.get("type") or ""


def uniprot_feature_rows(gene: str, accession: str, payload: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for feature in payload.get("features", []):
        ftype = feature.get("type", "")
        start, end = feature_range(feature)
        description = feature_description(feature)
        if start is None:
            continue
        relevant_exact_feature = ftype in {"Binding site", "Metal binding", "Site"}
        relevant_functional_effect = (
            ftype in {"Mutagenesis", "Natural variant"}
            and any(term in description.lower() for term in INTERACTION_TERMS)
        )
        if start == end and (relevant_exact_feature or relevant_functional_effect):
            rows.append(
                {
                    "gene": gene,
                    "protein_position": str(start),
                    "annotation_type": f"UniProt:{ftype}",
                    "annotation_label": description,
                    "evidence_source": "UniProtKB/Swiss-Prot",
                    "source_id": accession,
                    "structure_id": "",
                    "curation_status": "source_extracted",
                    "notes": json.dumps(feature.get("evidences", []), ensure_ascii=True),
                }
            )
    return rows


def uniprot_region_rows(gene: str, accession: str, payload: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for feature in payload.get("features", []):
        ftype = feature.get("type", "")
        start, end = feature_range(feature)
        description = feature_description(feature)
        if start is None or end is None or start == end:
            continue
        if ftype not in {"Region", "Domain", "Zinc finger", "Motif"}:
            continue
        if ftype == "Region" and not any(term in description.lower() for term in INTERACTION_TERMS):
            continue
        if ftype in {"Domain", "Zinc finger", "Motif"} and not any(
            term in description.lower()
            for term in ["brct", "ring", "zinc", "palb2", "rad51", "sem1", "dss1"]
        ):
            continue
        rows.append(
            {
                "gene": gene,
                "protein_start": str(start),
                "protein_end": str(end),
                "annotation_type": f"UniProt:{ftype}",
                "annotation_label": description,
                "evidence_source": "UniProtKB/Swiss-Prot",
                "source_id": accession,
                "structure_id": "",
                "curation_status": "source_extracted_region",
                "notes": json.dumps(feature.get("evidences", []), ensure_ascii=True),
            }
        )
    return rows


def seed_literature_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for pos, aa in [
        (24, "Cys"),
        (27, "Cys"),
        (39, "Cys"),
        (41, "His"),
        (44, "Cys"),
        (47, "Cys"),
        (61, "Cys"),
        (64, "Cys"),
    ]:
        rows.append(
            {
                "gene": "BRCA1",
                "protein_position": str(pos),
                "annotation_type": "zinc_coordination_core",
                "annotation_label": f"BRCA1 RING C3HC4 zinc-coordinating residue {aa}{pos}",
                "evidence_source": "UniProt/PDB/literature_seed",
                "source_id": "P38398",
                "structure_id": "",
                "curation_status": "seed_needs_source_line_review",
                "notes": "RING C3HC4 motif seed; keep as source-tracked seed until exact source lines are reviewed.",
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, str]] = []
    manifest = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "sources": [],
    }
    for gene, accession in UNIPROT.items():
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        payload = fetch_json(url)
        snapshot = DATA_DIR / f"uniprot_{gene}_{accession}.json"
        snapshot.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        all_rows.extend(uniprot_feature_rows(gene, accession, payload))
        region_rows = uniprot_region_rows(gene, accession, payload)
        write_csv(TABLE_DIR / f"uniprot_{gene.lower()}_interface_regions.csv", region_rows, REGION_FIELDNAMES)
        manifest["sources"].append(
            {
                "gene": gene,
                "accession": accession,
                "url": url,
                "snapshot": str(snapshot.relative_to(ROOT)),
                "entry_version": payload.get("entryAudit", {}).get("entryVersion"),
                "sequence_version": payload.get("entryAudit", {}).get("sequenceVersion"),
            }
        )
    all_rows.extend(seed_literature_rows())
    deduped = {
        (
            row["gene"],
            row["protein_position"],
            row["annotation_type"],
            row["annotation_label"],
            row["evidence_source"],
        ): row
        for row in all_rows
    }
    rows = sorted(deduped.values(), key=lambda row: (row["gene"], int(row["protein_position"]), row["annotation_type"]))
    write_csv(TABLE_DIR / "curated_active_site_interface_annotations.csv", rows, FIELDNAMES)
    (DATA_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
