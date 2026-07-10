"""Query public ClinVar and ClinGen/ENIGMA assertions for priority VUS."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
ANALYSIS_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
INPUT = ANALYSIS_DIR / "tables" / "vus_evidence_action_plan" / "vus_evidence_action_plan_top_candidates.csv"
OUT_DIR = ANALYSIS_DIR / "tables" / "public_classification_snapshot"
SOURCE_DIR = ANALYSIS_DIR / "external_sources" / "public_classification_snapshot_sources"
CACHE_PATH = SOURCE_DIR / "public_classification_snapshot_cache.json"
REPORT = ANALYSIS_DIR / "public_classification_snapshot_report.md"

from backend.lookups.clingen import clingen_erepo_lookup
from backend.lookups.clinvar import clinvar_lookup


PREFERRED_SUBMITTER_PATTERNS = [
    "ENIGMA",
    "Evidence-based Network",
    "ClinGen",
    "Ambry",
    "Invitae",
    "GeneDx",
    "Color",
    "Myriad",
    "Quest",
    "Baylor",
    "Laboratory for Molecular Medicine",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_cache() -> dict[str, Any]:
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    return {}


def save_cache(cache: dict[str, Any]) -> None:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    CACHE_PATH.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except ValueError:
        return default


def select_variants(limit: int) -> list[dict[str, str]]:
    rows = read_csv(INPUT)
    rows = sorted(
        rows,
        key=lambda row: (
            {"tier1_urgent": 0, "tier2_high": 1, "tier3_medium": 2, "tier4_low": 3}.get(row.get("priority_tier", ""), 9),
            -parse_int(row.get("priority_score", "")),
            row.get("gene", ""),
            row.get("c_notation", ""),
        ),
    )
    seen = set()
    selected = []
    for row in rows:
        key = (row["gene"], row["c_notation"])
        if key in seen:
            continue
        seen.add(key)
        selected.append(row)
        if len(selected) >= limit:
            break
    return selected


def snapshot_category(clinvar: dict[str, Any], clingen: dict[str, Any]) -> str:
    if clingen.get("status") == "ok":
        return "ClinGen/ENIGMA assertion"
    if clinvar.get("status") != "ok":
        return "no public assertion"
    aggregate = clinvar.get("aggregate", {})
    review_status = str(aggregate.get("review_status", "")).lower()
    aggregate_class = str(aggregate.get("classification", "")).lower()
    n_submitters = int(aggregate.get("n_submitters", 0) or 0)
    if "conflicting classifications" in review_status or aggregate_class.startswith("conflicting"):
        return "conflicting public assertions"
    if "practice guideline" in review_status or "expert panel" in review_status:
        return "panel-level public assertion"
    if n_submitters > 1:
        return "multi-submitter assertion"
    if n_submitters == 1:
        return "single-submitter assertion"
    return "submitter assertion"


def discordance_label(module1_class: str, public_class: str) -> str:
    public_l = public_class.lower()
    if not public_l or "not provided" in public_l:
        return "no_public_classification"
    if "conflicting classification" in public_l:
        return "public_conflict"
    if module1_class == "3" and any(term in public_l for term in ["pathogenic", "likely pathogenic"]):
        return "module1_vus_public_pathogenic_direction"
    if module1_class == "3" and any(term in public_l for term in ["benign", "likely benign"]):
        return "module1_vus_public_benign_direction"
    if module1_class == "3" and "uncertain" in public_l:
        return "module1_vus_public_vus"
    return "other_or_unclear"


def preferred_submitter_summary(submissions: list[dict[str, Any]]) -> str:
    names = []
    for submission in submissions:
        org = submission.get("org", "")
        if any(pattern.lower() in org.lower() for pattern in PREFERRED_SUBMITTER_PATTERNS):
            klass = submission.get("class", "")
            names.append(f"{org}: {klass}" if klass else org)
    return "; ".join(names[:8])


def query_one(row: dict[str, str], cache: dict[str, Any], refresh_errors: bool) -> dict[str, object]:
    key = f"{row['gene']}:{row['c_notation']}"
    cached = cache.get(key, {})
    clinvar_status = cached.get("clinvar", {}).get("status") if isinstance(cached.get("clinvar"), dict) else None
    clingen_status = cached.get("clingen", {}).get("status") if isinstance(cached.get("clingen"), dict) else None
    needs_refresh = key not in cache or (
        refresh_errors and (clinvar_status == "api_error" or clingen_status == "api_error")
    )
    if needs_refresh:
        clinvar = clinvar_lookup(row["gene"], row["c_notation"])
        time.sleep(0.2)
        clingen = clingen_erepo_lookup(row["gene"], row["c_notation"])
        cache[key] = {
            "clinvar": clinvar,
            "clingen": clingen,
            "queried_at": datetime.now().isoformat(timespec="seconds"),
        }
        save_cache(cache)
    item = cache[key]
    clinvar = item.get("clinvar", {})
    clingen = item.get("clingen", {})
    aggregate = clinvar.get("aggregate", {}) if isinstance(clinvar, dict) else {}
    submissions = clinvar.get("submissions", []) if isinstance(clinvar, dict) else []
    category = snapshot_category(clinvar, clingen)
    public_class = clingen.get("classification") or aggregate.get("classification", "")
    return {
        **row,
        "public_snapshot_category": category,
        "discordance_label": discordance_label(row.get("predicted_class", "3"), public_class),
        "clinvar_status": clinvar.get("status", ""),
        "clinvar_variation_id": clinvar.get("variation_id", ""),
        "clinvar_accession": clinvar.get("accession", ""),
        "clinvar_aggregate_classification": aggregate.get("classification", ""),
        "clinvar_review_status": aggregate.get("review_status", ""),
        "clinvar_n_submitters": aggregate.get("n_submitters", ""),
        "clinvar_has_conflict": "yes" if clinvar.get("has_conflict") else "no",
        "clinvar_unique_classes": ";".join(clinvar.get("unique_classes", [])),
        "preferred_submitters": preferred_submitter_summary(submissions),
        "clingen_status": clingen.get("status", ""),
        "clingen_classification": clingen.get("classification", ""),
        "clingen_caid": clingen.get("caid", ""),
        "clingen_evidence_codes": ";".join(
            f"{code.get('code', '')}:{code.get('status', '')}"
            for code in clingen.get("evidence_codes", [])
        ),
        "public_lookup_queried_at": item.get("queried_at", ""),
    }


def summarize(rows: list[dict[str, object]], key: str) -> list[dict[str, object]]:
    counts = Counter(str(row.get(key, "")) for row in rows)
    return [{key: label, "count": count} for label, count in counts.most_common()]


def cross_tab(rows: list[dict[str, object]], key1: str, key2: str) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], int] = defaultdict(int)
    for row in rows:
        grouped[(str(row.get(key1, "")), str(row.get(key2, "")))] += 1
    return [
        {key1: item[0][0], key2: item[0][1], "count": item[1]}
        for item in sorted(grouped.items(), key=lambda kv: (-kv[1], kv[0]))
    ]


def markdown_table(rows: list[dict[str, object]], columns: list[str], limit: int | None = None) -> str:
    selected = rows if limit is None else rows[:limit]
    if not selected:
        return "| none |\n"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in selected:
        lines.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=80)
    parser.add_argument("--refresh-errors", action="store_true")
    args = parser.parse_args()

    cache = load_cache()
    selected = select_variants(args.limit)
    rows = []
    for idx, row in enumerate(selected, start=1):
        print(f"[{idx}/{len(selected)}] {row['gene']} {row['c_notation']}")
        rows.append(query_one(row, cache, args.refresh_errors))

    fields = [
        "public_snapshot_category",
        "discordance_label",
        "priority_score",
        "priority_tier",
        "action_group",
        "bottleneck_category",
        "gene",
        "c_notation",
        "p_notation",
        "criteria_combo",
        "total_points",
        "spliceai_score",
        "clinvar_status",
        "clinvar_variation_id",
        "clinvar_accession",
        "clinvar_aggregate_classification",
        "clinvar_review_status",
        "clinvar_n_submitters",
        "clinvar_has_conflict",
        "clinvar_unique_classes",
        "preferred_submitters",
        "clingen_status",
        "clingen_classification",
        "clingen_caid",
        "clingen_evidence_codes",
        "public_lookup_queried_at",
    ]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUT_DIR / "public_classification_snapshot_variants.csv", rows, fields)
    by_category = summarize(rows, "public_snapshot_category")
    by_discordance = summarize(rows, "discordance_label")
    by_clinvar_status = summarize(rows, "clinvar_status")
    by_clingen_status = summarize(rows, "clingen_status")
    by_category_action = cross_tab(rows, "public_snapshot_category", "action_group")
    write_csv(OUT_DIR / "public_classification_snapshot_by_category.csv", by_category, ["public_snapshot_category", "count"])
    write_csv(OUT_DIR / "public_classification_snapshot_by_discordance.csv", by_discordance, ["discordance_label", "count"])
    write_csv(OUT_DIR / "public_classification_snapshot_by_clinvar_status.csv", by_clinvar_status, ["clinvar_status", "count"])
    write_csv(OUT_DIR / "public_classification_snapshot_by_clingen_status.csv", by_clingen_status, ["clingen_status", "count"])
    write_csv(
        OUT_DIR / "public_classification_snapshot_by_category_action_group.csv",
        by_category_action,
        ["public_snapshot_category", "action_group", "count"],
    )

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"""# Public Classification Snapshot

Generated: {generated}

## Purpose

This analysis adds a public assertion snapshot to the VUS review worklist. It
does not use public assertions as a validation target. ClinVar, ClinGen, and
ENIGMA-linked assertions are treated as public clinical/laboratory context that can help
prioritize manual review, identify conflicts, and avoid missing already
curated interpretations.

Queried variants: `{len(rows)}` selected from the highest-priority VUS action
plan.

## Snapshot Categories

{markdown_table(by_category, ["public_snapshot_category", "count"])}

## Discordance Worklist Labels

{markdown_table(by_discordance, ["discordance_label", "count"])}

`no public assertion` is a source category, while `no_public_classification` is
a worklist label. These counts can differ. In this snapshot, 7 variants had no
ClinVar record, but one additional variant had a ClinVar record with aggregate
classification `not provided`, giving 8 variants without a usable public
classification label.

## Lookup Status

ClinVar:

{markdown_table(by_clinvar_status, ["clinvar_status", "count"])}

ClinGen Evidence Repository:

{markdown_table(by_clingen_status, ["clingen_status", "count"])}

## ClinGen ERepo BRCA1/2 Coverage Context

ClinGen Evidence Repository is useful when a BRCA1/2 record is present, but it
is not comprehensive for these genes. A manual check of the public ClinGen
ERepo BRCA1/2 table showed the following small record counts:

| Gene | Total ERepo records | Pathogenic | Likely pathogenic | VUS | Likely benign | Benign |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCA1 | 75 | 9 | 1 | 16 | 33 | 16 |
| BRCA2 | 68 | 12 | 6 | 14 | 26 | 10 |

Public ClinGen ERepo entry point: https://erepo.clinicalgenome.org/evrepo/

Therefore, `clingen_status=not_found` should be interpreted as "not found in
this small ERepo subset", not as absence of public interpretation. ClinVar
expert-panel or ENIGMA-linked assertions may still exist for the same variant.

## Category By Evidence Action Group

{markdown_table(by_category_action, ["public_snapshot_category", "action_group", "count"], limit=30)}

## Highest Priority Public Context Rows

{markdown_table(rows, ["priority_score", "priority_tier", "public_snapshot_category", "discordance_label", "gene", "c_notation", "p_notation", "clinvar_aggregate_classification", "clinvar_review_status", "clinvar_n_submitters", "preferred_submitters", "clingen_classification"], limit=30)}

## Interpretation Boundary

This is a snapshot of public assertions, not a benchmark. A public pathogenic
or benign assertion is a reason to inspect the variant, submitter history,
evidence, disease context, and ClinGen/ENIGMA status. It should not be counted
as automatic ACMG/ENIGMA evidence unless the underlying criterion and source
requirements are met.

## Outputs

- `tables/public_classification_snapshot/public_classification_snapshot_variants.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_category.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_discordance.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_clinvar_status.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_clingen_status.csv`
- `tables/public_classification_snapshot/public_classification_snapshot_by_category_action_group.csv`
- `external_sources/public_classification_snapshot_sources/public_classification_snapshot_cache.json`
"""
    REPORT.write_text(text, encoding="utf-8")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
