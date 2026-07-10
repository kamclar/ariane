"""Build a practical curator review queue from VUS bottleneck analyses."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "Exploratory_analysis" / "precomputed_classification"
DEEP_DIVE_CSV = OUT_DIR / "tables" / "vus_bottleneck_deep_dive" / "deep_dive_action_queue.csv"
SANITY_CSV = OUT_DIR / "tables" / "criteria_sanity_audit" / "sanity_variant_level_audit.csv"
CASE_AUDIT_CSV = OUT_DIR / "tables" / "criteria_sanity_audit" / "pvs1_but_vus_case_audit.csv"
TABLE_DIR = OUT_DIR / "tables" / "curation_review_queue"
REPORT = OUT_DIR / "curation_review_queue.md"


def esc(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def as_float(value: object, default: float = 0.0) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def variant_key(row: dict[str, object]) -> tuple[str, str]:
    return (str(row.get("gene", "")), str(row.get("c_notation", "")))


def tier_for(row: dict[str, str]) -> tuple[str, str, str, str]:
    deep_bucket = row.get("deep_bucket", "")
    category = row.get("bottleneck_category", "")

    if deep_bucket == "PS3_PP3_splice_or_functional_followup":
        return (
            "A",
            "PS3+PP3 one-step-short pathogenic VUS",
            "Functional and computational/splice signals point toward pathogenicity, but ENIGMA combination is still incomplete.",
            "Look for RNA confirmation, PP1, PM3, PS4, PP4, curated PS1, or other accepted independent pathogenic evidence.",
        )
    if deep_bucket in {
        "PS3_PM2_needs_one_more_pathogenic_support",
        "PS3_only_needs_independent_support",
        "complex_strong_pathogenic_combo_review",
    }:
        return (
            "B",
            "Strong pathogenic evidence one step short",
            "The variant has strong pathogenic evidence but lacks the exact additional evidence needed for likely pathogenic.",
            "Search for one additional supporting pathogenic criterion or one moderate pathogenic criterion.",
        )
    if category == "mixed_splice_benign_and_pathogenic":
        return (
            "D",
            "Benign-leaning BP4/BP7 plus PM2 VUS",
            "Benign splice prediction is present, but PM2 background keeps the result in VUS.",
            "Look for RNA/functional benign evidence, BS2, BS4, or other benign support; treat pathogenic-density context cautiously.",
        )
    if category == "PM2_only":
        return (
            "E",
            "PM2-only background VUS",
            "Population absence is the only automated criterion; this is a data gap, not a pathogenic mechanism.",
            "Do not prioritize unless another independent signal or external evidence appears.",
        )
    return (
        "Z",
        "Other",
        "Not part of the primary curator queue.",
        "Manual review only if independently selected.",
    )


def conflict_tier(row: dict[str, str]) -> tuple[str, str, str, str]:
    reasons = row.get("sanity_reasons", "")
    if "PVS1_but_VUS" in reasons:
        return (
            "C",
            "PVS1 plus benign functional conflict",
            "PVS1 is present, but benign functional evidence keeps the variant from automatic upgrade.",
            "Expert adjudication of PVS1 boundary and BS3 functional evidence.",
        )
    if "PP3_with_BS3" in reasons:
        return (
            "C",
            "PP3 plus BS3 conflict",
            "Computational/splice signal and functional benign evidence point in opposite directions.",
            "Review functional assay details and consider RNA or independent clinical evidence.",
        )
    if "pathogenic_and_benign_evidence" in reasons:
        return (
            "C",
            "Mixed benign and pathogenic evidence",
            "The automated evidence contains both benign and pathogenic directions.",
            "Curator evidence-adjudication review before any class movement.",
        )
    return (
        "Z",
        "Other sanity flag",
        "Secondary sanity flag.",
        "Manual review if selected by curator.",
    )


def build_queue() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()

    for row in read_csv(DEEP_DIVE_CSV):
        tier, label, why, next_step = tier_for(row)
        if tier == "Z":
            continue
        priority = priority_rank(row, tier)
        out = {
            "tier": tier,
            "tier_label": label,
            "gene": row["gene"],
            "c_notation": row["c_notation"],
            "p_notation": row["p_notation"],
            "variant_type": row["variant_type"],
            "criteria_combo": row["criteria_combo"],
            "total_points": row["total_points"],
            "predicted_class": "3",
            "priority_score": row.get("priority_score", ""),
            "spliceai_score": row.get("spliceai_score", ""),
            "cds_exon": row.get("cds_exon", ""),
            "boundary_distance": row.get("boundary_distance", ""),
            "near_pathogenic_20bp": row.get("near_pathogenic_20bp", ""),
            "review_reason": why,
            "what_to_check": next_step,
            "source_bucket": row.get("deep_bucket", ""),
            "source_analysis": "vus_bottleneck_deep_dive",
            "sort_priority": priority,
        }
        key = (str(out["tier"]), str(out["gene"]), str(out["c_notation"]))
        if key not in seen:
            rows.append(out)
            seen.add(key)

    for row in read_csv(SANITY_CSV):
        tier, label, why, next_step = conflict_tier(row)
        if tier != "C":
            continue
        out = {
            "tier": tier,
            "tier_label": label,
            "gene": row["gene"],
            "c_notation": row["c_notation"],
            "p_notation": row["p_notation"],
            "variant_type": row["variant_type"],
            "criteria_combo": row["criteria_combo"],
            "total_points": row["total_points"],
            "predicted_class": row["predicted_class"],
            "priority_score": "",
            "spliceai_score": row["spliceai_score"],
            "cds_exon": "",
            "boundary_distance": "",
            "near_pathogenic_20bp": "",
            "review_reason": why,
            "what_to_check": next_step,
            "source_bucket": row["sanity_reasons"],
            "source_analysis": "criteria_sanity_audit",
            "sort_priority": priority_rank(row, tier),
        }
        key = (str(out["tier"]), str(out["gene"]), str(out["c_notation"]))
        if key not in seen:
            rows.append(out)
            seen.add(key)

    return sorted(
        rows,
        key=lambda row: (
            str(row["tier"]),
            -as_int(row["sort_priority"]),
            -as_int(row["priority_score"]),
            -as_float(row["spliceai_score"]),
            str(row["gene"]),
            str(row["c_notation"]),
        ),
    )


def priority_rank(row: dict[str, object], tier: str) -> int:
    if tier == "A":
        return 1000 + as_int(row.get("priority_score")) + round(as_float(row.get("spliceai_score")) * 100)
    if tier == "B":
        return 800 + as_int(row.get("priority_score")) + as_int(row.get("near_pathogenic_20bp"))
    if tier == "C":
        return 900 + round(as_float(row.get("spliceai_score")) * 100)
    if tier == "D":
        return 500 + as_int(row.get("priority_score")) + as_int(row.get("near_pathogenic_20bp"))
    if tier == "E":
        return 100 + as_int(row.get("priority_score"))
    return 0


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_tier: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_tier[(str(row["tier"]), str(row["tier_label"]))].append(row)
    output = []
    for (tier, label), items in sorted(by_tier.items()):
        output.append(
            {
                "tier": tier,
                "tier_label": label,
                "count": len(items),
                "top_source_bucket": Counter(row["source_bucket"] for row in items).most_common(1)[0][0],
                "max_spliceai": f"{max(as_float(row['spliceai_score']) for row in items):.2f}",
            }
        )
    return output


def report_table(rows: list[dict[str, object]], columns: list[str], limit: int = 12) -> str:
    if not rows:
        return "| none |\n| --- |\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows[:limit]:
        body.append("| " + " | ".join(esc(row.get(col, "")) for col in columns) + " |")
    return "\n".join([header, sep] + body)


def tier_subset(rows: list[dict[str, object]], tier: str) -> list[dict[str, object]]:
    return [row for row in rows if row["tier"] == tier]


def write_report(rows: list[dict[str, object]], summary_rows: list[dict[str, object]]) -> None:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    columns = ["tier", "gene", "c_notation", "p_notation", "criteria_combo", "total_points", "spliceai_score", "what_to_check"]
    text = f"""# Curator Review Queue

Generated: {generated}

Inputs:

- `tables/vus_bottleneck_deep_dive/deep_dive_action_queue.csv`
- `tables/criteria_sanity_audit/sanity_variant_level_audit.csv`

## Purpose

This checklist converts the exploratory VUS analyses into a practical curator review queue. It does not change classifications. It selects variants where manual evidence review is most likely to be useful, and separates them from large low-information groups.

## Tier Summary

{report_table(summary_rows, ["tier", "tier_label", "count", "top_source_bucket", "max_spliceai"], limit=20)}

## How To Read The Tiers

| tier | meaning |
| --- | --- |
| A | highest priority pathogenic-side VUS: PS3 and PP3 point in the same direction, but the ENIGMA combination is incomplete |
| B | strong pathogenic evidence one step short, often PS3 plus PM2 |
| C | evidence conflicts, such as PVS1 plus BS3 or PP3 plus BS3 |
| D | benign-leaning splice VUS, usually BP4/BP7 plus PM2 |
| E | PM2-only background VUS, lowest priority unless external evidence appears |

## Tier A: PS3+PP3 One-Step-Short

Why: functional and computational/splice evidence point toward pathogenicity, but another independent accepted evidence item is still needed.

{report_table(tier_subset(rows, "A"), columns, limit=20)}

## Tier B: Strong Pathogenic Evidence One Step Short

Why: these are often close to likely pathogenic, but still need one supporting or moderate pathogenic item under the ENIGMA combination rules.

{report_table(tier_subset(rows, "B"), columns, limit=25)}

## Tier C: Evidence Conflicts

Why: these are not automatic upgrades or downgrades. They need curator adjudication of conflicting evidence directions.

{report_table(tier_subset(rows, "C"), columns, limit=25)}

## Tier D: Benign-Leaning BP4/BP7 Plus PM2

Why: benign splice prediction is present, but PM2 keeps the automated result in VUS. RNA or stronger benign evidence would be most useful.

{report_table(tier_subset(rows, "D"), columns, limit=20)}

## Tier E: PM2-Only Background

Why: PM2 alone is not enough. Most of this tier should stay background unless another independent evidence source appears.

{report_table(tier_subset(rows, "E"), columns, limit=12)}

## Practical Use

Start with tiers A and B for pathogenic-side review. Use tier C as a quality-control and evidence-adjudication queue. Tier D is useful for benign resolution, especially if RNA or additional benign clinical evidence is available. Tier E should mostly stay as background.

## Outputs

- `tables/curation_review_queue/curation_review_queue.csv`
- `tables/curation_review_queue/curation_review_queue_summary.csv`
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> None:
    rows = build_queue()
    summary_rows = summarize(rows)
    fields = [
        "tier",
        "tier_label",
        "gene",
        "c_notation",
        "p_notation",
        "variant_type",
        "criteria_combo",
        "total_points",
        "predicted_class",
        "priority_score",
        "spliceai_score",
        "cds_exon",
        "boundary_distance",
        "near_pathogenic_20bp",
        "review_reason",
        "what_to_check",
        "source_bucket",
        "source_analysis",
    ]
    write_csv(TABLE_DIR / "curation_review_queue.csv", rows, fields)
    write_csv(TABLE_DIR / "curation_review_queue_summary.csv", summary_rows, ["tier", "tier_label", "count", "top_source_bucket", "max_spliceai"])
    write_report(rows, summary_rows)
    print(f"Wrote curator review queue to {REPORT}")


if __name__ == "__main__":
    main()
