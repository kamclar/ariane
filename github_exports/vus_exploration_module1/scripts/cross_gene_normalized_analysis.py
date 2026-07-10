"""Cross-gene normalized summaries for combined exploratory outputs."""

from __future__ import annotations

import csv
import html
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "Exploratory_analysis" / "precomputed_classification"
SNAPSHOT = ROOT / "variant_space_scan" / "outputs" / "brca_module1_full_snv_classification.csv"
VUS_PRIORITY = ANALYSIS / "tables" / "vus_prioritization" / "vus_priority_annotated.csv"
ACTION_PLAN = ANALYSIS / "tables" / "vus_evidence_action_plan" / "vus_evidence_action_plan_variants.csv"
QUEUE_SYNTHESIS = ANALYSIS / "tables" / "priority_queue_synthesis" / "priority_queue_synthesis_variants.csv"
OUT_DIR = ANALYSIS / "tables" / "cross_gene_normalized"
PLOT_DIR = ANALYSIS / "plots" / "31_cross_gene_normalized"
REPORT = ANALYSIS / "cross_gene_normalized_report.md"

GENES = ["BRCA1", "BRCA2"]
CLASS_LABELS = {
    "1": "Benign",
    "2": "Likely Benign",
    "3": "VUS",
    "4": "Likely Pathogenic",
    "5": "Pathogenic",
}
GROUP_LABELS = {
    "benign": "Benign/Likely Benign",
    "vus": "VUS",
    "pathogenic": "Likely Pathogenic/Pathogenic",
}
GROUP_ORDER = ["benign", "vus", "pathogenic"]
QUEUE_COLUMNS = [
    "near_threshold_queue",
    "splice_rna_queue",
    "functional_queue",
    "evidence_conflict_queue",
    "public_assertion_queue",
    "neighborhood_context_queue",
    "low_information_queue",
]
QUEUE_LABELS = {
    "near_threshold_queue": "near threshold",
    "splice_rna_queue": "splice/RNA",
    "functional_queue": "functional",
    "evidence_conflict_queue": "evidence conflict",
    "public_assertion_queue": "public assertion",
    "neighborhood_context_queue": "neighborhood context",
    "low_information_queue": "low information",
}


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


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def class_group(cls: str | object) -> str:
    value = str(cls)
    if value in {"1", "2"}:
        return "benign"
    if value in {"4", "5"}:
        return "pathogenic"
    return "vus"


def parse_criteria(criteria: str) -> set[str]:
    codes = set()
    for item in (criteria or "").split(";"):
        if item:
            codes.add(item.split(":", 1)[0])
    return codes


def rate(count: int, denominator: int) -> str:
    if not denominator:
        return "0.00"
    return f"{count / denominator * 100:.2f}"


def per_1000(count: int, denominator: int) -> str:
    if not denominator:
        return "0.00"
    return f"{count / denominator * 1000:.2f}"


def bar_svg(path: Path, title: str, rows: list[dict[str, object]], category_col: str, value_col: str, subtitle: str) -> None:
    categories = [str(row[category_col]) for row in rows]
    width = 1040
    height = 90 + len(categories) * 48
    left = 230
    top = 70
    plot_w = 700
    max_value = max(float(row[value_col]) for row in rows) if rows else 1.0
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.small{font-size:12px;fill:#43506a}.label{font-size:13px}</style>',
        f'<text x="30" y="30" font-size="20" font-weight="700">{esc(title)}</text>',
        f'<text x="30" y="50" class="small">{esc(subtitle)}</text>',
    ]
    colors = {"BRCA1": "#5578c2", "BRCA2": "#c26a55"}
    for idx, row in enumerate(rows):
        y = top + idx * 48
        value = float(row[value_col])
        gene = str(row.get("gene", ""))
        w = plot_w * value / max_value if max_value else 0
        label = str(row[category_col])
        body.append(f'<text x="{left - 10}" y="{y + 24}" text-anchor="end" class="label">{esc(label)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{w:.1f}" height="30" rx="4" fill="{colors.get(gene, "#5578c2")}"/>')
        body.append(f'<text x="{left + w + 8:.1f}" y="{y + 21}" class="label">{value:.2f}</text>')
    body.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(body), encoding="utf-8")


def paired_gene_rows(
    counts: Counter[tuple[str, str]],
    categories: list[str],
    denominators: dict[str, int],
    denominator_label: str,
) -> list[dict[str, object]]:
    rows = []
    for category in categories:
        for gene in GENES:
            count = counts[(gene, category)]
            rows.append(
                {
                    "gene": gene,
                    "category": category,
                    "count": count,
                    "denominator": denominators.get(gene, 0),
                    "denominator_label": denominator_label,
                    "percent_within_gene": rate(count, denominators.get(gene, 0)),
                    "per_1000": per_1000(count, denominators.get(gene, 0)),
                }
            )
    return rows


def main() -> None:
    snapshot = read_csv(SNAPSHOT)
    priority = read_csv(VUS_PRIORITY)
    action = read_csv(ACTION_PLAN)
    queues = read_csv(QUEUE_SYNTHESIS)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    snapshot_denoms = Counter(row["gene"] for row in snapshot)
    vus_denoms = Counter(row["gene"] for row in snapshot if row.get("predicted_class") == "3")

    class_counts = Counter((row["gene"], row["predicted_class"]) for row in snapshot)
    class_rows = paired_gene_rows(
        class_counts,
        sorted(CLASS_LABELS, key=int),
        dict(snapshot_denoms),
        "all generated coding SNVs in gene",
    )
    for row in class_rows:
        row["label"] = CLASS_LABELS[str(row["category"])]
    write_csv(
        OUT_DIR / "gene_normalized_class_distribution.csv",
        class_rows,
        ["gene", "category", "label", "count", "denominator", "denominator_label", "percent_within_gene", "per_1000"],
    )

    group_counts = Counter((row["gene"], class_group(row["predicted_class"])) for row in snapshot)
    group_rows = paired_gene_rows(group_counts, GROUP_ORDER, dict(snapshot_denoms), "all generated coding SNVs in gene")
    for row in group_rows:
        row["label"] = GROUP_LABELS[str(row["category"])]
    write_csv(
        OUT_DIR / "gene_normalized_grouped_class_distribution.csv",
        group_rows,
        ["gene", "category", "label", "count", "denominator", "denominator_label", "percent_within_gene", "per_1000"],
    )

    criterion_set = sorted({code for row in snapshot for code in parse_criteria(row.get("criteria", ""))})
    criteria_counts = Counter()
    for row in snapshot:
        for code in parse_criteria(row.get("criteria", "")):
            criteria_counts[(row["gene"], code)] += 1
    criteria_rows = paired_gene_rows(criteria_counts, criterion_set, dict(snapshot_denoms), "all generated coding SNVs in gene")
    write_csv(
        OUT_DIR / "gene_normalized_criteria_counts.csv",
        criteria_rows,
        ["gene", "category", "count", "denominator", "denominator_label", "percent_within_gene", "per_1000"],
    )

    tier_categories = sorted({row.get("priority_tier", "") for row in priority if row.get("priority_tier")})
    tier_counts = Counter((row["gene"], row["priority_tier"]) for row in priority)
    tier_rows = paired_gene_rows(tier_counts, tier_categories, dict(vus_denoms), "VUS in gene")
    write_csv(
        OUT_DIR / "gene_normalized_vus_priority_tiers.csv",
        tier_rows,
        ["gene", "category", "count", "denominator", "denominator_label", "percent_within_gene", "per_1000"],
    )

    action_categories = sorted({row.get("action_group", "") for row in action if row.get("action_group")})
    action_counts = Counter((row["gene"], row["action_group"]) for row in action)
    action_rows = paired_gene_rows(action_counts, action_categories, dict(vus_denoms), "VUS in gene")
    write_csv(
        OUT_DIR / "gene_normalized_evidence_action_groups.csv",
        action_rows,
        ["gene", "category", "count", "denominator", "denominator_label", "percent_within_gene", "per_1000"],
    )

    queue_counts = Counter()
    for row in queues:
        for queue in QUEUE_COLUMNS:
            if row.get(queue) == "yes":
                queue_counts[(row["gene"], QUEUE_LABELS[queue])] += 1
    queue_rows = paired_gene_rows(queue_counts, [QUEUE_LABELS[q] for q in QUEUE_COLUMNS], dict(vus_denoms), "VUS in gene")
    write_csv(
        OUT_DIR / "gene_normalized_priority_queues.csv",
        queue_rows,
        ["gene", "category", "count", "denominator", "denominator_label", "percent_within_gene", "per_1000"],
    )

    bar_svg(
        PLOT_DIR / "gene_normalized_grouped_class_per_1000.svg",
        "Grouped generated classes by gene",
        [
            {**row, "plot_label": f"{row['label']} {row['gene']}", "plot_value": float(row["per_1000"])}
            for row in group_rows
        ],
        "plot_label",
        "plot_value",
        "Rate per 1000 generated coding SNVs in each gene.",
    )
    bar_svg(
        PLOT_DIR / "gene_normalized_priority_queues_per_1000_vus.svg",
        "Priority queues by gene",
        [
            {**row, "plot_label": f"{row['category']} {row['gene']}", "plot_value": float(row["per_1000"])}
            for row in queue_rows
        ],
        "plot_label",
        "plot_value",
        "Rate per 1000 VUS in each gene.",
    )
    bar_svg(
        PLOT_DIR / "gene_normalized_evidence_action_groups_per_1000_vus.svg",
        "Evidence action groups by gene",
        [
            {**row, "plot_label": f"{row['category']} {row['gene']}", "plot_value": float(row["per_1000"])}
            for row in action_rows
        ],
        "plot_label",
        "plot_value",
        "Rate per 1000 VUS in each gene.",
    )

    generated = datetime.now().isoformat(timespec="seconds")
    report = f"""# Cross-Gene Normalized Analysis

Generated: {generated}

## Purpose

Several exploratory analyses initially used BRCA1 and BRCA2 together. That is
useful for a first overview, but raw combined counts can be misleading because
BRCA2 contributes more generated coding SNVs than BRCA1 in the current
snapshot.

This analysis adds per-gene normalized tables and plots. For full-snapshot
outputs, rates are normalized to all generated coding SNVs in the gene. For VUS
worklists, rates are normalized to VUS in the gene.

## Denominators

| Gene | All generated coding SNVs | VUS |
| --- | ---: | ---: |
| BRCA1 | {snapshot_denoms['BRCA1']} | {vus_denoms['BRCA1']} |
| BRCA2 | {snapshot_denoms['BRCA2']} | {vus_denoms['BRCA2']} |

## Outputs

- `tables/cross_gene_normalized/gene_normalized_class_distribution.csv`
- `tables/cross_gene_normalized/gene_normalized_grouped_class_distribution.csv`
- `tables/cross_gene_normalized/gene_normalized_criteria_counts.csv`
- `tables/cross_gene_normalized/gene_normalized_vus_priority_tiers.csv`
- `tables/cross_gene_normalized/gene_normalized_evidence_action_groups.csv`
- `tables/cross_gene_normalized/gene_normalized_priority_queues.csv`
- `plots/31_cross_gene_normalized/gene_normalized_grouped_class_per_1000.svg`
- `plots/31_cross_gene_normalized/gene_normalized_priority_queues_per_1000_vus.svg`
- `plots/31_cross_gene_normalized/gene_normalized_evidence_action_groups_per_1000_vus.svg`

## Interpretation Boundary

These normalized rates compare how the generated ARIANE Module 1 map behaves in
BRCA1 versus BRCA2. They do not prove biological differences by themselves.
Some differences reflect gene length, exon structure, available structured
functional data, curated lookup coverage, and how current Module 1 criteria are
implemented.

Findlay SGE-based analyses remain BRCA1-only in the current dataset and should
not be interpreted as a cross-gene comparison.
"""
    REPORT.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
