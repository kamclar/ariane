#!/usr/bin/env python3
"""Build the PDF report for the completed SpliceAI distance experiment."""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import shutil
import sqlite3

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
BASE_DATABASE = ROOT / "data/spliceai/build/distance-stratified/results.sqlite3"
EXTENSION_DATABASE = (
    ROOT / "data/spliceai/build/distance-defaults-server-v2/results.sqlite3"
)
OUTPUT = ROOT / "output/pdf/spliceai_distance_experiment_report.pdf"
WORK = ROOT / "tmp/pdfs/spliceai_distance_experiment"
FONT_REGULAR = Path(r"C:\Windows\Fonts\segoeui.ttf")
FONT_BOLD = Path(r"C:\Windows\Fonts\segoeuib.ttf")
TOLERANCE = 1e-6

PRIMARY = "#2563eb"
PRIMARY_DARK = "#1d4ed8"
INK = "#212529"
MUTED = "#6c757d"
BORDER = "#dee2e6"
PANEL = "#f8f9fa"
ORANGE = "#ea580c"
GRAY = "#94a3b8"
GREEN = "#059669"


def load_results(path: Path) -> dict[str, dict]:
    uri = f"file:{path.resolve().as_posix()}?mode=ro"
    with sqlite3.connect(uri, uri=True) as database:
        rows = database.execute(
            "SELECT case_id, result_json FROM cases WHERE status = 'ok'"
        )
        return {case_id: json.loads(payload) for case_id, payload in rows}


def first_runs(result: dict) -> dict[int, dict]:
    return {
        run["distance"]: run
        for run in result["runs"]
        if run.get("repeat") == 1
    }


def analyse() -> dict:
    base = load_results(BASE_DATABASE)
    extension = load_results(EXTENSION_DATABASE)
    if len(base) != 4047 or set(base) != set(extension):
        raise ValueError("The completed base and extension cohorts do not match")

    combined = {
        case_id: {
            "case": extension[case_id]["case"],
            "runs": {**first_runs(base[case_id]), **first_runs(extension[case_id])},
        }
        for case_id in extension
    }
    pair_definitions = [(50, 500), (500, 4999), (4999, 5000), (5000, 10000)]
    pair_results = []
    changed_rows = []
    transitions = Counter()
    genes = Counter()
    selections = Counter()
    kinds = Counter()
    changed_loci: set[str] = set()

    for smaller, larger in pair_definitions:
        differences = []
        above_tolerance = 0
        band_changes = 0
        for item in combined.values():
            left = item["runs"][smaller]
            right = item["runs"][larger]
            difference = abs(left["maximum_delta"] - right["maximum_delta"])
            differences.append(difference)
            if difference > TOLERANCE:
                above_tolerance += 1
            if left["band_rounded_3dp"] != right["band_rounded_3dp"]:
                band_changes += 1
                if (smaller, larger) == (50, 500):
                    case = item["case"]
                    transition = (
                        left["band_rounded_3dp"], right["band_rounded_3dp"]
                    )
                    transitions[transition] += 1
                    genes[case["gene"]] += 1
                    selections[case["selection"]] += 1
                    kinds[case["variant_kind"]] += 1
                    changed_loci.add(case["locus_id"])
                    changed_rows.append(
                        {
                            "case": case,
                            "d50": left["maximum_delta"],
                            "d500": right["maximum_delta"],
                        }
                    )
        pair_results.append(
            {
                "smaller": smaller,
                "larger": larger,
                "above_tolerance": above_tolerance,
                "band_changes": band_changes,
                "maximum_difference": max(differences),
                "median_difference": sorted(differences)[len(differences) // 2],
            }
        )

    return {
        "combined": combined,
        "pairs": pair_results,
        "changed_rows": changed_rows,
        "transitions": transitions,
        "genes": genes,
        "selections": selections,
        "kinds": kinds,
        "changed_loci": len(changed_loci),
    }


def configure_plot_font() -> None:
    font_manager.fontManager.addfont(str(FONT_REGULAR))
    font_manager.fontManager.addfont(str(FONT_BOLD))
    plt.rcParams.update(
        {
            "font.family": "Segoe UI",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "axes.edgecolor": BORDER,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def save_band_change_chart(analysis: dict, path: Path) -> None:
    pairs = analysis["pairs"]
    labels = [f"D={item['smaller']} vs D={item['larger']}" for item in pairs]
    values = [item["band_changes"] for item in pairs]
    figure, axis = plt.subplots(figsize=(8.8, 2.7))
    bars = axis.barh(
        list(reversed(labels)),
        list(reversed(values)),
        color=list(reversed([ORANGE, PRIMARY, PRIMARY, PRIMARY])),
        height=0.55,
    )
    axis.set_title("Variants crossing a SpliceAI decision threshold", loc="left", pad=12)
    axis.set_xlabel("Number of alleles")
    axis.set_xlim(0, 82)
    axis.grid(axis="x", color="#e5e7eb", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    for bar, value in zip(bars, reversed(values), strict=True):
        axis.text(
            max(value + 1.2, 1.2),
            bar.get_y() + bar.get_height() / 2,
            str(value),
            va="center",
            ha="left",
            color=INK,
            fontweight="bold",
        )
    figure.tight_layout()
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def save_scatter_chart(analysis: dict, path: Path) -> None:
    unchanged_x, unchanged_y, changed_x, changed_y = [], [], [], []
    for item in analysis["combined"].values():
        d50 = item["runs"][50]
        d500 = item["runs"][500]
        changed = d50["band_rounded_3dp"] != d500["band_rounded_3dp"]
        target_x = changed_x if changed else unchanged_x
        target_y = changed_y if changed else unchanged_y
        target_x.append(d50["maximum_delta"])
        target_y.append(d500["maximum_delta"])

    figure, axes = plt.subplots(1, 2, figsize=(9.1, 4.0))
    for axis, limit, title in (
        (axes[0], 1.0, "Full score range"),
        (axes[1], 0.4, "Threshold region"),
    ):
        axis.scatter(
            unchanged_x,
            unchanged_y,
            s=7,
            color=GRAY,
            alpha=0.25,
            linewidths=0,
            label="Same decision band",
        )
        axis.scatter(
            changed_x,
            changed_y,
            s=18,
            color=ORANGE,
            alpha=0.85,
            linewidths=0,
            label="Changed decision band",
        )
        axis.plot([0, limit], [0, limit], color=INK, linewidth=0.8, alpha=0.55)
        for threshold in (0.1, 0.2):
            axis.axhline(threshold, color=PRIMARY, linewidth=0.7, alpha=0.45)
            axis.axvline(threshold, color=PRIMARY, linewidth=0.7, alpha=0.45)
        axis.set_xlim(0, limit)
        axis.set_ylim(0, limit)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(title, loc="left")
        axis.set_xlabel("Maximum delta score at D=50")
        axis.set_ylabel("Maximum delta score at D=500")
        axis.grid(color="#eef2f7", linewidth=0.6)
        axis.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    figure.suptitle(
        "Effect of expanding the scored window from D=50 to D=500",
        x=0.06,
        y=1.01,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0.08, 1, 0.96))
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def band_label(value: str) -> str:
    return {
        "le_0.1": "0.10 or lower",
        "between_0.1_and_0.2": "above 0.10 and below 0.20",
        "ge_0.2": "0.20 or higher",
    }[value]


def save_transition_chart(analysis: dict, path: Path) -> None:
    ordered = [
        ("le_0.1", "between_0.1_and_0.2"),
        ("le_0.1", "ge_0.2"),
        ("between_0.1_and_0.2", "ge_0.2"),
    ]
    labels = [f"{band_label(left)} to {band_label(right)}" for left, right in ordered]
    values = [analysis["transitions"][pair] for pair in ordered]
    figure, axis = plt.subplots(figsize=(8.8, 2.6))
    bars = axis.barh(list(reversed(labels)), list(reversed(values)), color=ORANGE, height=0.55)
    axis.set_title("Threshold transitions observed between D=50 and D=500", loc="left", pad=12)
    axis.set_xlabel("Number of alleles")
    axis.set_xlim(0, 49)
    axis.grid(axis="x", color="#e5e7eb", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    for bar, value in zip(bars, reversed(values), strict=True):
        axis.text(
            value + 0.8,
            bar.get_y() + bar.get_height() / 2,
            str(value),
            va="center",
            ha="left",
            color=INK,
            fontweight="bold",
        )
    figure.tight_layout()
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def register_pdf_fonts() -> None:
    pdfmetrics.registerFont(TTFont("SegoeUI", str(FONT_REGULAR)))
    pdfmetrics.registerFont(TTFont("SegoeUI-Bold", str(FONT_BOLD)))
    pdfmetrics.registerFontFamily(
        "SegoeUI", normal="SegoeUI", bold="SegoeUI-Bold"
    )


def paragraph_styles() -> dict[str, ParagraphStyle]:
    styles = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "ReportTitle",
            parent=styles["Title"],
            fontName="SegoeUI-Bold",
            fontSize=24,
            leading=29,
            textColor=colors.HexColor(INK),
            alignment=TA_LEFT,
            spaceAfter=7 * mm,
        ),
        "subtitle": ParagraphStyle(
            "ReportSubtitle",
            parent=styles["Normal"],
            fontName="SegoeUI",
            fontSize=11,
            leading=16,
            textColor=colors.HexColor(MUTED),
            spaceAfter=8 * mm,
        ),
        "heading": ParagraphStyle(
            "ReportHeading",
            parent=styles["Heading2"],
            fontName="SegoeUI-Bold",
            fontSize=15,
            leading=19,
            textColor=colors.HexColor(INK),
            spaceBefore=2 * mm,
            spaceAfter=4 * mm,
        ),
        "body": ParagraphStyle(
            "ReportBody",
            parent=styles["BodyText"],
            fontName="SegoeUI",
            fontSize=9.2,
            leading=13.4,
            textColor=colors.HexColor(INK),
            spaceAfter=3 * mm,
        ),
        "small": ParagraphStyle(
            "ReportSmall",
            parent=styles["BodyText"],
            fontName="SegoeUI",
            fontSize=7.7,
            leading=10.5,
            textColor=colors.HexColor(MUTED),
        ),
        "card_value": ParagraphStyle(
            "CardValue",
            parent=styles["Normal"],
            fontName="SegoeUI-Bold",
            fontSize=19,
            leading=21,
            textColor=colors.HexColor(PRIMARY_DARK),
            alignment=TA_CENTER,
        ),
        "card_label": ParagraphStyle(
            "CardLabel",
            parent=styles["Normal"],
            fontName="SegoeUI",
            fontSize=7.5,
            leading=10,
            textColor=colors.HexColor(MUTED),
            alignment=TA_CENTER,
        ),
        "table_header": ParagraphStyle(
            "TableHeader",
            parent=styles["Normal"],
            fontName="SegoeUI-Bold",
            fontSize=7.6,
            leading=9.5,
            textColor=colors.white,
            alignment=TA_LEFT,
        ),
        "table_cell": ParagraphStyle(
            "TableCell",
            parent=styles["Normal"],
            fontName="SegoeUI",
            fontSize=7.4,
            leading=9.5,
            textColor=colors.HexColor(INK),
        ),
    }


def page_decoration(canvas, document) -> None:
    width, height = A4
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor(PRIMARY))
    canvas.setLineWidth(2.2)
    canvas.line(18 * mm, height - 14 * mm, width - 18 * mm, height - 14 * mm)
    canvas.setFont("SegoeUI-Bold", 9)
    canvas.setFillColor(colors.HexColor(PRIMARY))
    canvas.drawString(18 * mm, height - 10.5 * mm, "ARIANE")
    canvas.setFont("SegoeUI", 7.5)
    canvas.setFillColor(colors.HexColor(MUTED))
    canvas.drawRightString(width - 18 * mm, 10 * mm, f"Page {document.page}")
    canvas.drawString(18 * mm, 10 * mm, "SpliceAI search-distance experiment")
    canvas.restoreState()


def stat_card(value: str, label: str, styles: dict) -> Table:
    card = Table(
        [
            [Paragraph(value, styles["card_value"])],
            [Paragraph(label, styles["card_label"])],
        ],
        colWidths=[39 * mm],
        rowHeights=[12 * mm, 10 * mm],
    )
    card.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(PANEL)),
                ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(BORDER)),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3 * mm),
            ]
        )
    )
    return card


def format_difference(value: float) -> str:
    if value == 0:
        return "0"
    if value < 0.0001:
        return f"{value:.2e}"
    return f"{value:.6f}"


def build_pdf(analysis: dict, charts: dict[str, Path]) -> None:
    register_pdf_fonts()
    styles = paragraph_styles()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    document = SimpleDocTemplate(
        str(OUTPUT),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=21 * mm,
        bottomMargin=17 * mm,
        title="ARIANE SpliceAI search-distance experiment",
        author="ARIANE",
        subject="Comparison of SpliceAI D values in a stratified BRCA1/2 cohort",
    )

    story = [
        Spacer(1, 7 * mm),
        Paragraph("SpliceAI search-distance experiment", styles["title"]),
        Paragraph(
            "Completed stratified BRCA1/2 run using the model, reference sequence and "
            "transcript annotation pinned by ARIANE. Results dated 17 September 2026.",
            styles["subtitle"],
        ),
    ]
    cards = Table(
        [[
            stat_card("4,047", "alleles completed", styles),
            stat_card("0", "failed evaluations", styles),
            stat_card("75", "band changes at D=50 vs D=500", styles),
            stat_card("0", "band changes at D=500 vs D=4999", styles),
        ]],
        colWidths=[41.5 * mm] * 4,
    )
    cards.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP")]))
    story.extend(
        [
            cards,
            Spacer(1, 7 * mm),
            Paragraph("Result", styles["heading"]),
            Paragraph(
                "D=50 excluded splice signals that were visible at D=500. In 75 alleles "
                "at 52 distinct loci, the maximum score crossed one of the ARIANE "
                "decision thresholds at 0.10 or 0.20. No threshold-band change was "
                "observed between D=500 and D=4999, or among D=4999, D=5000 and "
                "D=10000, in this cohort.",
                styles["body"],
            ),
            Image(str(charts["band_changes"]), width=170 * mm, height=52 * mm),
            Spacer(1, 3 * mm),
            Paragraph(
                "The cohort is an engineering validation set with random and targeted "
                "cases. These counts describe the tested cohort and are not prevalence "
                "estimates.",
                styles["small"],
            ),
            PageBreak(),
            Spacer(1, 5 * mm),
            Paragraph("D=50 compared with D=500", styles["heading"]),
            Paragraph(
                "Each point represents one tested allele. Orange points changed decision "
                "band after the scored window was expanded. Blue guide lines mark the "
                "0.10 and 0.20 thresholds used by the ARIANE profile.",
                styles["body"],
            ),
            Image(str(charts["scatter"]), width=170 * mm, height=75 * mm),
            Spacer(1, 3 * mm),
            Image(str(charts["transitions"]), width=170 * mm, height=50 * mm),
            Spacer(1, 3 * mm),
        ]
    )
    breakdown = Table(
        [[
            Paragraph("52", styles["card_value"]),
            Paragraph("34 / 41", styles["card_value"]),
            Paragraph("23 / 51 / 1", styles["card_value"]),
            Paragraph("74 / 1", styles["card_value"]),
        ], [
            Paragraph("distinct loci", styles["card_label"]),
            Paragraph("BRCA1 / BRCA2 alleles", styles["card_label"]),
            Paragraph("random / targeted / control", styles["card_label"]),
            Paragraph("SNVs / deletions", styles["card_label"]),
        ]],
        colWidths=[41.5 * mm] * 4,
        rowHeights=[10 * mm, 9 * mm],
    )
    breakdown.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(PANEL)),
                ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(BORDER)),
                ("INNERGRID", (0, 0), (-1, -1), 0.4, colors.HexColor(BORDER)),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ]
        )
    )
    story.extend(
        [
            breakdown,
            PageBreak(),
            Spacer(1, 5 * mm),
            Paragraph("Pairwise results", styles["heading"]),
        ]
    )
    table_data = [[
        Paragraph("Comparison", styles["table_header"]),
        Paragraph("Maximum score changed by more than 1e-6", styles["table_header"]),
        Paragraph("Decision-band changes", styles["table_header"]),
        Paragraph("Largest absolute score difference", styles["table_header"]),
    ]]
    for pair in analysis["pairs"]:
        table_data.append(
            [
                Paragraph(
                    f"D={pair['smaller']} vs D={pair['larger']}",
                    styles["table_cell"],
                ),
                Paragraph(str(pair["above_tolerance"]), styles["table_cell"]),
                Paragraph(str(pair["band_changes"]), styles["table_cell"]),
                Paragraph(
                    format_difference(pair["maximum_difference"]),
                    styles["table_cell"],
                ),
            ]
        )
    results_table = Table(
        table_data,
        colWidths=[36 * mm, 57 * mm, 35 * mm, 42 * mm],
        repeatRows=1,
    )
    results_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(PRIMARY_DARK)),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor(PANEL)]),
                ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(BORDER)),
                ("INNERGRID", (0, 0), (-1, -1), 0.35, colors.HexColor(BORDER)),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING", (0, 0), (-1, -1), 2.5 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5 * mm),
            ]
        )
    )
    story.extend(
        [
            results_table,
            Spacer(1, 7 * mm),
            Paragraph("Experimental design", styles["heading"]),
            Paragraph(
                "The run evaluated 4,047 distinct alleles at 1,315 positions across all "
                "98 canonical BRCA1 and BRCA2 exon and intron regions. The set contained "
                "3,919 SNVs, 64 deletions and 64 insertions. It included 882 random SNVs "
                "at 294 loci plus targeted splice-boundary, distance-boundary and validation "
                "cases. Three alternate alleles at the same locus are correlated.",
                styles["body"],
            ),
            Paragraph(
                "Every evaluation used the same pinned SpliceAI 1.3.4 implementation, "
                "container image, GRCh38 reference and GENCODE v49 basic reference "
                "transcripts. The D=4999 results in the extension reproduced the saved "
                "base run. All 4,047 cases completed and no failed evaluation was counted "
                "as agreement.",
                styles["body"],
            ),
            Paragraph("Interpretation", styles["heading"]),
            Paragraph(
                "D=50 is not sufficient for this BRCA1/2 validation set. D=500 and D=4999 "
                "gave the same threshold bands in all tested cases, although 1,106 maximum "
                "scores differed by more than 1e-6 and the largest difference was 0.045821. "
                "The absence of band changes is therefore an empirical result, not proof "
                "that the settings are mathematically equivalent for every possible variant.",
                styles["body"],
            ),
            Paragraph(
                "D=4999, D=5000 and D=10000 produced no threshold-band changes. Their "
                "largest observed maximum-score difference was 7.75e-7, below the prespecified "
                "1e-6 numerical tolerance.",
                styles["body"],
            ),
            Paragraph("Production decision", styles["heading"]),
            Paragraph(
                "No ARIANE configuration or classification rule is changed by this experiment. "
                "The production profile remains pinned to D=10000, as recorded in "
                "data/spliceai/enigma_v1_2_spliceai_profile.json.",
                styles["body"],
            ),
            Spacer(1, 3 * mm),
            Paragraph(
                "Source records: data/spliceai/build/distance-stratified and "
                "data/spliceai/build/distance-defaults-server-v2. Analysis tolerance: 1e-6. "
                "Decision bands: 0.10 or lower, above 0.10 and below 0.20, and 0.20 or higher.",
                styles["small"],
            ),
        ]
    )
    document.build(story, onFirstPage=page_decoration, onLaterPages=page_decoration)


def main() -> None:
    if not FONT_REGULAR.is_file() or not FONT_BOLD.is_file():
        raise FileNotFoundError("Segoe UI font files are unavailable")
    if not BASE_DATABASE.is_file() or not EXTENSION_DATABASE.is_file():
        raise FileNotFoundError("Completed experiment databases are unavailable")
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    configure_plot_font()
    analysis = analyse()
    charts = {
        "band_changes": WORK / "band_changes.png",
        "scatter": WORK / "d50_d500_scatter.png",
        "transitions": WORK / "transitions.png",
    }
    save_band_change_chart(analysis, charts["band_changes"])
    save_scatter_chart(analysis, charts["scatter"])
    save_transition_chart(analysis, charts["transitions"])
    build_pdf(analysis, charts)
    print(OUTPUT)


if __name__ == "__main__":
    main()
