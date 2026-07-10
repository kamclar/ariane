from __future__ import annotations

import html
import re
import sys
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[2]
ARTICLE = ROOT / "Exploratory_analysis" / "precomputed_classification" / "vus_prioritization_module1_exploratory_report.md"
DEFAULT_OUT = ROOT / "Exploratory_analysis" / "precomputed_classification" / "vus_prioritization_module1_exploratory_report.docx"


FIGURES = [
    (
        "Figure 1A. VUS priority score distribution",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "08_vus_prioritization" / "vus_priority_score_histogram.svg",
    ),
    (
        "Figure 1B. Most common VUS priority reasons",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "08_vus_prioritization" / "vus_priority_reason_counts.svg",
    ),
    (
        "Figure 1C. Neighborhood contribution to VUS priority",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "27_vus_manuscript_critique_response" / "neighborhood_increment_summary.svg",
    ),
    (
        "Figure 1D. Priority tier change after removing neighborhood points",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "27_vus_manuscript_critique_response" / "tier_transition_without_neighborhood.svg",
    ),
    (
        "Figure 2A. VUS pathogenic-region summary",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "18_vus_pathogenic_regions" / "vus_pathogenic_region_summary.svg",
    ),
    (
        "Figure 2B. VUS local neighborhood: pathogenic versus benign, +/-20 bp",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "18_vus_pathogenic_regions" / "vus_pathogenic_vs_benign_neighborhood_20bp.svg",
    ),
    (
        "Figure 2C. BRCA1 BRCT1 mixed VUS cluster, c.5044-c.5143",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "23_brca1_brct1_mixed_cluster" / "brca1_brct1_mixed_cluster_class_profile.svg",
    ),
    (
        "Figure 3. Position is not variant: local context is not evidence",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "24_position_context_concept" / "position_is_not_variant_context_is_not_evidence.svg",
    ),
    (
        "Figure 4A. Regional signal categories",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "25_regional_driver_decomposition" / "regional_signal_categories.svg",
    ),
    (
        "Figure 4B. Regional pathogenic fraction after driver exclusions",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "25_regional_driver_decomposition" / "pathogenic_fraction_driver_ablation.svg",
    ),
    (
        "Figure 4C. BRCA1 regional driver tracks",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "25_regional_driver_decomposition" / "brca1_regional_driver_tracks.svg",
    ),
    (
        "Figure 4D. BRCA2 regional driver tracks",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "25_regional_driver_decomposition" / "brca2_regional_driver_tracks.svg",
    ),
    (
        "Figure 5A. Variant type and generated grouped class distribution",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "26_positional_context_followup" / "variant_type_grouped_class_stacked.svg",
    ),
    (
        "Figure 5B. Splice driver versus regional pathogenic fraction",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "26_positional_context_followup" / "splice_driver_vs_pathogenic_fraction.svg",
    ),
    (
        "Figure 5C. Functional evidence versus regional pathogenic fraction",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "26_positional_context_followup" / "functional_driver_vs_pathogenic_fraction.svg",
    ),
    (
        "Figure 5D. Mixed positions and codons after driver exclusions",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "26_positional_context_followup" / "mixed_position_codon_driver_ablation.svg",
    ),
    (
        "Figure 6. Module 1 VUS evidence bottlenecks",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "12_vus_bottleneck" / "vus_bottleneck_categories.svg",
    ),
    (
        "Figure 7A. Benign versus pathogenic feature rates",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "19_benign_structure_function" / "benign_vs_pathogenic_feature_comparison.svg",
    ),
    (
        "Figure 7B. Benign counterexample tiers",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "20_benign_counterexamples" / "tier_summary.svg",
    ),
    (
        "Figure 7C. BS3 benign variants by domain and mechanism",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "21_bs3_domain_conflicts" / "domain_mechanism_stacked.svg",
    ),
    (
        "Figure 8A. Findlay SGE function score by VUS priority tier",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "28_findlay_sge_vus_tier" / "findlay_sge_score_by_priority_tier.svg",
    ),
    (
        "Figure 8B. VUS evidence action groups",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "29_vus_evidence_action_plan" / "vus_evidence_action_groups.svg",
    ),
    (
        "Figure 9A. Overlapping VUS review queues",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "30_priority_queue_synthesis" / "priority_queue_counts.svg",
    ),
    (
        "Figure 9B. Priority queue overlap",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "30_priority_queue_synthesis" / "priority_queue_overlap.svg",
    ),
    (
        "Figure 10A. Gene-normalized grouped class distribution",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "01_overview" / "gene_normalized_grouped_class_distribution_heatmap.svg",
    ),
    (
        "Figure 10B. Gene-normalized criteria rates",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "01_overview" / "gene_normalized_criteria_rate_heatmap.svg",
    ),
    (
        "Figure 10C. Gene-normalized priority queues per 1000 VUS",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "31_cross_gene_normalized" / "gene_normalized_priority_queues_per_1000_vus.svg",
    ),
    (
        "Figure 10D. Gene-normalized evidence action groups per 1000 VUS",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "31_cross_gene_normalized" / "gene_normalized_evidence_action_groups_per_1000_vus.svg",
    ),
    (
        "Figure 11A. BRCA1 binned SpliceAI and boundary context",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "03_boundary_spliceai" / "brca1_boundary_spliceai_binned_heatmap.svg",
    ),
    (
        "Figure 11B. BRCA2 binned SpliceAI and boundary context",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "03_boundary_spliceai" / "brca2_boundary_spliceai_binned_heatmap.svg",
    ),
    (
        "Figure 11C. Gene-normalized SpliceAI high-rate by boundary and group",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "03_boundary_spliceai" / "gene_normalized_spliceai_high_rate_by_boundary_group_heatmap.svg",
    ),
    (
        "Figure 11D. Gene-normalized boundary distribution by group",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "03_boundary_spliceai" / "gene_normalized_boundary_distribution_by_group_heatmap.svg",
    ),
    (
        "Figure 12A. BRCA1 exon-level Module 1 signals",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "06_exon_vus_conflict" / "brca1_exon_signal_heatmap.svg",
    ),
    (
        "Figure 12B. BRCA2 exon-level Module 1 signals",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "06_exon_vus_conflict" / "brca2_exon_signal_heatmap.svg",
    ),
    (
        "Figure 12C. Conflict sanity-check patterns",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "06_exon_vus_conflict" / "conflict_sanity_check_patterns.svg",
    ),
    (
        "Figure 12D. Highest heuristic VUS score examples",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "06_exon_vus_conflict" / "heuristic_vus_score_examples.svg",
    ),
    (
        "Figure 13A. BRCA1 domains versus BRCA2 domains, within-domain percent",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "07_gene_comparison" / "domain_vs_domain_within_domain_percent_heatmap.svg",
    ),
    (
        "Figure 13B. BRCA1 domains versus BRCA2 domains, normalized density",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "07_gene_comparison" / "domain_vs_domain_normalized_density_heatmap.svg",
    ),
    (
        "Figure 13C. Domain and broad-region Module 1 signal heatmap",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "07_gene_comparison" / "domain_region_signal_heatmap.svg",
    ),
    (
        "Figure 13D. Domain and broad-region grouped class mix",
        ROOT / "Exploratory_analysis" / "precomputed_classification" / "plots" / "07_gene_comparison" / "domain_region_grouped_class_mix.svg",
    ),
]


NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "pic": "http://schemas.openxmlformats.org/drawingml/2006/picture",
}


def clean_inline(text: str) -> str:
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = re.sub(r"\*\*([^*]*)\*\*", r"\1", text)
    text = re.sub(r"\*([^*]*)\*", r"\1", text)
    return text


def text_run(text: str, bold: bool = False, italic: bool = False) -> str:
    props = ""
    if bold or italic:
        inner = ""
        if bold:
            inner += "<w:b/>"
        if italic:
            inner += "<w:i/>"
        props = f"<w:rPr>{inner}</w:rPr>"
    return f"<w:r>{props}<w:t xml:space=\"preserve\">{escape(text)}</w:t></w:r>"


def paragraph(text: str = "", style: str | None = None, align: str | None = None) -> str:
    ppr = ""
    if style or align:
        parts = []
        if style:
            parts.append(f'<w:pStyle w:val="{style}"/>')
        if align:
            parts.append(f'<w:jc w:val="{align}"/>')
        ppr = f"<w:pPr>{''.join(parts)}</w:pPr>"
    return f"<w:p>{ppr}{text_run(clean_inline(text))}</w:p>"


def bullet(text: str) -> str:
    return (
        "<w:p><w:pPr><w:pStyle w:val=\"ListBullet\"/>"
        "<w:ind w:left=\"720\" w:hanging=\"360\"/></w:pPr>"
        f"{text_run(clean_inline(text))}</w:p>"
    )


def numbered(text: str) -> str:
    return (
        "<w:p><w:pPr><w:pStyle w:val=\"ListNumber\"/>"
        "<w:ind w:left=\"720\" w:hanging=\"360\"/></w:pPr>"
        f"{text_run(clean_inline(text))}</w:p>"
    )


def table(rows: list[list[str]]) -> str:
    if not rows:
        return ""
    grid = "".join("<w:gridCol w:w=\"2400\"/>" for _ in rows[0])
    body = []
    for idx, row in enumerate(rows):
        cells = []
        for cell in row:
            shading = "<w:shd w:fill=\"EDEDED\"/>" if idx == 0 else ""
            cells.append(
                "<w:tc>"
                f"<w:tcPr><w:tcW w:w=\"2400\" w:type=\"dxa\"/>{shading}</w:tcPr>"
                f"{paragraph(clean_inline(cell))}"
                "</w:tc>"
            )
        body.append(f"<w:tr>{''.join(cells)}</w:tr>")
    return (
        "<w:tbl>"
        "<w:tblPr><w:tblW w:w=\"0\" w:type=\"auto\"/>"
        "<w:tblBorders><w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>"
        "<w:left w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>"
        "<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>"
        "<w:right w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>"
        "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"D9D9D9\"/>"
        "<w:insideV w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"D9D9D9\"/>"
        "</w:tblBorders></w:tblPr>"
        f"<w:tblGrid>{grid}</w:tblGrid>{''.join(body)}</w:tbl>"
    )


def is_table_separator(line: str) -> bool:
    stripped = line.strip()
    return bool(stripped) and all(ch in "|:- " for ch in stripped) and "---" in stripped


def parse_table(lines: list[str], start: int) -> tuple[str, int]:
    raw = []
    i = start
    while i < len(lines) and lines[i].strip().startswith("|"):
        raw.append(lines[i])
        i += 1
    rows = []
    for idx, line in enumerate(raw):
        if idx == 1 and is_table_separator(line):
            continue
        cells = [html.unescape(cell.strip()) for cell in line.strip().strip("|").split("|")]
        rows.append(cells)
    return table(rows), i


def markdown_to_body(markdown: str) -> list[str]:
    lines = markdown.splitlines()
    body = []
    para = []
    i = 0

    def flush_para() -> None:
        if para:
            body.append(paragraph(" ".join(para)))
            para.clear()

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            flush_para()
            i += 1
            continue
        if stripped.startswith("|"):
            flush_para()
            tbl, i = parse_table(lines, i)
            body.append(tbl)
            continue
        if stripped.startswith("#"):
            flush_para()
            level = len(stripped) - len(stripped.lstrip("#"))
            title = stripped[level:].strip()
            style = "Title" if level == 1 else f"Heading{min(level, 3)}"
            body.append(paragraph(title, style=style))
            i += 1
            continue
        if stripped.startswith("- "):
            flush_para()
            item_parts = [stripped[2:].strip()]
            i += 1
            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.strip()
                if not next_stripped:
                    break
                if (
                    next_stripped.startswith("- ")
                    or next_stripped.startswith("#")
                    or next_stripped.startswith("|")
                    or re.match(r"^\d+\.\s+", next_stripped)
                ):
                    break
                if next_line.startswith((" ", "\t")):
                    item_parts.append(next_stripped)
                    i += 1
                    continue
                break
            body.append(bullet(" ".join(item_parts)))
            continue
        match = re.match(r"^\d+\.\s+(.*)$", stripped)
        if match:
            flush_para()
            item_parts = [match.group(1)]
            i += 1
            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.strip()
                if not next_stripped:
                    break
                if (
                    next_stripped.startswith("- ")
                    or next_stripped.startswith("#")
                    or next_stripped.startswith("|")
                    or re.match(r"^\d+\.\s+", next_stripped)
                ):
                    break
                if next_line.startswith((" ", "\t")):
                    item_parts.append(next_stripped)
                    i += 1
                    continue
                break
            body.append(numbered(" ".join(item_parts)))
            continue
        para.append(stripped)
        i += 1
    flush_para()
    return body


def svg_size(path: Path) -> tuple[int, int]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(r"<svg[^>]*\bwidth=\"([0-9.]+)\"[^>]*\bheight=\"([0-9.]+)\"", text)
    if match:
        return int(float(match.group(1))), int(float(match.group(2)))
    match = re.search(r"viewBox=\"[^\"]*?\s+([0-9.]+)\s+([0-9.]+)\"", text)
    if match:
        return int(float(match.group(1))), int(float(match.group(2)))
    return 900, 500


def image_paragraph(rid: str, path: Path, doc_pr_id: int) -> str:
    width_px, height_px = svg_size(path)
    max_width_emu = int(6.3 * 914400)
    width_emu = width_px * 9525
    height_emu = height_px * 9525
    if width_emu > max_width_emu:
        scale = max_width_emu / width_emu
        width_emu = max_width_emu
        height_emu = int(height_emu * scale)
    name = escape(path.name)
    return f"""
<w:p>
  <w:pPr><w:jc w:val="center"/></w:pPr>
  <w:r>
    <w:drawing>
      <wp:inline distT="0" distB="0" distL="0" distR="0">
        <wp:extent cx="{width_emu}" cy="{height_emu}"/>
        <wp:effectExtent l="0" t="0" r="0" b="0"/>
        <wp:docPr id="{doc_pr_id}" name="{name}"/>
        <wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/></wp:cNvGraphicFramePr>
        <a:graphic>
          <a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
            <pic:pic>
              <pic:nvPicPr>
                <pic:cNvPr id="{doc_pr_id}" name="{name}"/>
                <pic:cNvPicPr/>
              </pic:nvPicPr>
              <pic:blipFill>
                <a:blip r:embed="{rid}"/>
                <a:stretch><a:fillRect/></a:stretch>
              </pic:blipFill>
              <pic:spPr>
                <a:xfrm><a:off x="0" y="0"/><a:ext cx="{width_emu}" cy="{height_emu}"/></a:xfrm>
                <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
              </pic:spPr>
            </pic:pic>
          </a:graphicData>
        </a:graphic>
      </wp:inline>
    </w:drawing>
  </w:r>
</w:p>
"""


def content_types(image_paths: list[Path]) -> str:
    return """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Default Extension="svg" ContentType="image/svg+xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
</Types>
"""


def root_rels() -> str:
    return """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>
"""


def document_rels(image_paths: list[Path]) -> str:
    rels = [
        '<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
    ]
    for idx, path in enumerate(image_paths, start=1):
        rels.append(
            f'<Relationship Id="rIdImg{idx}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{escape(path.name)}"/>'
        )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(rels)
        + "</Relationships>"
    )


def styles_xml() -> str:
    return """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/><w:rPr><w:sz w:val="22"/><w:szCs w:val="22"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:qFormat/><w:rPr><w:b/><w:sz w:val="34"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:qFormat/><w:rPr><w:b/><w:sz w:val="28"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Heading3"><w:name w:val="heading 3"/><w:qFormat/><w:rPr><w:b/><w:sz w:val="24"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="ListBullet"><w:name w:val="List Bullet"/><w:qFormat/></w:style>
  <w:style w:type="paragraph" w:styleId="ListNumber"><w:name w:val="List Number"/><w:qFormat/></w:style>
</w:styles>
"""


def document_xml(body_parts: list[str]) -> str:
    ns = " ".join(f'xmlns:{k}="{v}"' for k, v in NS.items())
    sect = (
        "<w:sectPr>"
        '<w:pgSz w:w="12240" w:h="15840"/>'
        '<w:pgMar w:top="1080" w:right="900" w:bottom="1080" w:left="900" w:header="720" w:footer="720" w:gutter="0"/>'
        "</w:sectPr>"
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f"<w:document {ns}><w:body>{''.join(body_parts)}{sect}</w:body></w:document>"
    )


def build_docx(out: Path = DEFAULT_OUT) -> None:
    markdown = ARTICLE.read_text(encoding="utf-8")
    body = markdown_to_body(markdown)
    image_paths = [path for _, path in FIGURES if path.exists()]
    if image_paths:
        body.append(paragraph("Embedded Figures", style="Heading2"))
    for idx, (caption, path) in enumerate(FIGURES, start=1):
        if not path.exists():
            body.append(paragraph(f"Missing figure: {path}", style="Heading3"))
            continue
        body.append(paragraph(caption, style="Heading3"))
        body.append(image_paragraph(f"rIdImg{len([p for _, p in FIGURES[:idx] if p.exists()])}", path, 100 + idx))

    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types(image_paths))
        zf.writestr("_rels/.rels", root_rels())
        zf.writestr("word/document.xml", document_xml(body))
        zf.writestr("word/styles.xml", styles_xml())
        zf.writestr("word/_rels/document.xml.rels", document_rels(image_paths))
        for path in image_paths:
            zf.write(path, f"word/media/{path.name}")
    print(out)


if __name__ == "__main__":
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT
    build_docx(output)
