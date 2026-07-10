"""Conceptual figure: position is not variant."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = (
    ROOT
    / "Exploratory_analysis"
    / "precomputed_classification"
    / "plots"
    / "24_position_context_concept"
    / "position_is_not_variant_context_is_not_evidence.svg"
)


def write_svg(path: Path, body: str, width: int, height: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
text {{ font-family: Arial, Helvetica, sans-serif; fill: #172033; }}
.title {{ font-size: 24px; font-weight: 700; }}
.subtitle {{ font-size: 13px; fill: #475569; }}
.panel-title {{ font-size: 16px; font-weight: 700; }}
.label {{ font-size: 13px; }}
.small {{ font-size: 11px; fill: #475569; }}
.tiny {{ font-size: 10px; fill: #475569; }}
.box {{ fill: #ffffff; stroke: #cbd5e1; stroke-width: 1.2; }}
.note {{ fill: #f8fafc; stroke: #cbd5e1; }}
.benign {{ fill: #16a34a; }}
.vus {{ fill: #64748b; }}
.pathogenic {{ fill: #b91c1c; }}
.context {{ fill: #7c3aed; }}
.axis {{ stroke: #334155; stroke-width: 2; }}
.thin {{ stroke: #64748b; stroke-width: 1.2; }}
</style>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def panel(x: int, y: int, w: int, h: int, title: str, letter: str) -> list[str]:
    return [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" class="box"/>',
        f'<circle cx="{x + 26}" cy="{y + 28}" r="14" fill="#e0e7ff" stroke="#818cf8"/>',
        f'<text x="{x + 21}" y="{y + 33}" class="panel-title">{letter}</text>',
        f'<text x="{x + 50}" y="{y + 34}" class="panel-title">{title}</text>',
    ]


def pill(x: float, y: float, text: str, color: str, w: int = 92) -> str:
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="24" rx="12" fill="{color}" opacity="0.88"/>'
        f'<text x="{x + 10}" y="{y + 16}" class="tiny" fill="#ffffff">{text}</text>'
    )


def build() -> None:
    width = 1220
    height = 860
    body: list[str] = [
        '<text x="34" y="38" class="title">Position is not variant: local context is not evidence</text>',
        '<text x="34" y="62" class="subtitle">A precomputed Module 1 map can prioritize review, but local enrichment is not an ACMG/ENIGMA criterion.</text>',
    ]

    px, py, pw, ph = 34, 92, 560, 310
    body += panel(px, py, pw, ph, "One position can produce different variants", "A")
    axis_y = py + 164
    pos_x = px + 220
    body.append(f'<line x1="{px + 70}" y1="{axis_y}" x2="{px + pw - 70}" y2="{axis_y}" class="axis"/>')
    body.append(f'<line x1="{pos_x}" y1="{axis_y - 78}" x2="{pos_x}" y2="{axis_y + 78}" stroke="#7c3aed" stroke-width="3"/>')
    body.append(f'<text x="{pos_x - 42}" y="{axis_y - 92}" class="label">same CDS position</text>')
    variants = [
        ("A>G", "synonymous", "class 1/2", "#16a34a", -92),
        ("A>C", "missense", "class 3", "#64748b", 0),
        ("A>T", "PTC", "class 4/5", "#b91c1c", 92),
    ]
    for ref_alt, consequence, klass, color, offset in variants:
        vx = pos_x + offset
        body.append(f'<line x1="{pos_x}" y1="{axis_y}" x2="{vx}" y2="{axis_y + 58}" class="thin"/>')
        body.append(f'<circle cx="{vx}" cy="{axis_y + 58}" r="8" fill="{color}" stroke="#172033"/>')
        body.append(f'<text x="{vx - 17}" y="{axis_y + 86}" class="small">{ref_alt}</text>')
        body.append(f'<text x="{vx - 35}" y="{axis_y + 104}" class="small">{consequence}</text>')
        body.append(pill(vx - 38, axis_y + 114, klass, color, 86))
    body.append(f'<rect x="{px + 348}" y="{py + 156}" width="172" height="48" rx="8" class="note"/>')
    body.append(f'<text x="{px + 360}" y="{py + 176}" class="small">Exact substitution</text>')
    body.append(f'<text x="{px + 360}" y="{py + 194}" class="small">matters.</text>')

    px, py = 626, 92
    body += panel(px, py, pw, ph, "One codon can contain a mixed class pattern", "B")
    codon_x = px + 92
    codon_y = py + 118
    for i, base in enumerate(["c.1", "c.2", "c.3"]):
        x = codon_x + i * 92
        body.append(f'<rect x="{x}" y="{codon_y}" width="70" height="48" rx="8" fill="#f1f5f9" stroke="#94a3b8"/>')
        body.append(f'<text x="{x + 20}" y="{codon_y + 30}" class="label">{base}</text>')
    points = [
        (codon_x + 10, codon_y + 88, "#16a34a", "benign"),
        (codon_x + 56, codon_y + 116, "#64748b", "VUS"),
        (codon_x + 118, codon_y + 86, "#b91c1c", "pathogenic"),
        (codon_x + 150, codon_y + 124, "#16a34a", "benign"),
        (codon_x + 218, codon_y + 90, "#64748b", "VUS"),
        (codon_x + 244, codon_y + 120, "#b91c1c", "pathogenic"),
    ]
    for x, y, color, label in points:
        body.append(f'<circle cx="{x}" cy="{y}" r="8" fill="{color}" opacity="0.85" stroke="#172033"/>')
        body.append(f'<text x="{x + 12}" y="{y + 4}" class="tiny">{label}</text>')
    body.append(f'<rect x="{px + 34}" y="{py + 270}" width="{pw - 68}" height="24" rx="8" class="note"/>')
    body.append(f'<text x="{px + 48}" y="{py + 286}" class="small">Different SNVs can activate different rules.</text>')

    px, py = 34, 436
    body += panel(px, py, pw, ph, "Local neighborhood is a context signal", "C")
    line_y = py + 166
    body.append(f'<line x1="{px + 70}" y1="{line_y}" x2="{px + pw - 70}" y2="{line_y}" class="axis"/>')
    vus_x = px + pw / 2
    body.append(f'<line x1="{vus_x}" y1="{line_y - 70}" x2="{vus_x}" y2="{line_y + 70}" stroke="#64748b" stroke-width="3"/>')
    body.append(f'<circle cx="{vus_x}" cy="{line_y}" r="13" fill="#64748b" stroke="#172033"/>')
    body.append(f'<text x="{vus_x - 15}" y="{line_y + 38}" class="label">VUS</text>')
    body.append(f'<line x1="{vus_x - 150}" y1="{line_y - 34}" x2="{vus_x - 150}" y2="{line_y + 34}" stroke="#334155"/>')
    body.append(f'<line x1="{vus_x + 150}" y1="{line_y - 34}" x2="{vus_x + 150}" y2="{line_y + 34}" stroke="#334155"/>')
    body.append(f'<text x="{vus_x - 170}" y="{line_y + 62}" class="small">-20 bp</text>')
    body.append(f'<text x="{vus_x + 130}" y="{line_y + 62}" class="small">+20 bp</text>')
    neighbors = [
        (-128, -28, "#16a34a"),
        (-98, 24, "#16a34a"),
        (-58, -22, "#b91c1c"),
        (62, 28, "#16a34a"),
        (96, -24, "#b91c1c"),
        (126, 22, "#16a34a"),
    ]
    for dx, dy, color in neighbors:
        body.append(f'<circle cx="{vus_x + dx}" cy="{line_y + dy}" r="8" fill="{color}" opacity="0.82" stroke="#172033"/>')
    body.append(pill(px + 92, py + 238, "context signal", "#7c3aed", 112))
    body.append(f'<text x="{px + 220}" y="{py + 254}" class="small">Prompt manual review; do not add points.</text>')
    body.append(f'<rect x="{px + 34}" y="{py + 274}" width="{pw - 68}" height="26" rx="8" fill="#fef2f2" stroke="#fecaca"/>')
    body.append(f'<text x="{px + 48}" y="{py + 292}" class="small" fill="#7f1d1d">Not a standalone ACMG/ENIGMA evidence criterion.</text>')

    px, py = 626, 436
    body += panel(px, py, pw, ph, "Example: BRCA1 BRCT1 mixed cluster", "D")
    chart_x = px + 52
    chart_y = py + 86
    chart_w = 410
    chart_h = 140
    body.append(f'<rect x="{chart_x}" y="{chart_y}" width="{chart_w}" height="{chart_h}" fill="#f8fafc" stroke="#cbd5e1"/>')
    bars = [
        ("benign", 139, "#16a34a"),
        ("VUS", 114, "#64748b"),
        ("pathogenic", 47, "#b91c1c"),
    ]
    max_v = 150
    for i, (label, value, color) in enumerate(bars):
        x = chart_x + 42 + i * 120
        h = chart_h * value / max_v
        y = chart_y + chart_h - h
        body.append(f'<rect x="{x}" y="{y}" width="64" height="{h}" fill="{color}" opacity="0.82"/>')
        body.append(f'<text x="{x + 17}" y="{y - 8}" class="small">{value}</text>')
        body.append(f'<text x="{x + 2}" y="{chart_y + chart_h + 20}" class="small">{label}</text>')
    body.append(f'<text x="{chart_x}" y="{chart_y - 16}" class="small">BRCA1 c.5044-c.5143 / p.1682-p.1715</text>')
    body.append(f'<rect x="{px + 34}" y="{py + 252}" width="{pw - 68}" height="46" rx="8" class="note"/>')
    body.append(f'<text x="{px + 52}" y="{py + 271}" class="small">Benign-enriched but mixed: exact variant mechanism still matters.</text>')
    body.append(f'<text x="{px + 52}" y="{py + 289}" class="small">Drivers include BRCT evidence, BS3/PS3, PP3/SpliceAI, and nearby PTCs.</text>')

    body.append('<rect x="34" y="788" width="1152" height="42" rx="12" fill="#eef2ff" stroke="#c7d2fe"/>')
    body.append('<text x="54" y="814" class="label" fill="#312e81">Take-home: positional context can prioritize review, but classification still belongs to accepted variant-level ACMG/ENIGMA evidence.</text>')

    write_svg(OUT, "\n".join(body), width, height)
    print(OUT)


if __name__ == "__main__":
    build()
