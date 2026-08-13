#!/usr/bin/env python3
"""Build or verify the BRCA pathogenic-founder exception snapshot.

The seven GeneReviews entries are selected from Table 6 by its explicit
"Founder variant" wording. BRCA1 c.181T>G is added only after the current
ClinVar page confirms both pathogenic/ENIGMA and founder assertions. The
ENIGMA v1.2 page is also checked for the BA1/BS1 exclusion text.
"""

from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import urllib.request


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "backend" / "data" / "brca_pathogenic_founder_variants.json"
CSPEC_URL = "https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2"
GENEREVIEWS_URL = (
    "https://www.ncbi.nlm.nih.gov/books/NBK1247/table/"
    "brca1.T.brca1_and_brca2associated_heredi/?report=objectonly"
)
CLINVAR_C181_URL = "https://www.ncbi.nlm.nih.gov/clinvar/variation/17661/"


class _TableRows(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_row = False
        self.in_cell = False
        self.cell_parts: list[str] = []
        self.row: list[str] = []
        self.rows: list[list[str]] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag == "tr":
            self.in_row = True
            self.row = []
        elif self.in_row and tag in {"td", "th"}:
            self.in_cell = True
            self.cell_parts = []

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.cell_parts.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag in {"td", "th"} and self.in_cell:
            value = re.sub(r"\s+", " ", " ".join(self.cell_parts)).strip()
            value = re.sub(r"[\u200b\u200c\u200d\u2060\ufeff]", "", value)
            self.row.append(value)
            self.in_cell = False
        elif tag == "tr" and self.in_row:
            if self.row:
                self.rows.append(self.row)
            self.in_row = False


def _fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "ARIANE-source-audit/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def _plain_text(html: bytes) -> str:
    parser = _TableRows()
    parser.feed(html.decode("utf-8", errors="replace"))
    return " ".join(" ".join(row) for row in parser.rows)


def _canonical_c(source_c: str) -> str:
    return re.sub(r"(del|dup)[ACGT]+$", r"\1", source_c)


def _protein_notation(value: str) -> str:
    token = re.search(r"p\.([A-Za-z0-9=*?]+)", value)
    if not token:
        raise RuntimeError(f"could not parse protein notation from {value!r}")
    return f"p.({token.group(1)})"


def _extract_genereviews_records(html: bytes) -> list[dict]:
    parser = _TableRows()
    parser.feed(html.decode("utf-8", errors="replace"))
    gene = None
    transcript = None
    records = []
    for row in parser.rows:
        if row and row[0] in {"BRCA1", "BRCA2"}:
            gene = row[0]
            transcript_match = re.search(r"NM_\d+\s*\.\s*\d+", " ".join(row))
            if transcript_match:
                transcript = re.sub(r"\s+", "", transcript_match.group(0))
        if not gene or "Founder variant" not in " ".join(row):
            continue
        row_text = " ".join(row)
        c_match = re.search(
            r"c\.(?:\d+_\d+|\d+)(?:del|dup)[ACGT]*|c\.\d+[ACGT]>[ACGT]",
            row_text,
        )
        p_match = re.search(r"p\.[A-Za-z0-9=*?]+", row_text)
        comment_match = re.search(r"Founder variant[^|]*", row[-1])
        if not c_match or not p_match:
            raise RuntimeError(f"could not parse founder row: {row!r}")
        source_c = c_match.group(0)
        canonical_c = _canonical_c(source_c)
        records.append(
            {
                "gene": gene,
                "transcript": transcript,
                "canonical_c_notation": canonical_c,
                "input_aliases": [source_c] if source_c != canonical_c else [],
                "protein_notation": _protein_notation(p_match.group(0)),
                "founder_context": comment_match.group(0) if comment_match else row[-1],
                "source_assertions": ["NCBI GeneReviews Table 6"],
            }
        )
    if len(records) != 7:
        raise RuntimeError(f"expected 7 GeneReviews founder rows, found {len(records)}")
    return records


def _c181_record(clinvar_html: bytes) -> dict:
    text = _plain_text(clinvar_html)
    required = ("c.181T>G", "Pathogenic", "founder", "ENIGMA")
    missing = [term for term in required if term.lower() not in text.lower()]
    if missing:
        raise RuntimeError(f"ClinVar c.181 source lacks required assertions: {missing}")
    return {
        "gene": "BRCA1",
        "transcript": "NM_007294.4",
        "canonical_c_notation": "c.181T>G",
        "input_aliases": [],
        "protein_notation": "p.(Cys61Gly)",
        "founder_context": "Polish and Eastern European pathogenic founder variant",
        "source_assertions": [
            "ENIGMA expert-panel pathogenic classification in ClinVar VCV000017661",
            "PMID:19594371",
        ],
    }


def _records_sha256(records: list[dict]) -> str:
    canonical = json.dumps(
        records, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def build_payload() -> dict:
    cspec = _fetch(CSPEC_URL)
    genereviews = _fetch(GENEREVIEWS_URL)
    clinvar = _fetch(CLINVAR_C181_URL)
    cspec_text = _plain_text(cspec).lower()
    if "do not apply to well-established pathogenic founder variants" not in cspec_text:
        raise RuntimeError("ENIGMA v1.2 BA1/BS1 founder exception text was not found")

    records = _extract_genereviews_records(genereviews)
    records.append(_c181_record(clinvar))
    records.sort(key=lambda item: (item["gene"], item["canonical_c_notation"]))
    sources = [
        {
            "name": "ENIGMA BRCA1/2 VCEP Criteria Specification v1.2",
            "url": CSPEC_URL,
            "content_sha256_at_access": hashlib.sha256(cspec).hexdigest(),
            "assertion": "BA1 and BS1 must not be applied to well-established pathogenic founder variants",
        },
        {
            "name": "NCBI GeneReviews Table 6: Notable Pathogenic Variants by Gene",
            "url": GENEREVIEWS_URL,
            "content_sha256_at_access": hashlib.sha256(genereviews).hexdigest(),
        },
        {
            "name": "ClinVar VCV000017661",
            "url": CLINVAR_C181_URL,
            "content_sha256_at_access": hashlib.sha256(clinvar).hexdigest(),
            "assertion": "BRCA1 c.181T>G is pathogenic by the ENIGMA expert panel and a documented founder variant",
        },
    ]
    return {
        "metadata": {
            "snapshot_version": "source-verified",
            "enigma_version": "1.2",
            "purpose": "BA1/BS1 pathogenic founder variant exception",
            "records_sha256": _records_sha256(records),
            "selection_policy": (
                "BRCA1/2 variants explicitly described as founder variants in NCBI "
                "GeneReviews Table 6, plus BRCA1 c.181T>G documented as a "
                "pathogenic founder variant by ENIGMA/ClinVar"
            ),
            "sources": sources,
        },
        "variants": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    if args.check:
        current = json.loads(args.output.read_text(encoding="utf-8"))
        def identity(record: dict) -> tuple:
            return (
                record.get("gene"),
                record.get("transcript"),
                record.get("canonical_c_notation"),
                record.get("protein_notation"),
                tuple(sorted(record.get("input_aliases") or [])),
            )

        current_identities = {identity(record) for record in current.get("variants", [])}
        source_identities = {identity(record) for record in payload["variants"]}
        if current_identities != source_identities:
            raise SystemExit("founder snapshot records differ from current official sources")
        print(f"OK: {len(payload['variants'])} founder records verified")
        return 0
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    print(payload["metadata"]["records_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
