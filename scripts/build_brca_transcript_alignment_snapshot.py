#!/usr/bin/env python3
"""Build the small BRCA transcript-alignment snapshot used for UTR SNVs.

GRCh37 exon coordinates are extracted from an explicitly supplied official
NCBI GFF. GRCh38 exon coordinates and CDS indices come from the checksum-pinned
cdot panel bundle already shipped with ARIANE. The large upstream GFF is never
copied into the repository.
"""

from __future__ import annotations

import argparse
from datetime import date
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PANEL_DIR = ROOT / "data" / "reference" / "panel"
OUTPUT = ROOT / "data" / "coordinates" / "brca_transcript_alignments.json"
METADATA = (
    ROOT / "data" / "coordinates" / "brca_transcript_alignments.metadata.json"
)
TRANSCRIPTS = {
    "BRCA1": "NM_007294.4",
    "BRCA2": "NM_000059.4",
}
GRCH37_GFF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/"
    "GCF_000001405.25-RS_2024_09/"
    "GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    accession = ""
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            accession = line[1:].split()[0]
            records[accession] = []
        elif accession:
            records[accession].append(line.strip().upper())
    return {key: "".join(value) for key, value in records.items()}


def normalize_exons(
    exons: list[tuple[int, int]],
    *,
    strand: str,
) -> list[dict[str, int]]:
    ordered = sorted(exons, reverse=strand == "-")
    transcript_position = 0
    normalized: list[dict[str, int]] = []
    for genomic_start, genomic_end in ordered:
        length = genomic_end - genomic_start + 1
        normalized.append(
            {
                "transcript_start": transcript_position,
                "transcript_end": transcript_position + length,
                "genomic_start": genomic_start,
                "genomic_end": genomic_end,
            }
        )
        transcript_position += length
    return normalized


def parse_grch37_gff(path: Path) -> dict[str, dict[str, Any]]:
    by_transcript: dict[str, dict[str, Any]] = {
        transcript: {"exons": []} for transcript in TRANSCRIPTS.values()
    }
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or "transcript_id=NM_" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "exon"}:
                continue
            attributes = {
                key: value
                for item in fields[8].split(";")
                if "=" in item
                for key, value in [item.split("=", 1)]
            }
            transcript = attributes.get("transcript_id", "")
            if transcript not in by_transcript:
                continue
            item = by_transcript[transcript]
            if fields[2] == "mRNA":
                item.update(
                    accession=fields[0],
                    chrom=fields[0].removeprefix("NC_0000").split(".", 1)[0].lstrip("0"),
                    strand=fields[6],
                )
            else:
                item["exons"].append((int(fields[3]), int(fields[4])))
    for transcript, item in by_transcript.items():
        if not item.get("accession") or not item["exons"]:
            raise RuntimeError(f"GRCh37 GFF has no complete alignment for {transcript}")
        item["exons"] = normalize_exons(item["exons"], strand=item["strand"])
    return by_transcript


def load_cdot() -> tuple[Path, dict[str, Any]]:
    candidates = list((PANEL_DIR / "cdot").glob("*.json.gz"))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one pinned cdot file, found {len(candidates)}")
    path = candidates[0]
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return path, json.load(handle)


def grch38_alignment(record: dict[str, Any]) -> dict[str, Any]:
    raw = record["genome_builds"]["GRCh38"]
    exons = [
        {
            "transcript_start": int(item[3]) - 1,
            "transcript_end": int(item[4]),
            "genomic_start": int(item[0]) + 1,
            "genomic_end": int(item[1]),
        }
        for item in sorted(raw["exons"], key=lambda exon: int(exon[3]))
    ]
    return {
        "accession": raw["contig"],
        "chrom": raw["contig"].removeprefix("NC_0000").split(".", 1)[0].lstrip("0"),
        "strand": raw["strand"],
        "exons": exons,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grch37-gff", type=Path, required=True)
    parser.add_argument("--version", default=date.today().strftime("%Y.%m.%d"))
    args = parser.parse_args()

    gff = args.grch37_gff.resolve()
    if not gff.is_file():
        raise SystemExit(f"GRCh37 GFF does not exist: {gff}")
    sequences = read_fasta(PANEL_DIR / "fasta" / "transcripts.fa")
    grch37 = parse_grch37_gff(gff)
    cdot_path, cdot = load_cdot()
    panel_metadata = json.loads((PANEL_DIR / "metadata.json").read_text(encoding="utf-8"))

    records: dict[str, Any] = {}
    for gene, transcript in TRANSCRIPTS.items():
        sequence = sequences[transcript]
        cdot_record = cdot["transcripts"][transcript]
        grch37_record = grch37[transcript]
        if grch37_record["exons"][-1]["transcript_end"] != len(sequence):
            raise RuntimeError(f"GRCh37 alignment length mismatch for {transcript}")
        grch38_record = grch38_alignment(cdot_record)
        if grch38_record["exons"][-1]["transcript_end"] != len(sequence):
            raise RuntimeError(f"GRCh38 alignment length mismatch for {transcript}")
        records[gene] = {
            "transcript": transcript,
            "transcript_length": len(sequence),
            "cds_start_index": int(cdot_record["start_codon"]),
            "cds_end_index": int(cdot_record["stop_codon"]),
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "assemblies": {
                "GRCh37": {
                    key: grch37_record[key]
                    for key in ("accession", "chrom", "strand", "exons")
                },
                "GRCh38": grch38_record,
            },
        }

    payload = {
        "schema_version": 1,
        "snapshot_id": "brca_transcript_alignments",
        "snapshot_version": args.version,
        "record_count": len(records) * 2,
        "genes": records,
    }
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    metadata = {
        "schema_version": 1,
        "snapshot_id": payload["snapshot_id"],
        "snapshot_version": payload["snapshot_version"],
        "sha256": sha256(OUTPUT),
        "records": payload["record_count"],
        "validation_status": "approved",
        "sources": {
            "GRCh37": {
                "assembly": "GCF_000001405.25 GRCh37.p13",
                "annotation_release": "GCF_000001405.25-RS_2024_09",
                "url": GRCH37_GFF_URL,
                "sha256": sha256(gff),
            },
            "GRCh38": {
                "assembly": "GRCh38",
                "url": panel_metadata["source"]["cdot_url"],
                "sha256": sha256(cdot_path),
            },
            "transcript_fasta": {
                "path": "data/reference/panel/fasta/transcripts.fa",
                "sha256": sha256(PANEL_DIR / "fasta" / "transcripts.fa"),
            },
        },
    }
    METADATA.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"Wrote {OUTPUT.relative_to(ROOT)} and {METADATA.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
