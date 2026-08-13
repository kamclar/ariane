#!/usr/bin/env python3
"""Build the versioned local transcript-reference bundle used by ARIANE.

The script does not invent a transcript data schema. It selects the approved
panel transcripts from a checksum-pinned cdot release, serializes them with
cdot's own typed model, and verifies each transcript/CDS against independently
downloaded NCBI nucleotide and protein records.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import io
import json
import os
import re
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import msgspec
from cdot import models
from cdot.hgvs.dataproviders import JSONDataProvider


CDOT_VERSION = "0.2.32"
CDOT_SOURCE_URL = (
    "https://github.com/SACGF/cdot/releases/download/data_v0.2.32/"
    "cdot-0.2.32.refseq.GRCh38.json.gz"
)
CDOT_SOURCE_SHA256 = "70fffc92d3178ae0de66efe97acd8835e5db7ae7d9e703b4b03e873ddf27f0f9"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
USER_AGENT = "ARIANE-reference-bundle-builder/1.0 (contact: repository maintainers)"


@dataclass(frozen=True)
class PanelTranscript:
    gene: str
    gene_id: str
    transcript: str
    protein: str
    selection_source: str


PANEL = (
    PanelTranscript("BRCA1", "672", "NM_007294.4", "NP_009225.1", "ENIGMA BRCA1/2 VCEP"),
    PanelTranscript("BRCA2", "675", "NM_000059.4", "NP_000050.3", "ENIGMA BRCA1/2 VCEP"),
)

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def download(url: str, destination: Path) -> Path:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        atomic_write(destination, response.read())
    return destination


def obtain_source(given: str | None, url: str, target: Path) -> Path:
    if given:
        source = Path(given).resolve()
        if not source.is_file():
            raise RuntimeError(f"Source file does not exist: {source}")
        return source
    return download(url, target)


def ncbi_url(database: str, accession: str, rettype: str, retmode: str) -> str:
    return (
        f"{NCBI_EFETCH}?db={database}&id={accession}&rettype={rettype}"
        f"&retmode={retmode}"
    )


def qualifier(feature: ET.Element, name: str) -> list[str]:
    values = []
    for item in feature.findall(".//GBQualifier"):
        if item.findtext("GBQualifier_name") == name:
            values.append(item.findtext("GBQualifier_value") or "")
    return values


def parse_ncbi_transcript(path: Path, expected: PanelTranscript) -> dict:
    record = ET.parse(path).find(".//GBSeq")
    if record is None:
        raise RuntimeError(f"NCBI transcript XML has no GBSeq record: {path}")
    accession = record.findtext("GBSeq_accession-version")
    sequence = (record.findtext("GBSeq_sequence") or "").upper()
    if accession != expected.transcript:
        raise RuntimeError(f"Expected {expected.transcript}, NCBI returned {accession}")
    if not sequence or re.search(r"[^ACGT]", sequence):
        raise RuntimeError(f"Transcript {accession} contains an unsupported nucleotide")
    features = [
        feature for feature in record.findall(".//GBFeature")
        if feature.findtext("GBFeature_key") == "CDS"
        and expected.gene in qualifier(feature, "gene")
        and expected.protein in qualifier(feature, "protein_id")
    ]
    if len(features) != 1:
        raise RuntimeError(f"Expected one matching CDS for {accession}, found {len(features)}")
    feature = features[0]
    location = feature.findtext("GBFeature_location") or ""
    match = re.fullmatch(r"(\d+)\.\.(\d+)", location)
    if not match:
        raise RuntimeError(f"Unsupported CDS location for {accession}: {location!r}")
    cds_start_1, cds_end_1 = map(int, match.groups())
    translation = "".join(qualifier(feature, "translation"))
    if qualifier(feature, "codon_start") != ["1"] or qualifier(feature, "transl_table") != ["1"]:
        raise RuntimeError(f"Unsupported translation settings for {accession}")
    return {
        "accession": accession,
        "sequence": sequence,
        "cds_start_1": cds_start_1,
        "cds_end_1": cds_end_1,
        "protein": expected.protein,
        "annotated_translation": translation,
    }


def parse_fasta(path: Path, expected_accession: str) -> str:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"Invalid FASTA file: {path}")
    returned = lines[0][1:].split()[0]
    if returned != expected_accession:
        raise RuntimeError(f"Expected {expected_accession}, FASTA returned {returned}")
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    if not sequence or re.search(r"[^ABCDEFGHIKLMNPQRSTVWXYZ*]", sequence):
        raise RuntimeError(f"Protein {expected_accession} contains an unsupported residue")
    return sequence.rstrip("*")


def translate_cds(sequence: str, start_1: int, end_1: int) -> str:
    cds = sequence[start_1 - 1:end_1]
    if len(cds) % 3:
        raise RuntimeError("Reference CDS length is not divisible by three")
    translated = "".join(CODON_TABLE[cds[index:index + 3]] for index in range(0, len(cds), 3))
    if not translated.endswith("*"):
        raise RuntimeError("Reference CDS does not end with a stop codon")
    if "*" in translated[:-1]:
        raise RuntimeError("Reference CDS contains an internal stop codon")
    return translated[:-1]


def fasta_bytes(records: Iterable[tuple[str, str]]) -> bytes:
    output: list[str] = []
    for accession, sequence in records:
        output.append(f">{accession}")
        output.extend(sequence[index:index + 70] for index in range(0, len(sequence), 70))
    return ("\n".join(output) + "\n").encode("ascii")


def deterministic_gzip(content: bytes) -> bytes:
    target = io.BytesIO()
    with gzip.GzipFile(fileobj=target, mode="wb", filename="", mtime=0) as handle:
        handle.write(content)
    return target.getvalue()


def build(args: argparse.Namespace) -> None:
    output = Path(args.output_dir).resolve()
    with tempfile.TemporaryDirectory(prefix="ariane-reference-build-") as temporary_name:
        temporary = Path(temporary_name)
        cdot_source = obtain_source(args.cdot_source, CDOT_SOURCE_URL, temporary / "cdot.json.gz")
        actual_source_sha = sha256(cdot_source)
        if actual_source_sha != CDOT_SOURCE_SHA256:
            raise RuntimeError(
                "cdot source checksum mismatch: "
                f"expected {CDOT_SOURCE_SHA256}, found {actual_source_sha}"
            )
        with gzip.open(cdot_source, "rb") as source_handle:
            source_data = models.loads(source_handle.read())
        selected = {}
        selected_genes = {}
        transcript_records = []
        protein_records = []
        manifest_transcripts = []

        for expected in PANEL:
            try:
                cdot_record = source_data.transcripts[expected.transcript]
            except KeyError as exc:
                raise RuntimeError(f"cdot release lacks {expected.transcript}") from exc
            if cdot_record.gene_name != expected.gene or cdot_record.protein != expected.protein:
                raise RuntimeError(f"cdot identity mismatch for {expected.transcript}")
            selected[expected.transcript] = cdot_record
            if expected.gene_id in source_data.genes:
                selected_genes[expected.gene_id] = source_data.genes[expected.gene_id]

            transcript_xml = obtain_source(
                str(Path(args.ncbi_source_dir) / f"{expected.transcript}.xml") if args.ncbi_source_dir else None,
                ncbi_url("nuccore", expected.transcript, "gb", "xml"),
                temporary / f"{expected.transcript}.xml",
            )
            protein_fasta = obtain_source(
                str(Path(args.ncbi_source_dir) / f"{expected.protein}.fa") if args.ncbi_source_dir else None,
                ncbi_url("protein", expected.protein, "fasta", "text"),
                temporary / f"{expected.protein}.fa",
            )
            transcript = parse_ncbi_transcript(transcript_xml, expected)
            protein = parse_fasta(protein_fasta, expected.protein)
            translated = translate_cds(
                transcript["sequence"], transcript["cds_start_1"], transcript["cds_end_1"]
            )
            if translated != protein or translated != transcript["annotated_translation"]:
                raise RuntimeError(f"Independent protein translation mismatch for {expected.transcript}")
            if cdot_record.start_codon != transcript["cds_start_1"] - 1:
                raise RuntimeError(f"cdot/NCBI CDS start mismatch for {expected.transcript}")
            if cdot_record.stop_codon != transcript["cds_end_1"]:
                raise RuntimeError(f"cdot/NCBI CDS end mismatch for {expected.transcript}")

            transcript_records.append((expected.transcript, transcript["sequence"]))
            protein_records.append((expected.protein, protein))
            manifest_transcripts.append({
                "gene": expected.gene,
                "gene_id": expected.gene_id,
                "transcript": expected.transcript,
                "protein": expected.protein,
                "selection_source": expected.selection_source,
                "cds_start_1_based_inclusive": transcript["cds_start_1"],
                "cds_end_1_based_inclusive": transcript["cds_end_1"],
                "transcript_length": len(transcript["sequence"]),
                "protein_length": len(protein),
                "ncbi_transcript_url": ncbi_url(
                    "nuccore", expected.transcript, "gb", "xml"
                ),
                "ncbi_protein_url": ncbi_url(
                    "protein", expected.protein, "fasta", "text"
                ),
                "ncbi_transcript_xml_sha256": sha256(transcript_xml),
                "ncbi_protein_fasta_sha256": sha256(protein_fasta),
            })

        subset = models.CdotData(
            cdot_version=source_data.cdot_version,
            genome_builds=source_data.genome_builds,
            transcripts=selected,
            genes=selected_genes,
            metadata={
                "ariane_bundle_schema": 1,
                "source_url": CDOT_SOURCE_URL,
                "source_sha256": CDOT_SOURCE_SHA256,
                "selected_transcripts": [item.transcript for item in PANEL],
            },
        )
        cdot_name = f"cdot-{CDOT_VERSION}.refseq.GRCh38.brca12.json.gz"
        cdot_content = deterministic_gzip(msgspec.json.encode(subset))
        transcript_content = fasta_bytes(transcript_records)
        protein_content = fasta_bytes(protein_records)
        files = {
            f"cdot/{cdot_name}": cdot_content,
            "fasta/transcripts.fa": transcript_content,
            "fasta/proteins.fa": protein_content,
        }
        for relative, content in files.items():
            atomic_write(output / relative, content)

        manifest = {
            "schema_version": 1,
            "bundle_id": "ariane-brca12-reference-v1",
            "capabilities": ["reference_transcript_c_to_p"],
            "genome_build_for_cdot_alignments": "GRCh38",
            "transcripts": manifest_transcripts,
        }
        atomic_write(
            output / "panel_manifest.json",
            (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )
        file_hashes = {
            relative: sha256(output / relative)
            for relative in [*files, "panel_manifest.json"]
        }
        metadata = {
            "schema_version": 1,
            "bundle_id": manifest["bundle_id"],
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "builder": "scripts/prepare_panel_reference_bundle.py",
            "source": {
                "cdot_release": CDOT_VERSION,
                "cdot_url": CDOT_SOURCE_URL,
                "cdot_sha256": CDOT_SOURCE_SHA256,
                "ncbi_eutils": NCBI_EFETCH,
            },
            "python_packages": {
                name: importlib.metadata.version(name) for name in ("cdot", "hgvs", "msgspec")
            },
            "files": file_hashes,
        }
        atomic_write(
            output / "metadata.json",
            (json.dumps(metadata, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )

        class MissingSequenceFetcher:
            def fetch_seq(self, accession, start_i=None, end_i=None):
                raise RuntimeError(f"Unexpected sequence request while validating {accession}")

        provider = JSONDataProvider([str(output / "cdot" / cdot_name)], seqfetcher=MissingSequenceFetcher())
        if {item["tx_ac"] for item in provider.get_tx_for_gene("BRCA1")} != {"NM_007294.4"}:
            raise RuntimeError("Serialized cdot bundle failed BRCA1 read-back validation")
        if {item["tx_ac"] for item in provider.get_tx_for_gene("BRCA2")} != {"NM_000059.4"}:
            raise RuntimeError("Serialized cdot bundle failed BRCA2 read-back validation")
        print(json.dumps({"output": str(output), "files": file_hashes}, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="data/reference/panel")
    parser.add_argument("--cdot-source", help="Use a previously downloaded checksum-pinned cdot file")
    parser.add_argument(
        "--ncbi-source-dir",
        help="Directory containing NM_*.xml and NP_*.fa files for an offline build",
    )
    return parser.parse_args()


if __name__ == "__main__":
    build(parse_args())
