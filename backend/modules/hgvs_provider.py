"""Load and validate the checksum-pinned local hgvs/cdot reference provider."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from cdot.hgvs.dataproviders import JSONDataProvider

from backend.config import PANEL_REFERENCE_DIR, TRANSCRIPTS
from backend.modules.panel_seqfetcher import PanelSequenceFetcher


EXPECTED_PACKAGES = {"hgvs": "1.5.7", "cdot": "0.2.30"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise RuntimeError(f"Required reference file cannot be read: {path}: {exc}") from exc
    return digest.hexdigest()


def _load_json(path: Path, label: str) -> dict:
    if not path.is_file():
        raise RuntimeError(f"Required {label} is missing: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Required {label} cannot be loaded: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"Required {label} is not a JSON object: {path}")
    return value


@dataclass(frozen=True)
class PanelProvider:
    data_provider: JSONDataProvider
    seqfetcher: PanelSequenceFetcher
    manifest: dict
    metadata: dict
    transcript_to_protein: dict[str, str]
    gene_to_transcript: dict[str, str]
    provenance: dict[str, str]


@lru_cache(maxsize=1)
def load_panel_provider(reference_dir: str | Path | None = None) -> PanelProvider:
    root = Path(reference_dir) if reference_dir is not None else PANEL_REFERENCE_DIR
    manifest_path = root / "panel_manifest.json"
    metadata_path = root / "metadata.json"
    manifest = _load_json(manifest_path, "panel reference manifest")
    metadata = _load_json(metadata_path, "panel reference metadata")
    if manifest.get("schema_version") != 1 or metadata.get("schema_version") != 1:
        raise RuntimeError("Unsupported panel reference bundle schema")
    if manifest.get("bundle_id") != metadata.get("bundle_id"):
        raise RuntimeError("Panel reference manifest/metadata bundle identity mismatch")

    expected_files = metadata.get("files")
    if not isinstance(expected_files, dict) or not expected_files:
        raise RuntimeError("Panel reference metadata contains no file checksums")
    for relative, expected_hash in expected_files.items():
        path = root / relative
        if not path.is_file():
            raise RuntimeError(f"Required panel reference file is missing: {path}")
        actual_hash = _sha256(path)
        if actual_hash != expected_hash:
            raise RuntimeError(
                f"Panel reference checksum mismatch for {relative}: "
                f"expected {expected_hash}, found {actual_hash}"
            )

    for package, expected_version in EXPECTED_PACKAGES.items():
        actual_version = importlib.metadata.version(package)
        if actual_version != expected_version:
            raise RuntimeError(
                f"Unsupported {package} version: expected {expected_version}, found {actual_version}"
            )

    transcripts = manifest.get("transcripts")
    if not isinstance(transcripts, list) or not transcripts:
        raise RuntimeError("Panel reference manifest contains no transcripts")
    gene_to_transcript: dict[str, str] = {}
    transcript_to_protein: dict[str, str] = {}
    for record in transcripts:
        gene = str(record.get("gene", ""))
        transcript = str(record.get("transcript", ""))
        protein = str(record.get("protein", ""))
        if gene in gene_to_transcript or transcript in transcript_to_protein:
            raise RuntimeError(f"Duplicate panel transcript assignment: {gene}/{transcript}")
        gene_to_transcript[gene] = transcript
        transcript_to_protein[transcript] = protein
    if gene_to_transcript != TRANSCRIPTS:
        raise RuntimeError(
            f"Panel transcript policy mismatch: expected {TRANSCRIPTS}, found {gene_to_transcript}"
        )

    seqfetcher = PanelSequenceFetcher(
        root / "fasta" / "transcripts.fa", root / "fasta" / "proteins.fa"
    )
    expected_accessions = set(transcript_to_protein) | set(transcript_to_protein.values())
    if seqfetcher.accessions() != expected_accessions:
        raise RuntimeError(
            "Panel FASTA accession set differs from manifest: "
            f"expected {sorted(expected_accessions)}, found {sorted(seqfetcher.accessions())}"
        )
    cdot_files = [root / relative for relative in expected_files if relative.startswith("cdot/")]
    if len(cdot_files) != 1:
        raise RuntimeError(f"Expected one panel cdot data file, found {len(cdot_files)}")
    try:
        data_provider = JSONDataProvider([str(cdot_files[0])], seqfetcher=seqfetcher)
    except Exception as exc:
        raise RuntimeError(f"Panel cdot provider cannot be loaded: {exc}") from exc
    for gene, transcript in gene_to_transcript.items():
        available = {item["tx_ac"] for item in data_provider.get_tx_for_gene(gene)}
        if available != {transcript}:
            raise RuntimeError(
                f"Panel cdot provider transcript mismatch for {gene}: {sorted(available)}"
            )

    provenance = {
        "normalization_engine": "biocommons.hgvs",
        "hgvs_version": EXPECTED_PACKAGES["hgvs"],
        "provider": "cdot.JSONDataProvider",
        "cdot_library_version": EXPECTED_PACKAGES["cdot"],
        "cdot_data_release": str(metadata.get("source", {}).get("cdot_release", "")),
        "cdot_source_sha256": str(metadata.get("source", {}).get("cdot_sha256", "")),
        "reference_bundle": str(manifest.get("bundle_id", "")),
        "reference_manifest_sha256": expected_files.get("panel_manifest.json", ""),
    }
    return PanelProvider(
        data_provider=data_provider,
        seqfetcher=seqfetcher,
        manifest=manifest,
        metadata=metadata,
        transcript_to_protein=transcript_to_protein,
        gene_to_transcript=gene_to_transcript,
        provenance=provenance,
    )


def validate_panel_provider() -> None:
    load_panel_provider()
