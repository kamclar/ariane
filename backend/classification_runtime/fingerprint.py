"""Deterministic identity of the active classifier implementation and data."""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path

from backend.config import (
    EXON_CNV_EVIDENCE_MANIFEST_PATH,
    EXON_CNV_EVIDENCE_PATH,
    GENE_POLICY_MANIFEST_PATH,
    GENE_POLICY_METADATA_PATH,
    PS1_PROTEIN_REGISTRY_PATH,
    ST2_SPLICE_EVIDENCE_PATH,
    ST7_PATH,
    TABLE4_PATH,
    TABLE9_PATH,
)
from backend.gene_policy import get_gene_policy
from backend.version import ARIANE_VERSION


PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SMALL_DATA_SOURCES = (
    TABLE4_PATH,
    TABLE9_PATH,
    ST7_PATH,
    PS1_PROTEIN_REGISTRY_PATH,
    ST2_SPLICE_EVIDENCE_PATH,
    EXON_CNV_EVIDENCE_PATH,
    EXON_CNV_EVIDENCE_MANIFEST_PATH,
    GENE_POLICY_MANIFEST_PATH,
    GENE_POLICY_METADATA_PATH,
    PROJECT_ROOT / "backend" / "data" / "clinically_important_residues.json",
    PROJECT_ROOT / "backend" / "data" / "brca_pathogenic_founder_variants.json",
    PROJECT_ROOT / "backend" / "data" / "gnomad" / "gnomad_panel_manifest.json",
    PROJECT_ROOT / "data" / "coordinates" / "coordinate_sources.manifest.json",
    PROJECT_ROOT / "data" / "coordinates" / "coordinate_sources.manifest.metadata.json",
    PROJECT_ROOT / "data" / "coordinates" / "brca_intronic_snv_coordinates.metadata.json",
    PROJECT_ROOT / "data" / "precomputed" / "brca_module1_snv_classification_snapshot.metadata.json",
    PROJECT_ROOT / "data" / "precomputed" / "brca_pp4_clinical_lr_snapshot.metadata.json",
    PROJECT_ROOT / "data" / "precomputed" / "brca_normalized_indel_snapshot.metadata.json",
    PROJECT_ROOT / "data" / "spliceai" / "enigma_v1_2_spliceai_profile.json",
    PROJECT_ROOT / "data" / "reference" / "panel" / "metadata.json",
)
_LARGE_VALIDATED_SOURCES = (
    PROJECT_ROOT / "backend" / "data" / "gnomad" / "gnomad_brca_frequency_snapshot.json",
    PROJECT_ROOT / "backend" / "data" / "gnomad" / "gnomad_brca_coverage_snapshot.json",
)
_CODE_DIRECTORIES = (
    PROJECT_ROOT / "backend" / "classification_dag",
    PROJECT_ROOT / "backend" / "lookups",
    PROJECT_ROOT / "backend" / "modules",
    PROJECT_ROOT / "backend" / "population_frequency",
    PROJECT_ROOT / "backend" / "services",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _small_source_identity(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {"path": str(path.relative_to(PROJECT_ROOT)), "status": "missing"}
    return {
        "path": str(path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "sha256": _sha256_file(path),
    }


def _large_source_identity(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {"path": str(path.relative_to(PROJECT_ROOT)), "status": "missing"}
    stat = path.stat()
    return {
        "path": str(path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _code_identity() -> str:
    digest = hashlib.sha256()
    paths = sorted(
        path
        for directory in _CODE_DIRECTORIES
        for path in directory.rglob("*.py")
        if path.is_file()
    )
    for path in paths:
        digest.update(str(path.relative_to(PROJECT_ROOT)).replace("\\", "/").encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


@lru_cache(maxsize=16)
def classification_fingerprint(gene: str, engine_mode: str) -> str:
    """Return the cache boundary for one active gene policy."""
    configured = get_gene_policy(gene)
    payload = {
        "schema_version": 1,
        "application_version": ARIANE_VERSION,
        "build_revision": os.getenv("ARIANE_BUILD_REVISION", "").strip(),
        "engine_mode": engine_mode,
        "policy_id": configured["policy"]["runtime_policy_id"],
        "policy_version": configured["policy"]["version"],
        "gene": gene,
        "reference_transcript": configured["gene_config"]["reference_transcript"],
        "code_sha256": _code_identity(),
        "data_sources": [_small_source_identity(path) for path in _SMALL_DATA_SOURCES],
        "large_validated_sources": [
            _large_source_identity(path) for path in _LARGE_VALIDATED_SOURCES
        ],
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()
