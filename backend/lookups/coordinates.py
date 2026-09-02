"""Resolve transcript HGVS variants with validated local coordinate sources.

Coordinate resolution is an evidence-provider concern. It does not assign any
classification criterion. Production requests never call a network coordinate
resolver and never read a mutable coordinate cache. A missing local record is
reported as unavailable so every dependent evidence lookup can fail closed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import logging
from pathlib import Path
import re
from typing import Any, Dict, Optional

from backend.config import TRANSCRIPTS
from backend.data_health import clear_issue, register_issue
from backend.lookups import indels, precomputed


logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[2]
COORDINATE_SOURCE_MANIFEST = (
    ROOT / "data" / "coordinates" / "coordinate_sources.manifest.json"
)
COORDINATE_SOURCE_METADATA = (
    ROOT / "data" / "coordinates" / "coordinate_sources.manifest.metadata.json"
)

_ALLOWED_SOURCE_TYPES = {
    "classification_snapshot",
    "normalized_indel_snapshot",
    "coordinate_map",
}
_SOURCE_MANIFEST: Optional[dict[str, Any]] = None
_COORDINATE_MAPS: dict[str, dict[str, Any]] = {}
_RESOLVER_CACHE: Dict[str, "ResolvedVariant"] = {}


@dataclass
class GenomicCoords:
    chrom: str
    pos: int
    ref: str
    alt: str
    assembly: str

    def variant_id(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"

    def hgvs_g(self) -> str:
        return f"chr{self.chrom}:g.{self.pos}{self.ref}>{self.alt}"

    def is_valid(self) -> bool:
        if not self.chrom or self.pos is None or self.pos < 1 or not self.ref or not self.alt:
            return False
        allowed = {"A", "C", "G", "T"}
        return (
            set(self.ref.upper()).issubset(allowed)
            and set(self.alt.upper()).issubset(allowed)
        )


@dataclass
class ResolvedVariant:
    gene: str
    transcript: str
    c_notation: str
    status: str
    source: str
    grch37: Optional[GenomicCoords] = None
    grch38: Optional[GenomicCoords] = None
    warnings: list[str] = field(default_factory=list)

    def has_grch37(self) -> bool:
        return self.grch37 is not None and self.grch37.is_valid()

    def has_grch38(self) -> bool:
        return self.grch38 is not None and self.grch38.is_valid()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_object(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise RuntimeError(f"Required {label} is missing: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Required {label} cannot be loaded: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"Required {label} must contain a JSON object: {path}")
    return value


def _source_path(relative_path: object, label: str) -> Path:
    if not isinstance(relative_path, str) or not relative_path:
        raise RuntimeError(f"Coordinate source {label} has no path")
    path = (ROOT / relative_path).resolve()
    try:
        path.relative_to(ROOT.resolve())
    except ValueError as exc:
        raise RuntimeError(f"Coordinate source {label} is outside the repository") from exc
    return path


def validate_coordinate_source_manifest() -> dict[str, Any]:
    """Validate the source registry and every directly loaded coordinate map."""
    global _SOURCE_MANIFEST, _COORDINATE_MAPS

    manifest = _read_object(COORDINATE_SOURCE_MANIFEST, "coordinate source manifest")
    metadata = _read_object(COORDINATE_SOURCE_METADATA, "coordinate source metadata")
    expected_sha = str(metadata.get("manifest_sha256") or "").lower()
    actual_sha = _sha256(COORDINATE_SOURCE_MANIFEST)
    if expected_sha != actual_sha:
        raise RuntimeError(
            "Coordinate source manifest checksum mismatch: "
            f"expected {expected_sha}, found {actual_sha}"
        )
    if (
        manifest.get("schema_version") != 1
        or manifest.get("status") != "active"
        or metadata.get("validation_status") != "approved"
        or metadata.get("manifest_id") != manifest.get("manifest_id")
        or metadata.get("manifest_version") != manifest.get("manifest_version")
    ):
        raise RuntimeError("Coordinate source manifest is not an approved schema version 1")
    runtime_policy = manifest.get("runtime_policy")
    if not isinstance(runtime_policy, dict):
        raise RuntimeError("Coordinate source manifest has no runtime policy")
    if runtime_policy.get("network_resolution_allowed") is not False:
        raise RuntimeError("Production coordinate network resolution must be disabled")
    if runtime_policy.get("missing_coordinate_status") != "failed":
        raise RuntimeError("Missing local coordinates must have failed status")

    sources = manifest.get("sources")
    if not isinstance(sources, list) or not sources:
        raise RuntimeError("Coordinate source manifest has no sources")
    seen_ids: set[str] = set()
    seen_priorities: set[int] = set()
    loaded_coordinate_maps: dict[str, dict[str, Any]] = {}
    for source in sources:
        if not isinstance(source, dict):
            raise RuntimeError("Coordinate source entries must be JSON objects")
        source_id = str(source.get("id") or "")
        source_type = source.get("source_type")
        priority = source.get("priority")
        if not source_id or source_id in seen_ids:
            raise RuntimeError(f"Invalid or duplicate coordinate source id: {source_id!r}")
        if source_type not in _ALLOWED_SOURCE_TYPES:
            raise RuntimeError(
                f"Coordinate source {source_id} has unsupported type {source_type!r}"
            )
        if isinstance(priority, bool) or not isinstance(priority, int) or priority in seen_priorities:
            raise RuntimeError(
                f"Coordinate source {source_id} has invalid or duplicate priority"
            )
        seen_ids.add(source_id)
        seen_priorities.add(priority)
        if not str(source.get("runtime_source_label") or "").strip():
            raise RuntimeError(f"Coordinate source {source_id} has no runtime label")
        if source.get("assemblies") != ["GRCh37", "GRCh38"]:
            raise RuntimeError(
                f"Coordinate source {source_id} must declare GRCh37 and GRCh38"
            )
        genes = source.get("genes")
        if not isinstance(genes, dict) or not genes:
            raise RuntimeError(f"Coordinate source {source_id} has no gene bindings")
        for gene, transcript in genes.items():
            if TRANSCRIPTS.get(gene) != transcript:
                raise RuntimeError(
                    f"Coordinate source {source_id} transcript does not match policy for {gene}"
                )
        data_path = _source_path(source.get("path"), source_id)
        metadata_path = _source_path(source.get("metadata_path"), source_id)
        if not data_path.is_file() or not metadata_path.is_file():
            raise RuntimeError(f"Coordinate source {source_id} data or metadata is missing")
        if source_type == "classification_snapshot" and (
            data_path != precomputed.CLASSIFICATION_SNAPSHOT_INDEX.resolve()
            or metadata_path != precomputed.CLASSIFICATION_SNAPSHOT_METADATA.resolve()
        ):
            raise RuntimeError(
                f"Coordinate source {source_id} does not match the configured snapshot loader"
            )
        if source_type == "normalized_indel_snapshot" and (
            data_path != indels.INDEX_PATH.resolve()
            or metadata_path != indels.METADATA_PATH.resolve()
        ):
            raise RuntimeError(
                f"Coordinate source {source_id} does not match the configured snapshot loader"
            )
        if source_type == "coordinate_map":
            source_metadata = _read_object(metadata_path, f"{source_id} metadata")
            if str(source_metadata.get("sha256") or "").lower() != _sha256(data_path):
                raise RuntimeError(f"Coordinate source {source_id} checksum mismatch")
            records = _read_object(data_path, source_id)
            if source_metadata.get("variants") != len(records):
                raise RuntimeError(f"Coordinate source {source_id} record count mismatch")
            for key, record in records.items():
                if not isinstance(record, dict):
                    raise RuntimeError(f"Coordinate source {source_id} has invalid record {key}")
                gene = str(record.get("gene") or "")
                if (
                    key != f"{gene}:{record.get('c_notation')}"
                    or gene not in genes
                    or record.get("transcript") != genes[gene]
                    or record.get("status") != "ok"
                    or _dict_to_coords(record.get("grch37"), "GRCh37") is None
                    or _dict_to_coords(record.get("grch38"), "GRCh38") is None
                ):
                    raise RuntimeError(
                        f"Coordinate source {source_id} has invalid or incomplete record {key}"
                    )
            loaded_coordinate_maps[source_id] = records

    precomputed.validate_classification_snapshot()
    indels.load_indel_snapshot()
    _COORDINATE_MAPS = loaded_coordinate_maps
    _SOURCE_MANIFEST = manifest
    clear_issue("Local coordinate sources")
    return manifest


def _parse_snapshot_coords(value: object, assembly: str) -> Optional[GenomicCoords]:
    if not isinstance(value, str):
        return None
    match = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([ACGT]+)>([ACGT]+)", value)
    if not match:
        return None
    result = GenomicCoords(
        chrom=match.group(1),
        pos=int(match.group(2)),
        ref=match.group(3),
        alt=match.group(4),
        assembly=assembly,
    )
    return result if result.is_valid() else None


def _dict_to_coords(value: object, assembly: str) -> Optional[GenomicCoords]:
    if not isinstance(value, dict):
        return None
    if value.get("assembly") not in (None, assembly):
        return None
    try:
        result = GenomicCoords(
            chrom=str(value["chrom"]),
            pos=int(value["pos"]),
            ref=str(value["ref"]),
            alt=str(value["alt"]),
            assembly=assembly,
        )
    except (KeyError, TypeError, ValueError):
        return None
    return result if result.is_valid() else None


def _resolved(
    *, gene: str, transcript: str, c_notation: str, source: dict[str, Any],
    grch37: Optional[GenomicCoords], grch38: Optional[GenomicCoords],
) -> Optional[ResolvedVariant]:
    if not (grch37 or grch38):
        return None
    return ResolvedVariant(
        gene=gene,
        transcript=transcript,
        c_notation=c_notation,
        status="ok" if grch37 and grch38 else "partial",
        source=str(source["runtime_source_label"]),
        grch37=grch37,
        grch38=grch38,
        warnings=[f"Coordinates loaded from validated local source {source['id']}"],
    )


def _resolve_source(
    source: dict[str, Any], gene: str, c_notation: str,
) -> Optional[ResolvedVariant]:
    transcript = str(source["genes"][gene])
    source_type = source["source_type"]
    if source_type == "normalized_indel_snapshot":
        record = indels.lookup_indel_snapshot(gene, c_notation)
        if not record:
            return None
        return _resolved(
            gene=gene,
            transcript=str(record["reference_transcript"]),
            c_notation=str(record["canonical_c_notation"]),
            source=source,
            grch37=_dict_to_coords(record.get("grch37"), "GRCh37"),
            grch38=_dict_to_coords(record.get("grch38"), "GRCh38"),
        )
    if source_type == "classification_snapshot":
        snapshot = precomputed.lookup_classification_snapshot(gene, c_notation)
        if not snapshot:
            return None
        record = snapshot.get("record", {})
        return _resolved(
            gene=gene,
            transcript=transcript,
            c_notation=c_notation,
            source=source,
            grch37=_parse_snapshot_coords(record.get("grch37"), "GRCh37"),
            grch38=_parse_snapshot_coords(record.get("grch38"), "GRCh38"),
        )
    records = _COORDINATE_MAPS[source["id"]]
    record = records.get(f"{gene}:{c_notation}")
    if not isinstance(record, dict):
        return None
    return _resolved(
        gene=gene,
        transcript=str(record.get("transcript") or transcript),
        c_notation=str(record.get("c_notation") or c_notation),
        source=source,
        grch37=_dict_to_coords(record.get("grch37"), "GRCh37"),
        grch38=_dict_to_coords(record.get("grch38"), "GRCh38"),
    )


def load_local_coordinate_sources() -> None:
    """Load and validate immutable local coordinate sources."""
    try:
        manifest = validate_coordinate_source_manifest()
    except Exception as exc:
        register_issue("Local coordinate sources", str(exc))
        raise
    records = sum(len(value) for value in _COORDINATE_MAPS.values())
    logger.info(
        "Loaded %s local coordinate map records from %s registered sources",
        records,
        len(manifest["sources"]),
    )


def resolve_variant(gene: str, c_notation: str) -> ResolvedVariant:
    """Resolve both assemblies from approved local sources, or fail closed."""
    normalized_gene = gene.strip().upper()
    normalized_c = c_notation.strip()
    key = f"{normalized_gene}:{normalized_c}"
    cached = _RESOLVER_CACHE.get(key)
    if cached is not None:
        return cached
    if _SOURCE_MANIFEST is None:
        load_local_coordinate_sources()
    assert _SOURCE_MANIFEST is not None

    sources = sorted(_SOURCE_MANIFEST["sources"], key=lambda item: item["priority"])
    for source in sources:
        if normalized_gene not in source["genes"]:
            continue
        result = _resolve_source(source, normalized_gene, normalized_c)
        if result is not None:
            _RESOLVER_CACHE[key] = result
            return result

    transcript = TRANSCRIPTS.get(normalized_gene, "")
    return ResolvedVariant(
        gene=normalized_gene,
        transcript=transcript,
        c_notation=normalized_c,
        status="failed",
        source="validated_local_sources",
        warnings=[
            f"No validated local coordinates are available for {normalized_gene} "
            f"{normalized_c}. Coordinate-dependent evidence was not evaluated."
        ],
    )


def resolve_all_variants(
    variants: list[dict[str, Any]], verbose: bool = True,
) -> Dict[str, ResolvedVariant]:
    resolved: Dict[str, ResolvedVariant] = {}
    for variant in variants:
        gene = str(variant["gene"])
        c_notation = str(variant["c_notation"])
        key = f"{gene}:{c_notation}"
        result = resolve_variant(gene, c_notation)
        resolved[key] = result
        if verbose:
            grch37 = (
                f"GRCh37 chr{result.grch37.chrom}:{result.grch37.pos}"
                if result.has_grch37() else "GRCh37 missing"
            )
            grch38 = (
                f"GRCh38 chr{result.grch38.chrom}:{result.grch38.pos}"
                if result.has_grch38() else "GRCh38 missing"
            )
            print(f"  {key}: [{result.status}] via {result.source}")
            print(f"    {grch37}  |  {grch38}")
            for warning in result.warnings:
                print(f"    warning: {warning}")
    return resolved


def get_grch37(
    resolved: Dict[str, ResolvedVariant], gene: str, c_notation: str,
) -> Optional[dict[str, Any]]:
    result = resolved.get(f"{gene}:{c_notation}")
    if result and result.has_grch37():
        return {
            "chrom": result.grch37.chrom,
            "pos": result.grch37.pos,
            "ref": result.grch37.ref,
            "alt": result.grch37.alt,
        }
    return None


def get_grch38(
    resolved: Dict[str, ResolvedVariant], gene: str, c_notation: str,
) -> Optional[dict[str, Any]]:
    result = resolved.get(f"{gene}:{c_notation}")
    if result and result.has_grch38():
        return {
            "chrom": result.grch38.chrom,
            "pos": result.grch38.pos,
            "ref": result.grch38.ref,
            "alt": result.grch38.alt,
        }
    return None


RESOLVED_VARIANTS: Dict[str, ResolvedVariant] = {}

load_local_coordinate_sources()
