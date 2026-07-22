# CoordinateResolver
#
# Translates HGVS c. notation to genomic coordinates in both GRCh37 and GRCh38.
# This is the central first step for all downstream lookups:
#   gnomAD v2.1  -> GRCh37
#   gnomAD v3.1  -> GRCh38
#   myvariant.info /v1/variant -> GRCh37 HGVS genomic ID
#   SpliceAI local BRCA subset -> GRCh38
#
# Resolution order:
#   1. VariantValidator REST API (rest.variantvalidator.org)
#      - designed for clinical use, supports BRCA1/2 explicitly
#      - returns both builds in one call
#      - best error messages, handles edge cases correctly
#   2. Mutalyzer normalize API (mutalyzer.nl) as a reported secondary resolver
#      - academic tool, good for most coding variants
#      - requires separate calls for GRCh37 vs GRCh38
#      - known issue: intronic notation (+/-) needs careful URL encoding
# There is no hand-maintained coordinate fallback.
#      - covers current 13 test variants only
#      - not a production solution
#
# Why not manual calculation?
# BRCA1 is on the minus strand. Converting c. position to genomic coordinate
# requires exact exon-intron boundaries and strand-aware arithmetic.
# Manual calculation is unreliable (previous versions had errors up to 12kb).

from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse
import urllib.error
import logging
from backend.data_health import clear_issue, register_issue
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)
_RESOLVER_FAILURES: Dict[str, list[str]] = {}

VARIANTVALIDATOR_API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"
MUTALYZER_API = "https://mutalyzer.nl/api"

TRANSCRIPTS = {
    "BRCA1": "NM_007294.4",
    "BRCA2": "NM_000059.4",
}

# NC accessions for parsing Mutalyzer output
NC_TO_CHROM = {
    "GRCh37": {"NC_000017.10": "17", "NC_000013.10": "13"},
    "GRCh38": {"NC_000017.11": "17", "NC_000013.11": "13"},
}

@dataclass
class GenomicCoords:
    chrom: str
    pos: int
    ref: str
    alt: str
    assembly: str  # "GRCh37" or "GRCh38"

    def variant_id(self) -> str:
        # format used by gnomAD and SpliceAI: "chrom-pos-ref-alt"
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"

    def hgvs_g(self) -> str:
        # hg19 genomic HGVS used by myvariant.info /v1/variant
        return f"chr{self.chrom}:g.{self.pos}{self.ref}>{self.alt}"

    def is_valid(self) -> bool:
        """
        Validate VCF-style coordinates.

        v1.5.2:
        Earlier versions accidentally accepted only single-base REF/ALT values.
        That made small insertions/deletions look as if they had no coordinates,
        even when VariantValidator returned valid VCF coordinates.

        For SpliceAI and gnomAD we need VCF-style alleles, so multi-base
        REF/ALT strings are valid as long as they contain only A/C/G/T.
        """
        if self.chrom is None or self.pos is None:
            return False
        if not self.ref or not self.alt:
            return False
        allowed = {"A", "C", "G", "T"}
        return set(self.ref.upper()).issubset(allowed) and set(self.alt.upper()).issubset(allowed)


@dataclass
class ResolvedVariant:
    gene: str
    transcript: str
    c_notation: str
    status: str          # "ok" | "partial" | "failed"
    source: str          # which resolver succeeded
    grch37: Optional[GenomicCoords] = None
    grch38: Optional[GenomicCoords] = None
    warnings: list = field(default_factory=list)

    def has_grch37(self) -> bool:
        return self.grch37 is not None and self.grch37.is_valid()

    def has_grch38(self) -> bool:
        return self.grch38 is not None and self.grch38.is_valid()


def _parse_snapshot_coords(value: object, assembly: str) -> Optional[GenomicCoords]:
    """Parse snapshot coordinate strings such as ``17:41201178:G>A``."""
    if not isinstance(value, str):
        return None
    match = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([ACGT]+)>([ACGT]+)", value)
    if not match:
        return None
    return GenomicCoords(
        chrom=match.group(1),
        pos=int(match.group(2)),
        ref=match.group(3),
        alt=match.group(4),
        assembly=assembly,
    )


def _resolve_precomputed_snapshot(gene: str, c_notation: str) -> Optional[ResolvedVariant]:
    """Use versioned coding SNV/indel snapshots before network requests."""
    from backend.lookups.precomputed import lookup_classification_snapshot
    from backend.lookups.indels import lookup_indel_snapshot

    indel = lookup_indel_snapshot(gene, c_notation)
    if indel:
        def convert(value: dict) -> GenomicCoords:
            return GenomicCoords(
                chrom=str(value["chrom"]), pos=int(value["pos"]),
                ref=str(value["ref"]), alt=str(value["alt"]),
                assembly=str(value["assembly"]),
            )
        return ResolvedVariant(
            gene=gene,
            transcript=str(indel["reference_transcript"]),
            c_notation=str(indel["canonical_c_notation"]),
            status="ok",
            source="normalized_indel_snapshot",
            grch37=convert(indel["grch37"]),
            grch38=convert(indel["grch38"]),
            warnings=[
                "Coordinates loaded from the versioned normalized BRCA indel snapshot"
            ],
        )

    snapshot = lookup_classification_snapshot(gene, c_notation)
    if not snapshot:
        return None
    record = snapshot.get("record", {})
    grch37 = _parse_snapshot_coords(record.get("grch37"), "GRCh37")
    grch38 = _parse_snapshot_coords(record.get("grch38"), "GRCh38")
    if not (grch37 or grch38):
        logger.warning(
            "Coordinate snapshot entry has no usable coordinates: key=%s:%s status=%s",
            gene,
            c_notation,
            record.get("coordinate_status"),
        )
        return None
    return ResolvedVariant(
        gene=gene,
        transcript=TRANSCRIPTS.get(gene, ""),
        c_notation=c_notation,
        status="ok" if grch37 and grch38 else "partial",
        source="precomputed_snapshot",
        grch37=grch37,
        grch38=grch38,
        warnings=[
            "Coordinates loaded from versioned precomputed snapshot "
            f"{snapshot.get('snapshot_created', '')} "
            f"(ARIANE {snapshot.get('snapshot_version', '')})"
        ],
    )


# ============================================================
# Resolver 1: VariantValidator
# ============================================================

def _resolve_variantvalidator(transcript: str, c_notation: str) -> Optional[Dict]:
    """
    Call VariantValidator to get genomic coordinates.
    Returns the raw JSON response or None on failure.

    Endpoint: /VariantValidator/variantvalidator/{genome}/{variant}/{select_transcripts}
    We request GRCh37 and GRCh38 in a single call using "all" for select_transcripts.
    """
    hgvs = f"{transcript}:{c_notation}"
    # VariantValidator requires genome_build; we call twice (37 and 38) or use "GRCh37|GRCh38"
    # The API supports comma-separated builds in some versions
    url = f"{VARIANTVALIDATOR_API}/GRCh37/{urllib.parse.quote(hgvs)}/mane_select"
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "BRCA-ACMG/1.7.0"},
    )

    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as exc:
        body = exc.read(1000).decode("utf-8", errors="replace")
        logger.warning(
            "VariantValidator HTTP error: status=%s reason=%s url=%s body=%r",
            exc.code, exc.reason, url, body,
        )
        _RESOLVER_FAILURES.setdefault(hgvs, []).append(
            f"VariantValidator HTTP {exc.code} {exc.reason}: {body[:300]}"
        )
        return None
    except Exception as exc:
        logger.warning(
            "VariantValidator request failed: error_type=%s error=%s url=%s",
            type(exc).__name__, exc, url,
        )
        _RESOLVER_FAILURES.setdefault(hgvs, []).append(
            f"VariantValidator {type(exc).__name__}: {exc}"
        )
        return None


def _parse_variantvalidator(data: Dict, transcript: str) -> tuple:
    """
    Parse VariantValidator response to extract GRCh37 and GRCh38 coords.
    Returns (grch37: Optional[GenomicCoords], grch38: Optional[GenomicCoords]).
    """
    grch37 = grch38 = None

    for key, val in data.items():
        if not isinstance(val, dict):
            continue
        if key in ("flag", "metadata", "validation_warning_1"):
            continue

        # primary_assembly_loci contains GRCh37 and GRCh38 coords
        loci = val.get("primary_assembly_loci", {})

        for assembly_key, locus in loci.items():
            vcf = locus.get("vcf", {})
            if not vcf:
                continue
            chrom = str(vcf.get("chr", "")).replace("chr", "")
            pos   = vcf.get("pos")
            ref   = vcf.get("ref", "")
            alt   = vcf.get("alt", "")

            if not (chrom and pos and ref and alt):
                continue

            coords = GenomicCoords(chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                                   assembly=assembly_key)

            if "GRCh37" in assembly_key or "hg19" in assembly_key:
                grch37 = coords
                grch37.assembly = "GRCh37"
            elif "GRCh38" in assembly_key or "hg38" in assembly_key:
                grch38 = coords
                grch38.assembly = "GRCh38"

    return grch37, grch38


# ============================================================
# Resolver 2: Mutalyzer
# ============================================================

def _resolve_mutalyzer(transcript: str, c_notation: str, assembly: str) -> Optional[GenomicCoords]:
    """
    Call Mutalyzer normalize to get genomic coordinates for one assembly.
    assembly: "GRCh37" or "GRCh38"
    """
    hgvs = f"{transcript}:{c_notation}"
    url = f"{MUTALYZER_API}/normalize/{urllib.parse.quote(hgvs, safe='')}"
    req = urllib.request.Request(url, headers={"User-Agent": "BRCA-ACMG/1.7.0"})

    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
    except urllib.error.HTTPError as e:
        body = e.read(1000).decode("utf-8", errors="replace")
        logger.warning(
            "Mutalyzer HTTP error: status=%s reason=%s assembly=%s url=%s body=%r",
            e.code, e.reason, assembly, url, body,
        )
        _RESOLVER_FAILURES.setdefault(hgvs, []).append(
            f"Mutalyzer {assembly} HTTP {e.code} {e.reason}: {body[:300]}"
        )
        return None
    except Exception as exc:
        logger.warning(
            "Mutalyzer request failed: error_type=%s error=%s assembly=%s url=%s",
            type(exc).__name__, exc, assembly, url,
        )
        _RESOLVER_FAILURES.setdefault(hgvs, []).append(
            f"Mutalyzer {assembly} {type(exc).__name__}: {exc}"
        )
        return None

    nc_map = NC_TO_CHROM[assembly]

    # try equivalent_descriptions first (standard Mutalyzer 3 response)
    for desc in data.get("equivalent_descriptions", {}).get("g", []):
        acc = desc.get("description", "")
        m = re.search(r'g\.(\d+)([ACGT])>([ACGT])', acc)
        if not m:
            continue
        pos, ref, alt = int(m.group(1)), m.group(2), m.group(3)
        chrom = next((c for nc, c in nc_map.items() if nc in acc), None)
        if chrom:
            return GenomicCoords(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=assembly)

    # try chromosomal_descriptions (older Mutalyzer versions)
    for desc in data.get("chromosomal_descriptions", []):
        acc = desc.get("genomic", "") or ""
        m = re.search(r'g\.(\d+)([ACGT])>([ACGT])', acc)
        if not m:
            continue
        pos, ref, alt = int(m.group(1)), m.group(2), m.group(3)
        chrom = next((c for nc, c in nc_map.items() if nc in acc), None)
        if chrom:
            return GenomicCoords(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=assembly)

    logger.warning(
        "Mutalyzer response had no usable %s genomic description: url=%s keys=%s",
        assembly, url, sorted(data.keys()),
    )
    return None


# ============================================================
# Main resolver + persistent cache
# ============================================================

import threading as _threading

_RESOLVER_CACHE: Dict[str, ResolvedVariant] = {}
_COORDS_CACHE_PATH = Path(__file__).resolve().parent.parent / "data" / "coordinates_cache.json"
_INTRONIC_COORDS_PATH = (
    Path(__file__).resolve().parent.parent.parent
    / "data" / "coordinates" / "brca_intronic_snv_coordinates.json"
)
_COORDS_FILE_LOCK  = _threading.Lock()


def _coords_to_dict(c: Optional[GenomicCoords]) -> Optional[dict]:
    if c is None:
        return None
    return {"chrom": c.chrom, "pos": c.pos, "ref": c.ref, "alt": c.alt, "assembly": c.assembly}


def _dict_to_coords(d: Optional[dict]) -> Optional[GenomicCoords]:
    if not d:
        return None
    return GenomicCoords(chrom=d["chrom"], pos=d["pos"], ref=d["ref"],
                         alt=d["alt"], assembly=d["assembly"])


def _load_coords_cache() -> None:
    loaded = 0
    # The generated, versioned map is loaded first. The small read-through cache
    # remains useful for variant classes outside the precomputed SNV spaces.
    for path in (_INTRONIC_COORDS_PATH, _COORDS_CACHE_PATH):
        if not path.exists():
            register_issue(f"Coordinate cache {path.name}", f"cache is missing: {path}")
            continue
        try:
            with open(path, encoding="utf-8") as fh:
                raw = json.load(fh)
            for key, entry in raw.items():
                if key in _RESOLVER_CACHE:
                    continue
                _RESOLVER_CACHE[key] = ResolvedVariant(
                    gene=entry["gene"], transcript=entry["transcript"],
                    c_notation=entry["c_notation"], status=entry["status"],
                    source=entry.get("source", "cache"),
                    grch37=_dict_to_coords(entry.get("grch37")),
                    grch38=_dict_to_coords(entry.get("grch38")),
                )
                loaded += 1
            clear_issue(f"Coordinate cache {path.name}")
        except Exception as exc:
            print(f"Warning: could not load coordinates cache {path}: {exc}")
            register_issue(
                f"Coordinate cache {path.name}",
                f"could not load {path}: {type(exc).__name__}: {exc}",
            )
    print(f"Loaded coordinates cache: {loaded} entries")


def _save_coords_cache() -> None:
    with _COORDS_FILE_LOCK:
        try:
            serialized = {
                k: {
                    "gene":       v.gene,
                    "transcript": v.transcript,
                    "c_notation": v.c_notation,
                    "status":     v.status,
                    "source":     v.source,
                    "grch37":     _coords_to_dict(v.grch37),
                    "grch38":     _coords_to_dict(v.grch38),
                }
                for k, v in _RESOLVER_CACHE.items()
                if v.status in ("ok", "partial")
            }
            with open(_COORDS_CACHE_PATH, "w", encoding="utf-8") as fh:
                json.dump(serialized, fh)
        except Exception as exc:
            print(f"Warning: could not save coordinates cache: {exc}")
            register_issue(
                "Coordinate read-through cache",
                f"could not save {_COORDS_CACHE_PATH}: {type(exc).__name__}: {exc}",
            )


def resolve_variant(gene: str, c_notation: str) -> ResolvedVariant:
    """
    Resolve a variant to genomic coordinates in both GRCh37 and GRCh38.

    Returns a ResolvedVariant with status:
      "ok"      - both builds resolved
      "partial" - at least one build resolved
      "failed"  - no coordinates available

    Results are cached per session.
    """
    key = f"{gene}:{c_notation}"
    if key in _RESOLVER_CACHE:
        return _RESOLVER_CACHE[key]

    transcript = TRANSCRIPTS.get(gene, "")
    result = ResolvedVariant(gene=gene, transcript=transcript, c_notation=c_notation,
                             status="failed", source="none")

    if not transcript:
        result.warnings.append(f"No transcript defined for {gene}")
        _RESOLVER_CACHE[key] = result
        return result
    resolver_key = f"{transcript}:{c_notation}"
    _RESOLVER_FAILURES.pop(resolver_key, None)

    # Stable, versioned coordinates for all coding SNVs are available locally.
    # They are tied to the transcript/build versions recorded in snapshot metadata.
    snapshot_result = _resolve_precomputed_snapshot(gene, c_notation)
    if snapshot_result is not None:
        _RESOLVER_CACHE[key] = snapshot_result
        return snapshot_result

    # --- Attempt 1: VariantValidator ---
    vv_data = _resolve_variantvalidator(transcript, c_notation)
    if vv_data:
        grch37, grch38 = _parse_variantvalidator(vv_data, transcript)
        if grch37 or grch38:
            result.grch37 = grch37
            result.grch38 = grch38
            result.source = "VariantValidator"
            result.status = "ok" if (grch37 and grch38) else "partial"
            if not grch37:
                result.warnings.append("VariantValidator: GRCh37 coords not returned")
            if not grch38:
                result.warnings.append("VariantValidator: GRCh38 coords not returned")
            _RESOLVER_CACHE[key] = result
            _save_coords_cache()
            return result

    result.warnings.append("VariantValidator: no usable response, trying Mutalyzer")

    # --- Attempt 2: Mutalyzer ---
    grch37 = _resolve_mutalyzer(transcript, c_notation, "GRCh37")
    time.sleep(0.1)
    grch38 = _resolve_mutalyzer(transcript, c_notation, "GRCh38")

    if grch37 or grch38:
        result.grch37 = grch37
        result.grch38 = grch38
        result.source = "Mutalyzer"
        result.status = "ok" if (grch37 and grch38) else "partial"
        if not grch37:
            result.warnings.append("Mutalyzer: GRCh37 not resolved")
        if not grch38:
            result.warnings.append("Mutalyzer: GRCh38 not resolved")
        result.warnings.extend(_RESOLVER_FAILURES.pop(resolver_key, []))
        _RESOLVER_CACHE[key] = result
        _save_coords_cache()
        return result

    result.warnings.append("Mutalyzer: no usable response")

    result.warnings.extend(_RESOLVER_FAILURES.pop(resolver_key, []))

    result.warnings.append("All resolvers failed - no coordinates available")
    # Do not cache transient total failures; a later request may succeed.
    return result


def resolve_all_variants(variants: list, verbose: bool = True) -> Dict[str, ResolvedVariant]:
    """
    Resolve coordinates for all variants in the test set.
    Returns dict keyed by 'gene:c_notation'.

    v1.5 note:
    We no longer skip all indels here. PM2 is still not applied to indels,
    but SpliceAI can score small indels if GRCh38 coordinates are available
    and the local SpliceAI indel subset has been prepared. Large exon-level
    events may still fail coordinate resolution and will be handled as
    unavailable evidence.
    """
    resolved = {}
    for v in variants:
        gene, c, vtype = v["gene"], v["c_notation"], v["variant_type"]
        key = f"{gene}:{c}"

        r = resolve_variant(gene, c)
        resolved[key] = r

        if verbose:
            grch37_str = f"GRCh37 chr{r.grch37.chrom}:{r.grch37.pos}" if r.has_grch37() else "GRCh37 missing"
            grch38_str = f"GRCh38 chr{r.grch38.chrom}:{r.grch38.pos}" if r.has_grch38() else "GRCh38 missing"
            print(f"  {key}: [{r.status}] via {r.source}")
            print(f"    {grch37_str}  |  {grch38_str}")
            for w in r.warnings:
                print(f"    warning: {w}")

        time.sleep(0.3)  # rate limiting

    return resolved


# ============================================================
# Convenience accessors used by downstream cells
# ============================================================

def get_grch37(resolved: Dict, gene: str, c_notation: str) -> Optional[Dict]:
    """Return GRCh37 coords as dict {chrom, pos, ref, alt} or None."""
    r = resolved.get(f"{gene}:{c_notation}")
    if r and r.has_grch37():
        return {"chrom": r.grch37.chrom, "pos": r.grch37.pos,
                "ref": r.grch37.ref, "alt": r.grch37.alt}
    return None


def get_grch38(resolved: Dict, gene: str, c_notation: str) -> Optional[Dict]:
    """Return GRCh38 coords as dict {chrom, pos, ref, alt} or None."""
    r = resolved.get(f"{gene}:{c_notation}")
    if r and r.has_grch38():
        return {"chrom": r.grch38.chrom, "pos": r.grch38.pos,
                "ref": r.grch38.ref, "alt": r.grch38.alt}
    return None


RESOLVED_VARIANTS = {}

_load_coords_cache()
