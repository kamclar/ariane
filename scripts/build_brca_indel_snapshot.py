"""Build a versioned normalized BRCA1/2 small-indel snapshot.

Input is the official BRCA Exchange release ``variants_output.tsv``. The
builder is deliberately fail-closed: records without the reference transcript,
canonical coding HGVS, protein consequence, or both genomic mappings are not
included in the runtime index.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data/sources/brca_exchange/output/variants_output.tsv"
DEFAULT_OUTPUT = ROOT / "data/precomputed/brca_normalized_indel_snapshot.index.json"
DEFAULT_METADATA = ROOT / "data/precomputed/brca_normalized_indel_snapshot.metadata.json"
TRANSCRIPTS = {"BRCA1": "NM_007294.4", "BRCA2": "NM_000059.4"}
STRANDS = {"BRCA1": "-", "BRCA2": "+"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
PLACEHOLDERS = {"", "-", ".", "?", "NA", "None"}
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in PLACEHOLDERS else text


def canonical_c(value: str, transcript: str) -> str:
    text = clean(value)
    match = re.fullmatch(r"(?:(NM_\d+\.\d+):)?(c\..+)", text)
    if not match:
        return ""
    stated = match.group(1)
    if stated and stated != transcript:
        return ""
    return match.group(2)


def canonical_p(value: str) -> str:
    text = clean(value)
    if not text:
        return ""
    if ":" in text:
        text = text.rsplit(":", 1)[-1]
    if text.startswith("p.") and not text.startswith("p.("):
        text = f"p.({text[2:]})"
    return re.sub(r"\*(?=\d|\))", "Ter", text)


def variant_type(c_notation: str, p_notation: str) -> str:
    c, p = c_notation.lower(), p_notation.lower()
    unknown_protein = p in {"p.?", "p.(?)"}
    if "delins" in c:
        return "frameshift" if "fs" in p else ("delins" if unknown_protein else "inframe_delins")
    if "dup" in c:
        return "frameshift" if "fs" in p else ("duplication" if unknown_protein else "inframe_duplication")
    if "ins" in c:
        return "frameshift" if "fs" in p else ("insertion" if unknown_protein else "inframe_insertion")
    if "del" in c:
        return "frameshift" if "fs" in p else ("deletion" if unknown_protein else "inframe_deletion")
    return ""


def parse_vcf(value: str, assembly: str) -> dict | None:
    text = clean(value)
    match = re.fullmatch(
        r"(?:chr)?(?P<chrom>\d+)(?::g\.)?:?(?P<pos>\d+)[:](?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)",
        text,
        re.IGNORECASE,
    )
    if not match:
        return None
    return {
        "assembly": assembly,
        "chrom": match.group("chrom"),
        "pos": int(match.group("pos")),
        "ref": match.group("ref").upper(),
        "alt": match.group("alt").upper(),
    }


def grch38_from_row(row: dict[str, str]) -> dict | None:
    chrom, pos = clean(row.get("chr")), clean(row.get("pos"))
    ref, alt = clean(row.get("ref")), clean(row.get("alt"))
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom and pos.isdigit() and re.fullmatch(r"[ACGT]+", ref) and re.fullmatch(r"[ACGT]+", alt):
        return {"assembly": "GRCh38", "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt}
    return parse_vcf(row.get("genomic_vcf38_source", ""), "GRCh38")


def aliases(row: dict[str, str], canonical: str, transcript: str, gene: str, grch38: dict) -> list[str]:
    found = {canonical}
    for field in ("synonyms", "hgvs_cdna_lovd", "hgvs_cdna_gnomadv2", "hgvs_cdna_gnomadv3"):
        value = clean(row.get(field))
        for token in re.split(r"[,;|\s]+", value):
            if token.startswith("c."):
                item = token
            elif token.startswith(transcript + ":c."):
                item = token
            elif token.startswith(transcript + ".c."):
                item = token[len(transcript) + 1:]
            else:
                continue
            normalized = canonical_c(item, transcript)
            if normalized:
                found.add(normalized)
    # A single-base duplication is also commonly submitted with the duplicated
    # base spelled out (c.5266dupC). Derive it from the normalized VCF allele.
    if re.fullmatch(r"c\.\d+dup", canonical):
        ref, alt = grch38["ref"], grch38["alt"]
        if alt.startswith(ref) and len(alt) == len(ref) + 1:
            inserted = alt[-1]
            if STRANDS[gene] == "-":
                inserted = inserted.translate(COMPLEMENT)
            found.add(canonical + inserted)
    return sorted(found)


def build(source: Path, output: Path, metadata_path: Path, release: str, release_date: str) -> dict:
    records: dict[str, dict] = {}
    conflicts: set[str] = set()
    reasons: Counter[str] = Counter()
    types: Counter[str] = Counter()
    rows_seen = 0

    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows_seen += 1
            gene = clean(row.get("gene_symbol"))
            if gene not in TRANSCRIPTS:
                reasons["unsupported_gene"] += 1
                continue
            transcript = clean(row.get("reference_sequence"))
            if transcript != TRANSCRIPTS[gene]:
                reasons["wrong_or_missing_transcript"] += 1
                continue
            c_notation = canonical_c(row.get("cdna", ""), transcript)
            if not c_notation or not re.search(r"delins|del|dup|ins", c_notation, re.IGNORECASE):
                reasons["not_small_indel"] += 1
                continue
            p_notation = canonical_p(row.get("protein", ""))
            if not p_notation:
                reasons["missing_protein_consequence"] += 1
                continue
            grch38 = grch38_from_row(row)
            grch37 = parse_vcf(row.get("genomic_vcf37_source", ""), "GRCh37")
            if not grch37 or not grch38:
                reasons["missing_genomic_mapping"] += 1
                continue
            kind = variant_type(c_notation, p_notation)
            if not kind:
                reasons["unrecognized_indel_type"] += 1
                continue
            key = f"{gene}:{c_notation}"
            record = {
                "gene": gene,
                "input_c_notations": aliases(row, c_notation, transcript, gene, grch38),
                "canonical_c_notation": c_notation,
                "p_notation": p_notation,
                "variant_type": kind,
                "grch37": grch37,
                "grch38": grch38,
                "reference_transcript": transcript,
                "source": {
                    "dataset": "BRCA Exchange",
                    "release": release,
                    "release_date": release_date,
                    "source_names": sorted(filter(None, clean(row.get("source")).split(","))),
                    "ca_id": clean(row.get("ca_id")),
                    "vrs_id": clean(row.get("vr_id")),
                },
            }
            previous = records.get(key)
            if previous:
                identity_fields = ("p_notation", "variant_type", "grch37", "grch38")
                if any(previous[field] != record[field] for field in identity_fields):
                    conflicts.add(key)
                    reasons["conflicting_duplicate"] += 1
                    continue
                previous["input_c_notations"] = sorted(set(previous["input_c_notations"] + record["input_c_notations"]))
            else:
                records[key] = record

    for key in conflicts:
        records.pop(key, None)
    alias_targets: dict[str, set[str]] = {}
    for key, record in records.items():
        gene = record["gene"]
        for notation in record["input_c_notations"]:
            alias_targets.setdefault(f"{gene}:{notation}", set()).add(key)
    ambiguous_aliases = {
        alias: sorted(targets)
        for alias, targets in alias_targets.items()
        if len(targets) > 1
    }
    # An ambiguous alternative spelling is unsafe for automated lookup. The
    # canonical record itself remains available by its exact canonical key.
    for record in records.values():
        gene = record["gene"]
        record["input_c_notations"] = [
            notation
            for notation in record["input_c_notations"]
            if f"{gene}:{notation}" not in ambiguous_aliases
        ]
    reasons["ambiguous_alias_removed"] += len(ambiguous_aliases)
    records = dict(sorted(records.items()))
    types.update(record["variant_type"] for record in records.values())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes((json.dumps(records, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    metadata = {
        "dataset": "Normalized BRCA1/2 small-indel reference snapshot",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "validated_reference_snapshot",
        "source": "BRCA Exchange official release archive",
        "source_release": release,
        "source_release_date": release_date,
        "source_file": source.name,
        "source_sha256": sha256(source),
        "reference_transcripts": TRANSCRIPTS,
        "rows_seen": rows_seen,
        "records": len(records),
        "variant_types": dict(sorted(types.items())),
        "excluded": dict(sorted(reasons.items())),
        "ambiguous_conflicts": len(conflicts),
        "ambiguous_aliases": len(ambiguous_aliases),
        "ambiguous_alias_examples": dict(list(sorted(ambiguous_aliases.items()))[:20]),
        "index_sha256": sha256(output),
    }
    metadata_path.write_bytes((json.dumps(metadata, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--release", default="70")
    parser.add_argument("--release-date", default="2026-03-08")
    args = parser.parse_args()
    print(json.dumps(build(args.source, args.output, args.metadata, args.release, args.release_date), indent=2))


if __name__ == "__main__":
    main()
