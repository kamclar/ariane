"""Strict in-memory sequence fetcher for the versioned ARIANE panel bundle."""
from __future__ import annotations

import re
from pathlib import Path

from hgvs.exceptions import HGVSDataNotAvailableError


class PanelSequenceFetcher:
    """Serve exact accession-version sequences to hgvs/cdot.

    The panel FASTA files are deliberately small. Loading them once per worker
    avoids platform-specific FASTA index dependencies while retaining strict
    bounds and accession checks.
    """

    def __init__(self, transcript_fasta: Path, protein_fasta: Path):
        self.source = f"{transcript_fasta};{protein_fasta}"
        self._sequences = {}
        self._sequences.update(self._read_fasta(transcript_fasta, r"[ACGT]+"))
        proteins = self._read_fasta(protein_fasta, r"[ABCDEFGHIKLMNPQRSTVWXYZ*]+")
        overlap = set(self._sequences) & set(proteins)
        if overlap:
            raise RuntimeError(f"Duplicate reference accessions: {sorted(overlap)}")
        self._sequences.update(proteins)
        self._data_provider = None

    @staticmethod
    def _read_fasta(path: Path, sequence_pattern: str) -> dict[str, str]:
        if not path.is_file():
            raise RuntimeError(f"Required reference FASTA is missing: {path}")
        records: dict[str, list[str]] = {}
        accession = None
        try:
            with path.open(encoding="ascii") as handle:
                for line_number, raw_line in enumerate(handle, 1):
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        accession = line[1:].split()[0]
                        if not accession or accession in records:
                            raise RuntimeError(
                                f"Invalid or duplicate FASTA accession at {path}:{line_number}"
                            )
                        records[accession] = []
                    elif accession is None:
                        raise RuntimeError(f"FASTA sequence precedes header at {path}:{line_number}")
                    else:
                        records[accession].append(line.upper())
        except OSError as exc:
            raise RuntimeError(f"Required reference FASTA cannot be read: {path}: {exc}") from exc
        sequences = {key: "".join(value) for key, value in records.items()}
        if not sequences:
            raise RuntimeError(f"Required reference FASTA has no records: {path}")
        for key, sequence in sequences.items():
            if not sequence or not re.fullmatch(sequence_pattern, sequence):
                raise RuntimeError(f"Invalid reference sequence for {key} in {path}")
        return sequences

    def set_data_provider(self, data_provider) -> None:
        self._data_provider = data_provider

    def fetch_seq(self, accession: str, start_i: int | None = None, end_i: int | None = None) -> str:
        try:
            sequence = self._sequences[accession]
        except KeyError as exc:
            raise HGVSDataNotAvailableError(
                f"Reference sequence is unavailable for exact accession {accession}"
            ) from exc
        start = 0 if start_i is None else start_i
        end = len(sequence) if end_i is None else end_i
        if not isinstance(start, int) or not isinstance(end, int):
            raise HGVSDataNotAvailableError(f"Non-integer sequence interval for {accession}")
        if start < 0 or end < start or end > len(sequence):
            raise HGVSDataNotAvailableError(
                f"Sequence interval {start}:{end} is outside {accession} length {len(sequence)}"
            )
        return sequence[start:end]

    def accessions(self) -> frozenset[str]:
        return frozenset(self._sequences)

    def sequence_length(self, accession: str) -> int:
        try:
            return len(self._sequences[accession])
        except KeyError as exc:
            raise HGVSDataNotAvailableError(f"Unknown reference accession {accession}") from exc
