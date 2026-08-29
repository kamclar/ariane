"""Production reference-transcript c. HGVS to p. HGVS normalization."""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import re

from hgvs.exceptions import (
    HGVSDataNotAvailableError,
    HGVSInvalidVariantError,
    HGVSNormalizationError,
    HGVSParseError,
    HGVSUnsupportedOperationError,
    HGVSVerifyFailedError,
)
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.variantmapper import VariantMapper

from backend.config import TRANSCRIPTS
from backend.gene_policy import active_genes, normalization_validation_variant
from backend.modules.hgvs import normalize_c_notation, normalize_protein_notation
from backend.modules.hgvs_provider import PanelProvider, load_panel_provider


_UNCERTAIN_STRUCTURAL_C = re.compile(
    r"^c\.\(\d+[+-]\d+_\d+[+-]\d+\)_"
    r"\(\d+[+-]\d+_\d+[+-]\d+\)(?:del|dup)$"
)


class VariantNormalizationError(ValueError):
    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class ProteinConsequence:
    gene: str
    transcript: str
    protein_accession: str
    submitted_c_notation: str
    canonical_c_notation: str
    p_notation: str
    consequence_status: str
    provenance: dict[str, str]


@dataclass(frozen=True)
class HGVSEngine:
    panel: PanelProvider
    parser: Parser
    normalizer: Normalizer
    mapper: VariantMapper

    def c_to_p(self, gene: str, c_notation: str) -> ProteinConsequence:
        gene = (gene or "").strip().upper()
        if gene not in self.panel.gene_to_transcript:
            raise VariantNormalizationError(
                "unsupported_gene", f"No approved reference transcript is configured for {gene}"
            )
        transcript = self.panel.gene_to_transcript[gene]
        if transcript != TRANSCRIPTS.get(gene):
            raise VariantNormalizationError(
                "transcript_policy_mismatch", f"Reference transcript policy mismatch for {gene}"
            )
        submitted_c = normalize_c_notation(c_notation)
        # Biocommons 1.5.7 does not parse the HGVS uncertain-breakpoint form
        # used for exon-level CNVs. Validate that form explicitly and preserve
        # it losslessly; no deterministic protein sequence can be derived.
        if _UNCERTAIN_STRUCTURAL_C.fullmatch(submitted_c):
            return ProteinConsequence(
                gene=gene,
                transcript=transcript,
                protein_accession=self.panel.transcript_to_protein[transcript],
                submitted_c_notation=submitted_c,
                canonical_c_notation=submitted_c,
                p_notation="p.?",
                consequence_status="protein_consequence_unknown",
                provenance=dict(self.panel.provenance),
            )
        try:
            variant_c = self.parser.parse_hgvs_variant(f"{transcript}:{submitted_c}")
        except HGVSParseError as exc:
            raise VariantNormalizationError(
                "invalid_hgvs_syntax", f"Invalid c. HGVS for {gene}: {submitted_c}: {exc}"
            ) from exc
        if variant_c.type != "c" or variant_c.ac != transcript:
            raise VariantNormalizationError(
                "invalid_hgvs_coordinate_system",
                f"Expected {transcript}:c. HGVS, received {variant_c}",
            )
        # Intronic offsets and uncertain structural breakpoints do not define a
        # unique altered coding sequence. They are valid c. descriptions, but
        # their protein consequence must remain unknown until RNA/breakpoint
        # evidence resolves it. Biocommons intentionally does not normalize
        # intronic c. variants without genomic sequence/alignment data.
        if re.search(r"\d[+-]\d|[()]", submitted_c):
            return ProteinConsequence(
                gene=gene,
                transcript=transcript,
                protein_accession=self.panel.transcript_to_protein[transcript],
                submitted_c_notation=submitted_c,
                canonical_c_notation=str(variant_c).split(":", 1)[-1],
                p_notation="p.?",
                consequence_status="protein_consequence_unknown",
                provenance=dict(self.panel.provenance),
            )
        try:
            canonical_variant = self.normalizer.normalize(variant_c)
            protein_variant = self.mapper.c_to_p(
                canonical_variant,
                pro_ac=self.panel.transcript_to_protein[transcript],
            )
        except HGVSVerifyFailedError as exc:
            raise VariantNormalizationError(
                "reference_allele_mismatch",
                f"Reference allele does not match {transcript} for {submitted_c}: {exc}",
            ) from exc
        except HGVSDataNotAvailableError as exc:
            raise VariantNormalizationError(
                "reference_data_unavailable",
                f"Required local reference data are unavailable for {transcript}:{submitted_c}: {exc}",
            ) from exc
        except HGVSUnsupportedOperationError as exc:
            raise VariantNormalizationError(
                "unsupported_hgvs_operation",
                f"The variant cannot be mapped safely to protein HGVS: {transcript}:{submitted_c}: {exc}",
            ) from exc
        except HGVSInvalidVariantError as exc:
            message = str(exc)
            if "reference" in message.lower() and (
                "does not agree" in message.lower() or "does not match" in message.lower()
            ):
                raise VariantNormalizationError(
                    "reference_allele_mismatch",
                    f"Reference allele does not match {transcript} for {submitted_c}: {exc}",
                ) from exc
            raise VariantNormalizationError(
                "invalid_hgvs_variant",
                f"The variant is invalid on {transcript}: {submitted_c}: {exc}",
            ) from exc
        except HGVSNormalizationError as exc:
            raise VariantNormalizationError(
                "invalid_hgvs_variant",
                f"The variant is invalid on {transcript}: {submitted_c}: {exc}",
            ) from exc
        except Exception as exc:
            raise VariantNormalizationError(
                "normalization_engine_failure",
                f"Local HGVS normalization failed for {transcript}:{submitted_c}: "
                f"{type(exc).__name__}: {exc}",
            ) from exc

        canonical_c = str(canonical_variant).split(":", 1)[-1]
        raw_p = str(protein_variant).split(":", 1)[-1]
        p_notation = normalize_protein_notation(raw_p)
        if not p_notation:
            raise VariantNormalizationError(
                "protein_consequence_missing",
                f"HGVS engine returned no protein consequence for {transcript}:{canonical_c}",
            )
        if p_notation == "p.?":
            status = "protein_consequence_unknown"
        elif p_notation.endswith("=)"):
            status = "sequence_derived_synonymous"
        else:
            status = "sequence_derived"
        return ProteinConsequence(
            gene=gene,
            transcript=transcript,
            protein_accession=self.panel.transcript_to_protein[transcript],
            submitted_c_notation=submitted_c,
            canonical_c_notation=canonical_c,
            p_notation=p_notation,
            consequence_status=status,
            provenance=dict(self.panel.provenance),
        )


@lru_cache(maxsize=1)
def load_hgvs_engine() -> HGVSEngine:
    panel = load_panel_provider()
    mapper = VariantMapper(
        panel.data_provider,
        replace_reference=True,
        prevalidation_level="EXTRINSIC",
    )
    normalizer = Normalizer(
        panel.data_provider,
        # Transcript c. variants may legitimately span exon/exon or UTR/CDS
        # boundaries. The local cdot alignment data are checksum-verified, so
        # the normalizer is allowed to traverse those annotated boundaries.
        cross_boundaries=True,
        shuffle_direction=3,
        validate=True,
        variantmapper=mapper,
    )
    return HGVSEngine(panel=panel, parser=Parser(), normalizer=normalizer, mapper=mapper)


def derive_protein_consequence(gene: str, c_notation: str) -> ProteinConsequence:
    return load_hgvs_engine().c_to_p(gene, c_notation)


def validate_hgvs_engine() -> None:
    engine = load_hgvs_engine()
    for gene in active_genes():
        validation = normalization_validation_variant(gene)
        c_notation = validation["c_notation"]
        expected_p = validation["p_notation"]
        result = engine.c_to_p(gene, c_notation)
        if result.p_notation != expected_p:
            raise RuntimeError(
                f"HGVS engine startup validation failed for {gene} {c_notation}: "
                f"expected {expected_p}, found {result.p_notation}"
            )
