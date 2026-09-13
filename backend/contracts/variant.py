"""Variant input and normalization HTTP schemas."""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from backend.policy.gene import active_genes
from backend.variant_processing.variant_input import normalize_variant_input


class VariantRequest(BaseModel):
    gene: str
    c_notation: str
    p_notation: str = ""
    dup_type: str = "Unknown"
    assembly: Literal["GRCh37", "GRCh38"] | None = None
    submitted_notation: str = ""
    reference_transcript: str = ""
    normalization_source: str = ""
    consequence_status: str = ""
    normalization_provenance: dict[str, str] = Field(default_factory=dict)
    protein_consequence_explanation: str = ""

    @model_validator(mode="before")
    @classmethod
    def normalize_hgvs_fields(cls, data):
        if not isinstance(data, dict):
            return data
        gene = str(data.get("gene", "")).strip().upper()
        raw_notation = str(data.get("c_notation", "")).strip()
        normalized = normalize_variant_input(
            gene,
            raw_notation,
            assembly=data.get("assembly"),
            p_notation=data.get("p_notation"),
        )
        return {
            **data,
            "gene": normalized.gene,
            "submitted_notation": raw_notation,
            "c_notation": normalized.c_notation,
            "p_notation": normalized.p_notation,
            "assembly": normalized.assembly or None,
            "reference_transcript": normalized.reference_transcript,
            "normalization_source": normalized.normalization_source,
            "consequence_status": normalized.consequence_status,
            "normalization_provenance": normalized.normalization_provenance or {},
            "protein_consequence_explanation": normalized.protein_consequence_explanation,
        }

    @field_validator("gene")
    @classmethod
    def validate_gene(cls, value):
        gene = value.strip().upper()
        supported = set(active_genes())
        if gene not in supported:
            raise ValueError(f"Gene must be one of: {', '.join(sorted(supported))}")
        return gene

    @field_validator("c_notation")
    @classmethod
    def validate_c_notation(cls, value):
        notation = value.strip()
        if not notation.startswith("c."):
            raise ValueError("c. notation must start with 'c.', for example c.4185G>A")
        substitution = r"^c\.[-*]?\d+(?:[+-]\d+)?[ACGT]>[ACGT]$"
        equality = r"^c\.[-*]?\d+(?:[+-]\d+)?[ACGT]?=$"
        sequence_change = r"^c\.[0-9*+_?()\-]+(?:delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+)$"
        if not any(
            re.fullmatch(pattern, notation, re.IGNORECASE)
            for pattern in (substitution, equality, sequence_change)
        ):
            raise ValueError(
                "Unrecognised c. notation. Use HGVS format, for example "
                "c.4185G>A, c.68_69delAG, or c.212+1G>T"
            )
        return notation

    @field_validator("p_notation")
    @classmethod
    def validate_p_notation(cls, value):
        notation = value.strip()
        if not notation:
            raise ValueError(
                "p. notation is required, for example p.(Gln1395=); "
                "use p.? when the protein consequence is unknown"
            )
        if not notation.startswith("p."):
            raise ValueError("p. notation must start with 'p.', for example p.(Gln1395=)")
        protein = r"^(?:p\.\?|p\.\([A-Za-z0-9_?=*]+\))$"
        if not re.fullmatch(protein, notation):
            raise ValueError(
                "Unrecognised p. notation. Use HGVS format, for example "
                "p.(Cys61Gly), p.(Gln1395=), or p.(Glu23ValfsTer17)"
            )
        return notation

    @field_validator("dup_type")
    @classmethod
    def validate_dup_type(cls, value):
        dup_type = value.strip().title()
        if dup_type not in ("Unknown", "Tandem"):
            raise ValueError("dup_type must be Unknown or Tandem")
        return dup_type


class VariantNormalizationResponse(BaseModel):
    gene: str
    submitted_notation: str
    c_notation: str
    p_notation: str
    reference_transcript: str
    normalization_source: str
    consequence_status: str = ""
    normalization_provenance: dict[str, str] = Field(default_factory=dict)
    protein_consequence_explanation: str = ""
    assembly: Literal["GRCh37", "GRCh38"] | None = None
