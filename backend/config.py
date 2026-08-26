# ============================================================
# ARIANE configuration
# Thresholds, paths, API URLs
# ============================================================
from pathlib import Path
import os

from backend.gene_policy import (
    GENE_POLICY_MANIFEST_PATH,
    GENE_POLICY_METADATA_PATH,
    active_genes,
    domains_by_gene,
    transcripts_by_gene,
)

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
PROJECT_DIR = BASE_DIR.parent
PANEL_REFERENCE_DIR = PROJECT_DIR / "data" / "reference" / "panel"

TABLE4_PATH      = DATA_DIR / "enigma_table4.json"
TABLE9_PATH      = DATA_DIR / "enigma_table9.json"
ENIGMA_RULE_CATALOG_PATH = DATA_DIR / "enigma_rule_catalog.json"
ENIGMA_RULE_DIAGRAMS_PATH = DATA_DIR / "enigma_rule_diagrams.json"
ENIGMA_REFERENCE_TABLES_PATH = DATA_DIR / "enigma_reference_tables.json"
ST7_PATH         = DATA_DIR / "st7_reference_set.json"
PS1_PROTEIN_REGISTRY_PATH = DATA_DIR / "ps1_protein_reference_registry.json"
ST2_SPLICE_EVIDENCE_PATH = DATA_DIR / "enigma_st2_splice_evidence.json"
EXON_CNV_EVIDENCE_PATH = DATA_DIR / "exon_cnv_evidence.json"
EXON_CNV_EVIDENCE_MANIFEST_PATH = (
    PROJECT_DIR / "data" / "sources" / "enigma" / "exon_cnv_evidence_manifest.json"
)
RESIDUES_PATH    = DATA_DIR / "clinically_important_residues.json"

# ── ENIGMA VCEP v1.2 thresholds ───────────────────────────────────────────
_ACTIVE_GENES = active_genes()
# ── Functional domains (ENIGMA Appendix Tables 3/4) ───────────────────────
# RING starts at aa 2 per Appendix Table 3 (AA start=2)
FUNCTIONAL_DOMAINS = domains_by_gene()

# ── External API URLs ──────────────────────────────────────────────────────
SPLICEAI_API_URL     = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/"
SPLICEAI_API_TIMEOUT = 25
SPLICEAI_API_SLEEP   = 1.5

MYVARIANT_URL    = "https://myvariant.info/v1/variant"
CLINVAR_EUTILS   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR_API_SLEEP = 0.4
EREPO_BASE       = "https://erepo.clinicalgenome.org/evrepo/api"

TRANSCRIPTS = transcripts_by_gene()
ALLOWED_GENES = set(_ACTIVE_GENES)
