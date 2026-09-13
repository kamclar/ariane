"""Immutable paths for versioned reference datasets."""

from pathlib import Path


BACKEND_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = BACKEND_DIR.parent
DATA_DIR = BACKEND_DIR / "data"
PANEL_REFERENCE_DIR = PROJECT_DIR / "data" / "reference" / "panel"

TABLE4_PATH = DATA_DIR / "enigma_table4.json"
TABLE9_PATH = DATA_DIR / "enigma_table9.json"
ENIGMA_RULE_CATALOG_PATH = DATA_DIR / "enigma_rule_catalog.json"
ENIGMA_RULE_DIAGRAMS_PATH = DATA_DIR / "enigma_rule_diagrams.json"
ENIGMA_REFERENCE_TABLES_PATH = DATA_DIR / "enigma_reference_tables.json"
ST7_PATH = DATA_DIR / "st7_reference_set.json"
PS1_PROTEIN_REGISTRY_PATH = DATA_DIR / "ps1_protein_reference_registry.json"
ENIGMA_EREPO_VCEP_REGISTRY_PATH = DATA_DIR / "enigma_erepo_vcep_registry.json"
ENIGMA_EREPO_VCEP_METADATA_PATH = DATA_DIR / "enigma_erepo_vcep_registry.metadata.json"
ST2_SPLICE_EVIDENCE_PATH = DATA_DIR / "enigma_st2_splice_evidence.json"
EXON_CNV_EVIDENCE_PATH = DATA_DIR / "exon_cnv_evidence.json"
EXON_CNV_EVIDENCE_MANIFEST_PATH = (
    PROJECT_DIR / "data" / "sources" / "enigma" / "exon_cnv_evidence_manifest.json"
)
RESIDUES_PATH = DATA_DIR / "clinically_important_residues.json"
