# ============================================================
# ARIANE configuration
# Thresholds, paths, API URLs
# ============================================================
from pathlib import Path
import os

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

TABLE4_PATH      = DATA_DIR / "enigma_table4.json"
TABLE9_PATH      = DATA_DIR / "enigma_table9.json"
ST7_PATH         = DATA_DIR / "st7_reference_set.json"
RESIDUES_PATH    = DATA_DIR / "clinically_important_residues.json"
GNOMAD_CACHE_DIR = DATA_DIR / "gnomad_brca_cache"
SPLICEAI_CACHE       = DATA_DIR / "spliceai_api_cache.json"
BAYESDEL_CACHE_FILE  = DATA_DIR / "bayesdel_cache.json"
COORDINATES_CACHE_FILE = DATA_DIR / "coordinates_cache.json"

# ── ENIGMA VCEP v1.2 thresholds ───────────────────────────────────────────
BA1_FAF_THRESHOLD        = 0.001
BS1_STRONG_THRESHOLD     = 0.0001
BS1_SUPPORTING_THRESHOLD = 0.00002
PM2_MIN_COVERAGE         = 25

BAYESDEL_THRESHOLDS = {
    "BRCA1": {"pp3": 0.28, "bp4": 0.15},
    "BRCA2": {"pp3": 0.30, "bp4": 0.18},
}

SPLICEAI_PP3_THRESHOLD = 0.2
SPLICEAI_BP4_THRESHOLD = 0.1

PP4_LR_THRESHOLDS = {
    "Very Strong": 350,
    "Strong":      18.7,
    "Moderate":    4.3,
    "Supporting":  2.08,
}
BP5_LR_THRESHOLDS = {
    "Very Strong": 0.00285,
    "Strong":      0.05,
    "Moderate":    0.23,
    "Supporting":  0.48,
}
MULTIFACTORIAL_PRIOR = 0.10

CLASS_THRESHOLDS = {5: 10, 4: 6, 3: -1, 2: -6}

# ── Functional domains (ENIGMA Appendix Tables 3/4) ───────────────────────
# RING starts at aa 2 per Appendix Table 3 (AA start=2)
FUNCTIONAL_DOMAINS = {
    "BRCA1": {
        "RING":        (2, 101),
        "coiled_coil": (1391, 1424),
        "BRCT":        (1650, 1857),
    },
    "BRCA2": {
        "PALB2_binding": (10, 40),
        "DBD":           (2481, 3186),
    },
}

# ── External API URLs ──────────────────────────────────────────────────────
SPLICEAI_API_URL     = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/"
SPLICEAI_API_TIMEOUT = 60
SPLICEAI_API_SLEEP   = 1.5

MYVARIANT_URL    = "https://myvariant.info/v1/variant"
CLINVAR_EUTILS   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR_API_SLEEP = 0.4
EREPO_BASE       = "https://erepo.clinicalgenome.org/evrepo/api"

TRANSCRIPTS = {
    "BRCA1": "NM_007294.4",
    "BRCA2": "NM_000059.4",
}

ALLOWED_GENES = {"BRCA1", "BRCA2"}
