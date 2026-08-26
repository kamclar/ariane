# ============================================================
# ClinGen Evidence Repository lookup
# ENIGMA VCEP classifications with evidence codes
# Docs: https://erepo.clinicalgenome.org/evrepo/api
# ============================================================
from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.parse
import urllib.request

from backend.gene_policy import external_evidence_config, reference_transcript
from backend.version import ARIANE_VERSION

EREPO_BASE  = 'https://erepo.clinicalgenome.org/evrepo/api'
EREPO_CACHE: Dict[str, dict] = {}


def clingen_erepo_lookup(gene: str, c_notation: str) -> dict:
    """
    Look up ENIGMA VCEP classification in ClinGen Evidence Repository.
    Returns classification + evidence codes if found, else status=not_found.

    Coverage depends on the VCEP affiliate configured for the active gene.
    """
    key = f'{gene}:{c_notation}'
    if key in EREPO_CACHE:
        return EREPO_CACHE[key]

    tx = reference_transcript(gene)
    hgvs = f'{tx}:{c_notation}'
    affiliate = external_evidence_config(gene)["clingen_erepo_affiliate"]

    url = (
        f"{EREPO_BASE}/classifications"
        f"?hgvs={urllib.parse.quote(hgvs)}"
        f"&affiliate={urllib.parse.quote(affiliate)}"
    )
    try:
        req = urllib.request.Request(
            url,
            headers={'Accept': 'application/json', 'User-Agent': f'ARIANE/{ARIANE_VERSION}'}
        )
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        result = {'status': 'api_error', 'error': str(e)}
        # Network/service failures are transient and must not become sticky.
        return result

    items = data.get('variantInterpretations', [])
    if not items:
        result = {'status': 'not_found', 'hgvs': hgvs}
        EREPO_CACHE[key] = result
        return result

    if len(items) > 1:
        result = {
            'status': 'ambiguous',
            'error': (
                f"ClinGen ERepo returned {len(items)} interpretations for the exact HGVS/affiliate query; "
                "no record was selected"
            ),
            'candidate_caids': [item.get('caid', '') for item in items],
        }
        EREPO_CACHE[key] = result
        return result

    # Exactly one interpretation is safe to use.
    item = items[0]
    guidelines = item.get('guidelines', [])
    classification = ''
    evidence_codes = []
    summary_text = ''

    for g in guidelines:
        classification = g.get('outcome', {}).get('label', '')
        for code in g.get('evidenceCodes', []):
            evidence_codes.append({
                'code':   code.get('label', ''),
                'status': code.get('status', ''),
            })
        summary_text = g.get('description', '')

    result = {
        'status':         'ok',
        'caid':           item.get('caid', ''),
        'classification': classification,
        'evidence_codes': evidence_codes,
        'summary_text':   summary_text,
    }
    EREPO_CACHE[key] = result
    return result
