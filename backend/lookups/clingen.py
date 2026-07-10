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
import urllib.request
import urllib.parse

EREPO_BASE  = 'https://erepo.clinicalgenome.org/evrepo/api'
EREPO_CACHE: Dict[str, dict] = {}


def clingen_erepo_lookup(gene: str, c_notation: str) -> dict:
    """
    Look up ENIGMA VCEP classification in ClinGen Evidence Repository.
    Returns classification + evidence codes if found, else status=not_found.

    Coverage is limited - ENIGMA VCEP has ~140 BRCA1/2 submissions in ERepo.
    Most variants will return not_found; use ClinVar for broader coverage.
    """
    key = f'{gene}:{c_notation}'
    if key in EREPO_CACHE:
        return EREPO_CACHE[key]

    tx = 'NM_007294.4' if gene == 'BRCA1' else 'NM_000059.4'
    hgvs = f'{tx}:{c_notation}'

    url = (
        f"{EREPO_BASE}/classifications"
        f"?hgvs={urllib.parse.quote(hgvs)}"
        f"&affiliate={urllib.parse.quote('ENIGMA BRCA1 and BRCA2 VCEP')}"
    )
    try:
        req = urllib.request.Request(
            url,
            headers={'Accept': 'application/json', 'User-Agent': 'BRCA-ACMG/1.7.0'}
        )
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        result = {'status': 'api_error', 'error': str(e)}
        EREPO_CACHE[key] = result
        return result

    items = data.get('variantInterpretations', [])
    if not items:
        result = {'status': 'not_found', 'hgvs': hgvs}
        EREPO_CACHE[key] = result
        return result

    # parse first result
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


if __name__ == "__main__":
    print('ClinGen ERepo lookup functions loaded.')
