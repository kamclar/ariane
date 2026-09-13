# ============================================================
# ClinGen Evidence Repository lookup
# ENIGMA VCEP classifications with evidence codes
# Docs: https://erepo.clinicalgenome.org/evrepo/api
# ============================================================
from typing import Dict
import json
import urllib.parse
import urllib.request

from backend.policy.gene import external_evidence_config, reference_transcript
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
    guideline_versions = []
    cspec_ids = []

    for g in guidelines:
        classification = g.get('outcome', {}).get('label', '')
        if g.get('version'):
            guideline_versions.append(str(g['version']))
        if g.get('cspecId'):
            cspec_ids.append(str(g['cspecId']))
        agents = g.get('agents', [])
        if not isinstance(agents, list):
            return {
                'status': 'schema_error',
                'error': 'ClinGen ERepo guideline agents field is not a list',
            }
        for agent in agents:
            codes = agent.get('evidenceCodes', []) if isinstance(agent, dict) else []
            if not isinstance(codes, list):
                return {
                    'status': 'schema_error',
                    'error': 'ClinGen ERepo evidenceCodes field is not a list',
                }
            for code in codes:
                if not isinstance(code, dict):
                    continue
                evidence_codes.append({
                    'code': code.get('label', ''),
                    'status': code.get('status', ''),
                })

    result = {
        'status':         'ok',
        'caid':           item.get('caid', ''),
        'classification': classification,
        'evidence_codes': evidence_codes,
        'guideline_versions': list(dict.fromkeys(guideline_versions)),
        'cspec_ids': list(dict.fromkeys(cspec_ids)),
        'assertion_id': item.get('@id', ''),
    }
    EREPO_CACHE[key] = result
    return result
