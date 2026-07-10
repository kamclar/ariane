# ============================================================
# ClinVar lookup - aggregate classification + all submitters
# Source: NCBI eutils efetch VCV XML
# Docs: https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/
# ============================================================
from typing import Optional, Dict, List, Tuple
from pathlib import Path
import json
import re
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from typing import Optional, Dict, List

from backend.lookups.coordinates import resolve_variant, get_grch38

CLINVAR_EUTILS   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR_CACHE: Dict[str, dict] = {}  # gene:c_notation -> parsed result
CLINVAR_API_SLEEP = 0.4  # NCBI rate limit: max 3 requests/sec without API key


def clinvar_review_stars(review_status: str) -> int:
    """Map the official ClinVar aggregate review status to its star level."""
    status = (review_status or "").strip().lower()
    if "no assertion criteria" in status or "no classification" in status:
        return 0
    if "practice guideline" in status:
        return 4
    if "expert panel" in status:
        return 3
    if "multiple submitters" in status and "no conflicts" in status:
        return 2
    if "criteria provided" in status:
        return 1
    return 0


def clinvar_search_variation_id(gene: str, c_notation: str) -> Optional[str]:
    """
    Search ClinVar by SPDI (preferred) or HGVS fallback.
    SPDI is unambiguous; HGVS esearch can return multiple hits for the same c. notation.

    NOTE: n_submitters is taken from the VCV aggregate (NumberOfSubmitters),
    which counts unique submitters across ALL disease associations (RCV records).
    This may differ from submitter counts shown in Excel/tutorial spreadsheets,
    which typically reflect a single RCV (one disease condition).

    To verify manually:
      1. Go to https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/
      2. In "Submitted interpretations and evidence", check count per condition (RCV)
         vs. total across all conditions (VCV).
      3. The VCV page header shows the total; individual RCV tabs show per-disease counts.
    """
    # ── SPDI lookup (preferred) ────────────────────────────────────────────
    resolved = {}
    resolved_variant = resolve_variant(gene, c_notation)
    if resolved_variant and resolved_variant.status != "failed":
        resolved[f"{gene}:{c_notation}"] = resolved_variant
    coords = get_grch38(resolved, gene, c_notation)
    if coords:
        chrom     = str(coords["chrom"])
        pos       = int(coords["pos"]) - 1  # SPDI is 0-based
        ref       = coords["ref"]
        alt       = coords["alt"]
        chrom_acc = {"17": "NC_000017.11", "13": "NC_000013.11"}.get(chrom)
        if chrom_acc:
            spdi = f"{chrom_acc}:{pos}:{ref}:{alt}"
            url  = (
                f"{CLINVAR_EUTILS}/esearch.fcgi"
                f"?db=clinvar&retmode=json"
                f"&term={urllib.parse.quote(spdi)}"
            )
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "BRCA-ACMG/1.7.0"})
                with urllib.request.urlopen(req, timeout=20) as resp:
                    data = json.loads(resp.read())
                ids = data.get("esearchresult", {}).get("idlist", [])
                if len(ids) == 1:
                    return ids[0]
            except Exception:
                pass

    # ── HGVS fallback - filter by exact c. notation in title ──────────────
    tx   = "NM_007294.4" if gene == "BRCA1" else "NM_000059.4"
    hgvs = f"{tx}:{c_notation}"
    url  = (
        f"{CLINVAR_EUTILS}/esearch.fcgi"
        f"?db=clinvar&retmode=json&retmax=10"
        f"&term={urllib.parse.quote(hgvs)}"
    )
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "BRCA-ACMG/1.7.0"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read())
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None
        if len(ids) == 1:
            return ids[0]

        # multiple hits - fetch summaries and match exact c. notation in title
        time.sleep(CLINVAR_API_SLEEP)
        sum_url = (
            f"{CLINVAR_EUTILS}/esummary.fcgi"
            f"?db=clinvar&id={','.join(ids)}&retmode=json"
        )
        req2 = urllib.request.Request(sum_url, headers={"User-Agent": "BRCA-ACMG/1.7.0"})
        with urllib.request.urlopen(req2, timeout=20) as resp2:
            sum_data = json.loads(resp2.read())
        result = sum_data.get("result", {})
        for vid in ids:
            title = result.get(vid, {}).get("title", "")
            if c_notation in title:
                return vid
        return ids[0]

    except Exception:
        return None


def clinvar_fetch_vcv(variation_id: str) -> Optional[ET.Element]:
    """Fetch VCV XML for a variation ID and return the VariationArchive element."""
    url = (
        f"{CLINVAR_EUTILS}/efetch.fcgi"
        f"?db=clinvar&rettype=vcv&is_variationid"
        f"&id={variation_id}&release_date=current"
    )
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "BRCA-ACMG/1.7.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            root = ET.fromstring(resp.read())
        return root.find("VariationArchive")
    except Exception:
        return None


def clinvar_parse_vcv(va: ET.Element) -> dict:
    """Parse a VariationArchive element into a structured dict."""
    cr = va.find("ClassifiedRecord")
    if cr is None:
        return {"status": "no_classified_record"}

    agg      = cr.find("Classifications/GermlineClassification")
    no_class = cr.find("Classifications/NoClassification")

    if agg is not None:
        agg_desc  = agg.find("Description")
        aggregate = {
            "classification": agg.findtext("Description")           if agg is not None else "",
            "review_status":  agg.findtext("ReviewStatus")          if agg is not None else "",
            "date_evaluated": agg_desc.get("DateLastEvaluated", "") if agg_desc is not None else "",
            "n_submitters":   int(agg.get("NumberOfSubmitters", 0)) if agg is not None else 0,
        }
    elif no_class is not None:
        aggregate = {
            "classification": "no classification provided",
            "review_status":  no_class.findtext("ReviewStatus") or "",
            "date_evaluated": "",
            "n_submitters":   int(no_class.get("NumberOfSubmitters", 0)),
        }
    else:
        aggregate = {
            "classification": "",
            "review_status":  "",
            "date_evaluated": "",
            "n_submitters":   0,
        }

    submissions       = []
    enigma_submission = None
    for ca in cr.findall("ClinicalAssertionList/ClinicalAssertion"):
        scv_elem    = ca.find("ClinVarAccession")
        scv         = scv_elem.get("Accession", "")            if scv_elem is not None else ""
        org         = scv_elem.get("SubmitterName", "")        if scv_elem is not None else ""
        org_cat     = scv_elem.get("OrganizationCategory", "") if scv_elem is not None else ""
        cl_elem     = ca.find("Classification")
        classif     = cl_elem.findtext("GermlineClassification") if cl_elem is not None else ""
        date_eval   = cl_elem.get("DateLastEvaluated", "")      if cl_elem is not None else ""
        review      = cl_elem.findtext("ReviewStatus")          if cl_elem is not None else ""
        comment     = cl_elem.findtext("Comment")               if cl_elem is not None else ""
        contributes = ca.get("ContributesToAggregateClassification", "")

        entry = {
            "scv":          scv,
            "org":          org,
            "org_cat":      org_cat,
            "class":        classif,
            "date_eval":    date_eval,
            "review":       review,
            "comment":      comment or "",
            "contributes":  contributes == "true",
            "is_enigma_ep": "ENIGMA" in org or "Evidence-based Network" in org,
        }
        submissions.append(entry)
        if entry["is_enigma_ep"]:
            enigma_submission = entry

    classes      = set(s["class"] for s in submissions if s["class"])
    has_conflict = len(classes) > 1

    return {
        "status":            "ok",
        "variation_id":      va.get("VariationID", ""),
        "accession":         va.get("Accession", ""),
        "aggregate":         aggregate,
        "submissions":       submissions,
        "enigma_submission": enigma_submission,
        "has_conflict":      has_conflict,
        "unique_classes":    sorted(classes),
    }


def clinvar_lookup(gene: str, c_notation: str) -> dict:
    """
    Main entry point. Searches ClinVar by SPDI or HGVS and returns parsed VCV.
    Results are cached in-memory per session.
    """
    key = f"{gene}:{c_notation}"
    if key in CLINVAR_CACHE:
        return CLINVAR_CACHE[key]

    tx   = "NM_007294.4" if gene == "BRCA1" else "NM_000059.4"
    hgvs = f"{tx}:{c_notation}"

    time.sleep(CLINVAR_API_SLEEP)
    variation_id = clinvar_search_variation_id(gene, c_notation)
    if not variation_id:
        result = {"status": "not_found", "hgvs": hgvs}
        CLINVAR_CACHE[key] = result
        return result

    time.sleep(CLINVAR_API_SLEEP)
    va = clinvar_fetch_vcv(variation_id)
    if va is None:
        result = {"status": "fetch_failed", "variation_id": variation_id}
        CLINVAR_CACHE[key] = result
        return result

    result = clinvar_parse_vcv(va)
    CLINVAR_CACHE[key] = result
    return result


if __name__ == "__main__":
    print("ClinVar lookup functions loaded.")
