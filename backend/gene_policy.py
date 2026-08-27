"""Strict, versioned gene and VCEP policy registry.

The registry is the single runtime source for active genes, reference
transcripts, rule applicability and numeric decision thresholds. Unknown genes,
missing fields and checksum drift are fatal. There are no inherited policies or
generic threshold fallbacks.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping


DATA_DIR = Path(__file__).resolve().parent / "data"
GENE_POLICY_MANIFEST_PATH = DATA_DIR / "gene_policy_manifest.json"
GENE_POLICY_METADATA_PATH = DATA_DIR / "gene_policy_manifest.metadata.json"

_TRANSCRIPT_RE = re.compile(r"^N[MR]_\d+\.\d+$")
_PROTEIN_RE = re.compile(r"^NP_\d+\.\d+$")
_STANDARD_CRITERION_RE = re.compile(
    r"^(?:PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])$"
)


class GenePolicyError(RuntimeError):
    """Raised when the configured gene policy cannot be used safely."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise GenePolicyError(f"Required {label} is missing: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GenePolicyError(f"Required {label} cannot be loaded: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise GenePolicyError(f"Required {label} must contain a JSON object: {path}")
    return value


def _require_number(value: Any, location: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise GenePolicyError(f"Gene policy {location} must be numeric")
    return float(value)


def _require_keys(value: Mapping[str, Any], keys: set[str], location: str) -> None:
    missing = sorted(keys - set(value))
    if missing:
        raise GenePolicyError(f"Gene policy {location} is incomplete: missing {missing}")


def validate_gene_policy_payload(
    manifest: Mapping[str, Any], metadata: Mapping[str, Any] | None = None,
    *, manifest_bytes: bytes | None = None,
) -> None:
    """Validate the complete registry and, when supplied, its file checksum."""
    _require_keys(
        manifest,
        {"schema_version", "manifest_id", "manifest_version", "status", "policies", "genes"},
        "manifest",
    )
    if manifest["schema_version"] != 1 or manifest["status"] != "active":
        raise GenePolicyError("Gene policy manifest must be active schema version 1")
    policies = manifest["policies"]
    genes = manifest["genes"]
    if not isinstance(policies, dict) or not policies:
        raise GenePolicyError("Gene policy manifest has no policies")
    if not isinstance(genes, dict) or not genes:
        raise GenePolicyError("Gene policy manifest has no genes")

    active_policies = 0
    for policy_id, policy in policies.items():
        if not isinstance(policy, dict):
            raise GenePolicyError(f"Gene policy policies/{policy_id} must be an object")
        _require_keys(
            policy,
            {
                "name", "version", "runtime_policy_id", "implementation_profile",
                "external_evidence", "thresholds", "supported_rule_codes",
                "criterion_applicability",
            },
            f"policies/{policy_id}",
        )
        if not re.fullmatch(r"[a-z0-9_]+", str(policy["implementation_profile"])):
            raise GenePolicyError(
                f"Policy {policy_id} has an invalid implementation profile"
            )
        external_evidence = policy["external_evidence"]
        if not isinstance(external_evidence, dict):
            raise GenePolicyError(
                f"Policy {policy_id} external_evidence must be an object"
            )
        _require_keys(
            external_evidence,
            {"clingen_erepo_affiliate"},
            f"policies/{policy_id}/external_evidence",
        )
        if not str(external_evidence["clingen_erepo_affiliate"]).strip():
            raise GenePolicyError(
                f"Policy {policy_id} has no ClinGen ERepo affiliate"
            )
        active_policies += 1
        thresholds = policy["thresholds"]
        _require_keys(
            thresholds,
            {"spliceai", "clinical_likelihood_ratio", "population_frequency", "mixed_evidence_points"},
            f"policies/{policy_id}/thresholds",
        )
        splice = thresholds["spliceai"]
        _require_keys(
            splice,
            {"no_impact_max_inclusive", "impact_min_inclusive"},
            f"policies/{policy_id}/thresholds/spliceai",
        )
        no_impact = _require_number(splice["no_impact_max_inclusive"], "spliceai/no_impact")
        impact = _require_number(splice["impact_min_inclusive"], "spliceai/impact")
        if not 0 <= no_impact < impact <= 1:
            raise GenePolicyError("SpliceAI policy must satisfy 0 <= no impact < impact <= 1")

        clinical = thresholds["clinical_likelihood_ratio"]
        _require_keys(clinical, {"pp4", "bp5"}, f"policies/{policy_id}/clinical_likelihood_ratio")
        pp4 = clinical["pp4"]
        bp5 = clinical["bp5"]
        pp4_values = [
            _require_number(pp4[key], f"clinical_likelihood_ratio/pp4/{key}")
            for key in (
                "supporting_min_inclusive", "moderate_min_inclusive",
                "strong_min_inclusive", "very_strong_min_inclusive",
            )
        ]
        bp5_values = [
            _require_number(bp5[key], f"clinical_likelihood_ratio/bp5/{key}")
            for key in (
                "supporting_max_inclusive", "moderate_max_inclusive",
                "strong_max_inclusive", "very_strong_max_inclusive",
            )
        ]
        if pp4_values != sorted(pp4_values):
            raise GenePolicyError("PP4 likelihood-ratio thresholds must increase with strength")
        if bp5_values != sorted(bp5_values, reverse=True):
            raise GenePolicyError("BP5 likelihood-ratio thresholds must decrease with strength")

        population = thresholds["population_frequency"]
        _require_keys(
            population,
            {
                "gnomad_policy_id", "ba1_faf95_min_exclusive",
                "bs1_strong_faf95_min_exclusive", "bs1_supporting_faf95_min_exclusive",
                "ba1_bs1_minimum_mean_depth", "pm2_minimum_mean_depth",
            },
            f"policies/{policy_id}/thresholds/population_frequency",
        )
        ba1 = _require_number(population["ba1_faf95_min_exclusive"], "population/ba1")
        bs1s = _require_number(population["bs1_strong_faf95_min_exclusive"], "population/bs1_strong")
        bs1p = _require_number(population["bs1_supporting_faf95_min_exclusive"], "population/bs1_supporting")
        if not 0 < bs1p < bs1s < ba1:
            raise GenePolicyError("Population thresholds must satisfy BS1 Supporting < BS1 Strong < BA1")
        point_thresholds = thresholds["mixed_evidence_points"]
        _require_keys(
            point_thresholds,
            {
                "pathogenic_min_inclusive", "likely_pathogenic_min_inclusive",
                "vus_min_inclusive", "likely_benign_min_inclusive",
            },
            f"policies/{policy_id}/thresholds/mixed_evidence_points",
        )
        point_values = [
            _require_number(point_thresholds[key], f"mixed_evidence_points/{key}")
            for key in (
                "pathogenic_min_inclusive", "likely_pathogenic_min_inclusive",
                "vus_min_inclusive", "likely_benign_min_inclusive",
            )
        ]
        if point_values != sorted(point_values, reverse=True) or len(set(point_values)) != 4:
            raise GenePolicyError("Mixed-evidence class thresholds must strictly decrease")
        rules = policy["supported_rule_codes"]
        if not isinstance(rules, list) or not rules or len(rules) != len(set(rules)):
            raise GenePolicyError(f"Policy {policy_id} supported_rule_codes must be a unique non-empty list")
        applicability = policy["criterion_applicability"]
        if not isinstance(applicability, dict):
            raise GenePolicyError(
                f"Policy {policy_id} criterion_applicability must be an object"
            )
        _require_keys(
            applicability,
            {"not_used"},
            f"policies/{policy_id}/criterion_applicability",
        )
        not_used = applicability["not_used"]
        if not isinstance(not_used, list):
            raise GenePolicyError(
                f"Policy {policy_id} criterion_applicability/not_used must be a list"
            )
        identities: set[tuple[str, str]] = set()
        for index, item in enumerate(not_used):
            location = f"policies/{policy_id}/criterion_applicability/not_used/{index}"
            if not isinstance(item, dict):
                raise GenePolicyError(f"Gene policy {location} must be an object")
            _require_keys(item, {"code", "reason", "source_section"}, location)
            code = str(item["code"]).strip().upper()
            scope = str(item.get("scope") or "").strip()
            if not _STANDARD_CRITERION_RE.fullmatch(code):
                raise GenePolicyError(f"Gene policy {location} has invalid criterion code")
            identity = (code, scope)
            if identity in identities:
                raise GenePolicyError(
                    f"Policy {policy_id} has a duplicate not-used criterion {identity}"
                )
            identities.add(identity)
            if code in set(rules):
                raise GenePolicyError(
                    f"Policy {policy_id} marks supported rule {code} as not used"
                )
            if not str(item["reason"]).strip() or not str(item["source_section"]).strip():
                raise GenePolicyError(
                    f"Gene policy {location} requires a reason and source section"
                )

    active_genes = 0
    for raw_gene, gene in genes.items():
        if raw_gene != raw_gene.upper() or not re.fullmatch(r"[A-Z0-9-]+", raw_gene):
            raise GenePolicyError(f"Invalid gene symbol in policy manifest: {raw_gene!r}")
        if not isinstance(gene, dict):
            raise GenePolicyError(f"Gene policy genes/{raw_gene} must be an object")
        _require_keys(
            gene,
            {
                "activation_status", "reference_transcript", "reference_protein",
                "vcep_policy_id", "vcep_specification", "normalization_validation",
                "decision_assets", "thresholds",
                "functional_domains", "applicable_rules", "required_rule_data",
            },
            f"genes/{raw_gene}",
        )
        if gene["activation_status"] != "active":
            continue
        active_genes += 1
        if not _TRANSCRIPT_RE.fullmatch(str(gene["reference_transcript"])):
            raise GenePolicyError(f"Invalid reference transcript for {raw_gene}")
        if not _PROTEIN_RE.fullmatch(str(gene["reference_protein"])):
            raise GenePolicyError(f"Invalid reference protein for {raw_gene}")
        normalization_validation = gene["normalization_validation"]
        if not isinstance(normalization_validation, dict):
            raise GenePolicyError(
                f"Gene {raw_gene} normalization_validation must be an object"
            )
        _require_keys(
            normalization_validation,
            {"c_notation", "p_notation"},
            f"genes/{raw_gene}/normalization_validation",
        )
        if not str(normalization_validation["c_notation"]).startswith("c."):
            raise GenePolicyError(
                f"Gene {raw_gene} normalization smoke-test c. notation is invalid"
            )
        if not str(normalization_validation["p_notation"]).startswith("p."):
            raise GenePolicyError(
                f"Gene {raw_gene} normalization smoke-test p. notation is invalid"
            )
        decision_assets = gene["decision_assets"]
        if not isinstance(decision_assets, dict) or "PVS1" not in decision_assets:
            raise GenePolicyError(f"Gene {raw_gene} has no PVS1 decision assets")
        for branch in ("splice", "nonsplice"):
            asset = decision_assets["PVS1"].get(branch)
            if not isinstance(asset, dict):
                raise GenePolicyError(
                    f"Gene {raw_gene} has no PVS1 {branch} decision asset"
                )
            _require_keys(
                asset,
                {"figure_number", "figure_url"},
                f"genes/{raw_gene}/decision_assets/PVS1/{branch}",
            )
        policy_id = gene["vcep_policy_id"]
        if policy_id not in policies:
            raise GenePolicyError(f"Gene {raw_gene} refers to unknown policy {policy_id!r}")
        applicable = gene["applicable_rules"]
        if not isinstance(applicable, list) or not applicable or len(applicable) != len(set(applicable)):
            raise GenePolicyError(f"Gene {raw_gene} applicable_rules must be a unique non-empty list")
        unsupported = sorted(set(applicable) - set(policies[policy_id]["supported_rule_codes"]))
        if unsupported:
            raise GenePolicyError(f"Gene {raw_gene} enables rules not supported by its policy: {unsupported}")
        bayesdel = gene["thresholds"].get("bayesdel_noaf")
        if not isinstance(bayesdel, dict):
            raise GenePolicyError(f"Gene {raw_gene} has no BayesDel_noAF threshold set")
        _require_keys(
            bayesdel, {"bp4_max_inclusive", "pp3_min_inclusive"},
            f"genes/{raw_gene}/thresholds/bayesdel_noaf",
        )
        bp4 = _require_number(bayesdel["bp4_max_inclusive"], f"{raw_gene}/BayesDel/BP4")
        pp3 = _require_number(bayesdel["pp3_min_inclusive"], f"{raw_gene}/BayesDel/PP3")
        if not 0 <= bp4 < pp3 <= 1:
            raise GenePolicyError(f"Gene {raw_gene} BayesDel thresholds must satisfy 0 <= BP4 < PP3 <= 1")
        pvs1 = gene["thresholds"].get("pvs1")
        if not isinstance(pvs1, dict):
            raise GenePolicyError(f"Gene {raw_gene} has no PVS1 threshold set")
        _require_keys(
            pvs1, {"nmd_boundary_c_first_not_predicted"},
            f"genes/{raw_gene}/thresholds/pvs1",
        )
        nmd_boundary = _require_number(
            pvs1["nmd_boundary_c_first_not_predicted"],
            f"{raw_gene}/PVS1/NMD boundary",
        )
        if int(nmd_boundary) != nmd_boundary or nmd_boundary < 1:
            raise GenePolicyError(f"Gene {raw_gene} has an invalid PVS1 NMD boundary")
        domains = gene["functional_domains"]
        if not isinstance(domains, dict):
            raise GenePolicyError(f"Gene {raw_gene} functional_domains must be an object")
        for domain_name, interval in domains.items():
            _require_keys(
                interval,
                {"start", "end", "description"},
                f"genes/{raw_gene}/functional_domains/{domain_name}",
            )
            start = _require_number(interval["start"], f"{raw_gene}/{domain_name}/start")
            end = _require_number(interval["end"], f"{raw_gene}/{domain_name}/end")
            if int(start) != start or int(end) != end or start < 1 or end < start:
                raise GenePolicyError(f"Invalid functional domain interval for {raw_gene}/{domain_name}")
            if not str(interval["description"]).strip():
                raise GenePolicyError(
                    f"Functional domain description is missing for {raw_gene}/{domain_name}"
                )

    if metadata is not None:
        _require_keys(
            metadata,
            {
                "schema_version", "manifest_id", "manifest_version", "manifest_sha256",
                "active_gene_count", "active_policy_count", "validation_status",
            },
            "metadata",
        )
        if metadata["validation_status"] != "approved":
            raise GenePolicyError("Gene policy metadata is not approved")
        if metadata["manifest_id"] != manifest["manifest_id"] or metadata["manifest_version"] != manifest["manifest_version"]:
            raise GenePolicyError("Gene policy metadata identity differs from the manifest")
        if metadata["active_gene_count"] != active_genes or metadata["active_policy_count"] != active_policies:
            raise GenePolicyError("Gene policy metadata record counts differ from the manifest")
        if manifest_bytes is None:
            raise GenePolicyError("Raw manifest bytes are required for checksum validation")
        digest = hashlib.sha256(manifest_bytes).hexdigest()
        if metadata["manifest_sha256"] != digest:
            raise GenePolicyError("Gene policy manifest checksum mismatch")


@lru_cache(maxsize=1)
def load_gene_policy_manifest() -> dict[str, Any]:
    try:
        manifest_bytes = GENE_POLICY_MANIFEST_PATH.read_bytes()
    except OSError as exc:
        raise GenePolicyError(f"Required gene policy manifest cannot be read: {GENE_POLICY_MANIFEST_PATH}: {exc}") from exc
    manifest = _read_json(GENE_POLICY_MANIFEST_PATH, "gene policy manifest")
    metadata = _read_json(GENE_POLICY_METADATA_PATH, "gene policy metadata")
    validate_gene_policy_payload(manifest, metadata, manifest_bytes=manifest_bytes)
    return manifest


def active_genes() -> tuple[str, ...]:
    manifest = load_gene_policy_manifest()
    return tuple(
        gene for gene, record in manifest["genes"].items()
        if record["activation_status"] == "active"
    )


def get_gene_policy(gene: str) -> dict[str, Any]:
    symbol = str(gene).strip().upper()
    manifest = load_gene_policy_manifest()
    record = manifest["genes"].get(symbol)
    if not isinstance(record, dict) or record.get("activation_status") != "active":
        raise GenePolicyError(f"No active gene policy is configured for {symbol or '<empty>'}")
    policy = manifest["policies"].get(record["vcep_policy_id"])
    if not isinstance(policy, dict):
        raise GenePolicyError(f"Active gene {symbol} refers to an unavailable VCEP policy")
    return {"gene": symbol, "gene_config": deepcopy(record), "policy": deepcopy(policy)}


def reference_transcript(gene: str) -> str:
    return str(get_gene_policy(gene)["gene_config"]["reference_transcript"])


def vcep_specification(gene: str) -> dict[str, str]:
    return deepcopy(get_gene_policy(gene)["gene_config"]["vcep_specification"])


def normalization_validation_variant(gene: str) -> dict[str, str]:
    return deepcopy(get_gene_policy(gene)["gene_config"]["normalization_validation"])


def decision_asset(gene: str, criterion: str, branch: str) -> dict[str, str]:
    assets = get_gene_policy(gene)["gene_config"]["decision_assets"]
    try:
        return deepcopy(assets[criterion.upper()][branch])
    except (KeyError, TypeError) as exc:
        raise GenePolicyError(
            f"No {criterion.upper()} {branch} decision asset is configured for {gene}"
        ) from exc


def runtime_policy_id(gene: str) -> str:
    return str(get_gene_policy(gene)["policy"]["runtime_policy_id"])


def implementation_profile(gene: str) -> str:
    return str(get_gene_policy(gene)["policy"]["implementation_profile"])


def policy_version(gene: str) -> str:
    return str(get_gene_policy(gene)["policy"]["version"])


def policy_name(gene: str) -> str:
    return str(get_gene_policy(gene)["policy"]["name"])


def external_evidence_config(gene: str) -> dict[str, str]:
    return deepcopy(get_gene_policy(gene)["policy"]["external_evidence"])


def not_used_criteria(gene: str) -> list[dict[str, str]]:
    """Return policy-level ACMG/AMP uses explicitly rejected by the VCEP."""
    values = get_gene_policy(gene)["policy"]["criterion_applicability"]["not_used"]
    return deepcopy(values)


def resolve_policy_identity(gene: str | None = None) -> tuple[str, str]:
    """Resolve a policy without silently selecting one in a multi-policy panel."""
    if gene:
        return runtime_policy_id(gene), policy_version(gene)
    identities = {
        (runtime_policy_id(symbol), policy_version(symbol))
        for symbol in active_genes()
    }
    if len(identities) != 1:
        raise GenePolicyError(
            "Gene is required because active genes use different VCEP policies"
        )
    return next(iter(identities))


def resolve_policy_gene(gene: str | None = None) -> str:
    """Return a policy-equivalent gene only when that choice is unambiguous."""
    if gene:
        get_gene_policy(gene)
        return str(gene).strip().upper()
    resolve_policy_identity(None)
    return active_genes()[0]


def rule_is_applicable(gene: str, code: str) -> bool:
    normalized = str(code).strip().upper()
    applicable = set(get_gene_policy(gene)["gene_config"]["applicable_rules"])
    if normalized in applicable:
        return True
    for suffix in ("_VERY_STRONG", "_STRONG", "_MODERATE", "_SUPPORTING"):
        if normalized.endswith(suffix) and normalized[: -len(suffix)] in applicable:
            return True
    return False


def bayesdel_thresholds(gene: str) -> dict[str, float]:
    values = get_gene_policy(gene)["gene_config"]["thresholds"]["bayesdel_noaf"]
    return {
        "bp4": float(values["bp4_max_inclusive"]),
        "pp3": float(values["pp3_min_inclusive"]),
    }


def spliceai_thresholds(gene: str) -> dict[str, float]:
    values = get_gene_policy(gene)["policy"]["thresholds"]["spliceai"]
    return {
        "bp4": float(values["no_impact_max_inclusive"]),
        "pp3": float(values["impact_min_inclusive"]),
    }


def pvs1_thresholds(gene: str) -> dict[str, int]:
    values = get_gene_policy(gene)["gene_config"]["thresholds"]["pvs1"]
    return {
        "nmd_boundary_c_first_not_predicted": int(
            values["nmd_boundary_c_first_not_predicted"]
        )
    }


def functional_domains(gene: str) -> dict[str, tuple[int, int]]:
    values = get_gene_policy(gene)["gene_config"]["functional_domains"]
    return {
        name: (int(interval["start"]), int(interval["end"]))
        for name, interval in values.items()
    }


def functional_domain_descriptions(gene: str) -> dict[str, str]:
    values = get_gene_policy(gene)["gene_config"]["functional_domains"]
    return {
        name: str(interval["description"])
        for name, interval in values.items()
    }


def clinical_lr_thresholds(gene: str) -> dict[str, dict[str, float]]:
    values = get_gene_policy(gene)["policy"]["thresholds"]["clinical_likelihood_ratio"]
    return deepcopy(values)


def population_thresholds(gene: str) -> dict[str, Any]:
    values = get_gene_policy(gene)["policy"]["thresholds"]["population_frequency"]
    return deepcopy(values)


def mixed_evidence_point_thresholds(gene: str | None = None) -> dict[str, int]:
    policy_gene = resolve_policy_gene(gene)
    values = get_gene_policy(policy_gene)["policy"]["thresholds"]["mixed_evidence_points"]
    return {key: int(value) for key, value in values.items()}


def transcripts_by_gene() -> dict[str, str]:
    return {gene: reference_transcript(gene) for gene in active_genes()}


def domains_by_gene() -> dict[str, dict[str, tuple[int, int]]]:
    return {gene: functional_domains(gene) for gene in active_genes()}


def validate_policy_source_bindings(project_root: Path | None = None) -> None:
    """Prove that source-specific manifests agree with the runtime policy."""
    root = project_root or Path(__file__).resolve().parent.parent
    gnomad = _read_json(
        Path(__file__).resolve().parent / "data" / "gnomad" / "gnomad_panel_manifest.json",
        "gnomAD panel manifest",
    )
    spliceai = _read_json(
        root / "data" / "spliceai" / "enigma_v1_2_spliceai_profile.json",
        "SpliceAI scoring profile",
    )
    reference = _read_json(
        root / "data" / "reference" / "panel" / "panel_manifest.json",
        "reference transcript panel manifest",
    )
    reference_by_gene = {
        str(item.get("gene", "")).upper(): item
        for item in reference.get("transcripts", [])
    }
    gnomad_targets = {
        str(item.get("gene", "")).upper(): item
        for item in gnomad.get("targets", [])
        if item.get("activation_status") == "active"
    }

    for gene in active_genes():
        configured = get_gene_policy(gene)
        gene_config = configured["gene_config"]
        policy = configured["policy"]
        transcript = gene_config["reference_transcript"]
        protein = gene_config["reference_protein"]
        ref_record = reference_by_gene.get(gene)
        if not ref_record or ref_record.get("transcript") != transcript or ref_record.get("protein") != protein:
            raise GenePolicyError(f"Reference panel identity differs from the active policy for {gene}")
        splice_record = (spliceai.get("reference_transcripts") or {}).get(gene)
        if not splice_record or splice_record.get("refseq") != transcript:
            raise GenePolicyError(f"SpliceAI transcript differs from the active policy for {gene}")
        splice_threshold = policy["thresholds"]["spliceai"]
        profile_threshold = spliceai.get("thresholds") or {}
        if (
            profile_threshold.get("bp4_max_inclusive") != splice_threshold["no_impact_max_inclusive"]
            or profile_threshold.get("pp3_min_inclusive") != splice_threshold["impact_min_inclusive"]
        ):
            raise GenePolicyError(f"SpliceAI thresholds differ from the active policy for {gene}")

        target = gnomad_targets.get(gene)
        if not target or target.get("reference_transcript") != transcript:
            raise GenePolicyError(f"gnomAD target differs from the active policy for {gene}")
        if (target.get("vcep_specification") or {}).get("id") != gene_config["vcep_specification"]["id"]:
            raise GenePolicyError(f"gnomAD VCEP identity differs from the active policy for {gene}")
        frequency_id = policy["thresholds"]["population_frequency"]["gnomad_policy_id"]
        if target.get("classification_policy_id") != frequency_id:
            raise GenePolicyError(f"gnomAD policy binding differs from the active policy for {gene}")
        frequency = (gnomad.get("classification_policies") or {}).get(frequency_id, {}).get("frequency_criteria", {})
        top = policy["thresholds"]["population_frequency"]
        if (
            (frequency.get("ba1") or {}).get("threshold") != top["ba1_faf95_min_exclusive"]
            or ((frequency.get("bs1") or {}).get("strong") or {}).get("threshold") != top["bs1_strong_faf95_min_exclusive"]
            or ((frequency.get("bs1") or {}).get("supporting") or {}).get("lower_threshold") != top["bs1_supporting_faf95_min_exclusive"]
            or (frequency.get("ba1") or {}).get("minimum_mean_depth") != top["ba1_bs1_minimum_mean_depth"]
            or (frequency.get("pm2") or {}).get("minimum_mean_depth") != top["pm2_minimum_mean_depth"]
        ):
            raise GenePolicyError(f"gnomAD thresholds differ from the active policy for {gene}")
