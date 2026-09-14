import importlib.util
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from backend.lookups import spliceai


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "validate_spliceai_service.py"
SPEC = importlib.util.spec_from_file_location("validate_spliceai_service", SCRIPT_PATH)
assert SPEC and SPEC.loader
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


def load_inputs():
    profile = validator.load_json_object(validator.PROFILE_PATH)
    fixture = validator.load_json_object(validator.CASES_PATH)
    return profile, fixture["cases"][0]


def valid_payload(profile, case):
    row = {
        "t_id": case["transcript"],
        "t_refseq_ids": [case["refseq"]],
        **case["delta_scores"],
        **case["reference_scores"],
        **case["alternate_scores"],
    }
    return {
        "hg": "38",
        "genomeVersion": "38",
        "bc": profile["annotation_subset"],
        "distance": profile["max_distance"],
        "mask": profile["mask"],
        "scores": [row],
    }


def test_validator_accepts_exact_profile_and_reference_transcript_scores():
    profile, case = load_inputs()
    result = validator.validate_payload(valid_payload(profile, case), profile, case)
    assert result == {
        "id": "BRCA1_c.4185G>A",
        "variant": "chr17-43090944-C-T",
        "transcript": "ENST00000357654.9",
        "score": 0.93,
    }


def test_validator_rejects_score_difference():
    profile, case = load_inputs()
    payload = valid_payload(profile, case)
    payload["scores"][0]["DS_DL"] = 0.929
    with pytest.raises(validator.ValidationError, match="Score mismatch"):
        validator.validate_payload(payload, profile, case)


def test_validator_rejects_wrong_transcript_version():
    profile, case = load_inputs()
    payload = valid_payload(profile, case)
    payload["scores"][0]["t_id"] = "ENST00000357654.8"
    with pytest.raises(validator.ValidationError, match="Expected one row"):
        validator.validate_payload(payload, profile, case)


def test_service_validation_uses_versioned_cases_without_network():
    profile = validator.load_json_object(validator.PROFILE_PATH)
    fixture = validator.load_json_object(validator.CASES_PATH)
    payloads = [valid_payload(profile, case) for case in fixture["cases"]]
    with patch.object(validator, "fetch_json", side_effect=payloads) as fetch:
        result = validator.validate_service("http://127.0.0.1:8082/spliceai/")
    assert result["status"] == "ok"
    assert result["profile_id"].endswith("-v3")
    assert "@sha256:" in result["docker_image"]
    assert len(result["cases"]) == 2
    requested_url = fetch.call_args_list[0].args[0]
    assert "distance=10000" in requested_url
    assert "mask=0" in requested_url
    assert "bc=basic" in requested_url
    assert "show-ref-alt=1" in requested_url


def test_validator_rejects_additional_transcript_rows():
    profile, case = load_inputs()
    payload = valid_payload(profile, case)
    payload["scores"].append({**payload["scores"][0], "t_id": "ENST_OTHER.1"})
    with pytest.raises(validator.ValidationError, match="more than one transcript"):
        validator.validate_payload(payload, profile, case)


def test_reference_transcript_annotation_matches_profile_and_checksum():
    profile = validator.load_json_object(validator.PROFILE_PATH)
    engine = profile["approved_engine"]
    annotation_path = PROJECT_ROOT / engine["reference_transcript_annotation"]
    assert hashlib.sha256(annotation_path.read_bytes()).hexdigest() == engine[
        "reference_transcript_annotation_sha256"
    ]
    transcript_ids = {
        line.split("\t", 1)[0]
        for line in annotation_path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }
    assert transcript_ids == {
        transcript["ensembl"]
        for transcript in profile["reference_transcripts"].values()
    }


def test_concurrent_validation_checks_every_configured_worker_slot():
    expected = {
        "status": "ok",
        "profile_id": "profile",
        "docker_image": "image",
        "cases": [],
    }
    with patch.object(validator, "validate_service", return_value=expected) as validate:
        result = validator.validate_service_copies(
            "http://127.0.0.1:8082/spliceai/",
            copies=3,
        )
    assert validate.call_count == 3
    assert result["concurrent_validation_copies"] == 3


def test_local_service_image_and_operational_scripts_are_pinned():
    profile = json.loads(validator.PROFILE_PATH.read_text(encoding="utf-8"))
    image = profile["approved_engine"]["docker_image"]
    assert image == (
        "docker.io/weisburd/spliceai-38@sha256:"
        "c2bb3d5c65eee01087b5dd130d72233fc812d158dce9779ae63117e826ec8399"
    )
    installer = (
        PROJECT_ROOT / "scripts" / "server-ops" / "install-spliceai-service.sh"
    ).read_text(encoding="utf-8")
    assert "spliceai-38:latest" not in installer
    assert "127.0.0.1:${SPLICEAI_PORT}:8080" in installer
    assert 'SPLICEAI_PORT="${SPLICEAI_PORT:-8082}"' in installer
    assert "SpliceAI port is already used by another service" in installer
    assert "restore_service()" in installer
    assert "DATABASE_ENABLED=0" in installer
    assert "--env WORKERS=${SPLICEAI_WORKERS}" in installer
    assert "reference_transcript_annotation_sha256" in installer
    assert "dst=/gencode.v49.basic.annotation.txt.gz,readonly" in installer
    assert "ExecStartPost=" in installer
    assert "wait-and-validate-spliceai-service.sh" in installer
    assert "python3 \"$VALIDATOR\"" in installer
    restart = (
        PROJECT_ROOT / "scripts" / "server-ops" / "restart-ariane.sh"
    ).read_text(encoding="utf-8")
    assert "validate_spliceai_service.py" in restart
    assert "systemctl is-active --quiet ariane-spliceai.service" in restart
    assert "ARIANE_RUNTIME_DATA_DIR" in restart
    assert "--concurrent-copies 3" in restart
    assert "--workers 1" in restart
    assert "install-ariane-service.sh" in restart
    assert "$(seq 1 60)" in restart
    assert "systemctl status ariane --no-pager -l" in restart
    assert "journalctl -u ariane -n 50 --no-pager" in restart

    for script_name in ("install-ariane-service.sh", "deploy-ariane.sh"):
        service_script = (
            PROJECT_ROOT / "scripts" / "server-ops" / script_name
        ).read_text(encoding="utf-8")
        assert "ReadWritePaths=$ARIANE_HOME/backend/data" not in service_script
        assert "ReadWritePaths=/var/lib/ariane/runtime-cache" in service_script
        assert "$(seq 1 60)" in service_script
        assert "health_ok=1" in service_script
        assert "systemctl status ariane --no-pager -l" in service_script
        assert "journalctl -u ariane -n 50 --no-pager" in service_script

    assert "must keep versioned backend/data reference files read-only" in restart

    startup_validator = (
        PROJECT_ROOT
        / "scripts"
        / "server-ops"
        / "wait-and-validate-spliceai-service.sh"
    ).read_text(encoding="utf-8")
    assert '--concurrent-copies "$SPLICEAI_WORKERS"' in startup_validator
    assert "http://127.0.0.1:${SPLICEAI_PORT}/" in startup_validator
    assert "Installed reference-transcript SpliceAI annotation checksum" in startup_validator


def test_local_runtime_health_checks_only_the_configured_private_service():
    connection = patch.object(spliceai.socket, "create_connection")
    with patch.object(
        spliceai, "SPLICEAI_API_URL", "http://127.0.0.1:8082/spliceai/"
    ), connection as create_connection:
        create_connection.return_value.__enter__.return_value = object()
        health = spliceai.spliceai_runtime_health()
    assert health["status"] == "ok"
    assert health["local"] is True
    create_connection.assert_called_once_with(("127.0.0.1", 8082), timeout=0.25)


def test_local_runtime_health_reports_unavailable_service():
    with patch.object(
        spliceai, "SPLICEAI_API_URL", "http://127.0.0.1:8082/spliceai/"
    ), patch.object(
        spliceai.socket,
        "create_connection",
        side_effect=ConnectionRefusedError("not listening"),
    ):
        health = spliceai.spliceai_runtime_health()
    assert health["status"] == "unavailable"
    assert "ConnectionRefusedError" in health["reason"]
