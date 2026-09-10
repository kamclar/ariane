from pathlib import Path
import sqlite3

from fastapi import FastAPI
from fastapi.testclient import TestClient

from backend.review_api import router
from backend.review_records import ReviewRecordRepository


def _repository_draft(repository: ReviewRecordRepository):
    return repository.create_draft(
        variant_context={
            "gene": "BRCA1",
            "c_notation": "c.5366C>T",
            "p_notation": "p.(Ala1789Val)",
        },
        module1_result={"predicted_class": 3, "criteria": []},
        manual_evidence={"criteria": [{"code": "PP1"}]},
        amended_result={"predicted_class": 4, "criteria": [{"code": "PP1"}]},
        reviewer_name="Test reviewer",
        reviewer_role="Variant curator",
        review_date="2026-09-02",
        authenticated_account="test-admin",
    )


def test_review_repository_creates_immutable_versions(tmp_path: Path):
    repository = ReviewRecordRepository(tmp_path / "reviews.sqlite3")
    draft = _repository_draft(repository)
    assert draft["status"] == "draft"
    assert draft["version"] == 1
    assert len(draft["hashes"]["content_sha256"]) == 64

    approved = repository.approve(
        draft["record_id"],
        approver_name="Approving expert",
        approver_role="Clinical geneticist",
        approval_comment="Evidence and interactions reviewed.",
        authenticated_account="test-admin",
    )
    assert approved["status"] == "approved"
    assert approved["version"] == 2
    assert approved["parent_record_id"] == draft["record_id"]
    assert approved["hashes"]["content_sha256"] == draft["hashes"]["content_sha256"]
    assert repository.get(draft["record_id"])["status"] == "draft"

    next_draft = _repository_draft(repository)
    assert next_draft["version"] == 3
    assert next_draft["parent_record_id"] == approved["record_id"]


def test_approved_record_cannot_be_approved_again(tmp_path: Path):
    repository = ReviewRecordRepository(tmp_path / "reviews.sqlite3")
    draft = _repository_draft(repository)
    approved = repository.approve(
        draft["record_id"],
        approver_name="Approving expert",
        approver_role="Clinical geneticist",
        approval_comment="",
        authenticated_account="test-admin",
    )
    try:
        repository.approve(
            approved["record_id"],
            approver_name="Approving expert",
            approver_role="Clinical geneticist",
            approval_comment="",
            authenticated_account="test-admin",
        )
    except ValueError as exc:
        assert "Only a draft" in str(exc)
    else:
        raise AssertionError("Approving an approved immutable version must fail")


def test_record_read_fails_closed_after_database_tampering(tmp_path: Path):
    database = tmp_path / "reviews.sqlite3"
    repository = ReviewRecordRepository(database)
    draft = _repository_draft(repository)
    with sqlite3.connect(database) as connection:
        connection.execute(
            "UPDATE review_records SET amended_result_json = ? WHERE record_id = ?",
            ('{"predicted_class":1}', draft["record_id"]),
        )
        connection.commit()
    try:
        repository.get(draft["record_id"])
    except RuntimeError as exc:
        assert "failed amended_result_sha256 verification" in str(exc)
    else:
        raise AssertionError("A modified review record must fail integrity verification")


def _api_payload():
    return {
        "module1_result": {
            "variant": "BRCA1 c.5366C>T p.(Ala1789Val)",
            "gene": "BRCA1",
            "c_notation": "c.5366C>T",
            "p_notation": "p.(Ala1789Val)",
            "predicted_class": 3,
            "predicted_label": "VUS",
            "total_points": 0,
            "criteria": [],
        },
        "manual_evidence": {
            "base_criteria": [],
            "variant_context": {
                "gene": "BRCA1",
                "c_notation": "c.5366C>T",
                "p_notation": "p.(Ala1789Val)",
            },
            "manual_criteria": [{
                "code": "PP1",
                "enabled": True,
                "evidence": {"likelihood_ratio": 18.7},
                "notes": "Quantitative segregation LR reviewed.",
                "references": ["PMID:12345678"],
            }],
            "assessor": "Test reviewer",
            "assessed_at": "2026-09-02",
        },
        "reviewer_role": "Variant curator",
    }


def test_review_api_requires_authentication_and_recomputes_result(tmp_path, monkeypatch):
    monkeypatch.setenv("ARIANE_ADMIN_USER", "review-admin")
    monkeypatch.setenv("ARIANE_ADMIN_PASSWORD", "test-password")
    monkeypatch.setenv("ARIANE_RUNTIME_DATA_DIR", str(tmp_path))
    app = FastAPI()
    app.include_router(router)
    client = TestClient(app)

    assert client.post("/api/review-records", json=_api_payload()).status_code == 401
    response = client.post(
        "/api/review-records",
        json=_api_payload(),
        auth=("review-admin", "test-password"),
    )
    assert response.status_code == 200, response.text
    draft = response.json()
    assert draft["status"] == "draft"
    assert draft["amended_result"]["predicted_class"] == 3
    assert draft["amended_result"]["total_points"] == 4

    approval = client.post(
        f"/api/review-records/{draft['record_id']}/approve",
        json={
            "approver_name": "Approving expert",
            "approver_role": "Clinical geneticist",
            "approval_comment": "Reviewed.",
            "attestation_confirmed": True,
        },
        auth=("review-admin", "test-password"),
    )
    assert approval.status_code == 200, approval.text
    assert approval.json()["status"] == "approved"
