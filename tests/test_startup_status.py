from __future__ import annotations

from backend.bootstrap import StartupStatus


def test_required_startup_failure_is_structured_and_not_raised() -> None:
    status = StartupStatus()

    result = status.run(
        "ENIGMA Table 9",
        lambda: (_ for _ in ()).throw(RuntimeError("checksum mismatch")),
        required=True,
    )

    assert result is None
    assert status.ready is False
    assert status.report()["failures"] == [
        {
            "component": "ENIGMA Table 9",
            "required": True,
            "status": "failed",
            "detail": "RuntimeError: checksum mismatch",
        }
    ]


def test_optional_startup_failure_degrades_without_blocking_readiness() -> None:
    status = StartupStatus()
    status.run("required data", lambda: object(), required=True)
    status.run(
        "runtime cache",
        lambda: (_ for _ in ()).throw(OSError("read failed")),
        required=False,
    )

    assert status.ready is True
    assert status.degraded is True
