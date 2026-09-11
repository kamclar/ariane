from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_fresh_nginx_configuration_limits_classification_connections():
    script = (PROJECT_ROOT / "scripts/server-ops/deploy-ariane.sh").read_text(
        encoding="utf-8"
    )

    assert "limit_conn_status 429;" in script
    assert "limit_conn ariane_classification_ip_conn 4;" in script
    assert "limit_conn ariane_classification_key_conn 2;" in script
    assert "/(?:ui-api/classify|api/(?:v1/)?classify(?:/batch)?)" in script


def test_existing_nginx_update_adds_and_validates_connection_limits():
    script = (
        PROJECT_ROOT / "scripts/server-ops/update-nginx-api-settings.sh"
    ).read_text(encoding="utf-8")

    assert "limit_conn_status 429;" in script
    assert "limit_conn ariane_classification_ip_conn 4;" in script
    assert "limit_conn ariane_classification_key_conn 2;" in script
    assert "nginx -t" in script
    assert "Restoring the previous configuration" in script
