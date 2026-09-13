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


def test_fresh_nginx_configuration_does_not_rate_limit_frontend_assets():
    script = (PROJECT_ROOT / "scripts/server-ops/deploy-ariane.sh").read_text(
        encoding="utf-8"
    )

    assert "map $uri $ariane_api_ip_key" in script
    assert "~^/api/ $binary_remote_addr;" in script
    assert "limit_req_zone $ariane_api_ip_key zone=ariane_public_api" in script
    assert "limit_req_zone $binary_remote_addr zone=ariane_api" not in script
    assert "limit_req zone=ariane_public_api burst=20 nodelay;" in script
    assert "limit_req zone=ariane_api burst=20 nodelay;" not in script


def test_existing_nginx_update_adds_and_validates_connection_limits():
    script = (
        PROJECT_ROOT / "scripts/server-ops/update-nginx-api-settings.sh"
    ).read_text(encoding="utf-8")

    assert "limit_conn_status 429;" in script
    assert "limit_conn ariane_classification_ip_conn 4;" in script
    assert "limit_conn ariane_classification_key_conn 2;" in script
    assert "nginx -t" in script
    assert "Restoring the previous configuration" in script


def test_existing_nginx_update_removes_page_wide_request_limit():
    script = (
        PROJECT_ROOT / "scripts/server-ops/update-nginx-api-settings.sh"
    ).read_text(encoding="utf-8")

    assert 'old = "limit_req_zone $binary_remote_addr zone=ariane_api:10m rate=5r/s;"' in script
    assert "limit_req_zone $ariane_api_ip_key zone=ariane_public_api:10m rate=5r/s;" in script
    assert "The obsolete page-wide request limit is still present" in script
