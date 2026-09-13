#!/bin/bash

# Apply the versioned API proxy settings to an existing ARIANE installation.

set -euo pipefail

NGINX_CONFIG="${NGINX_CONFIG:-/etc/nginx/sites-available/ariane}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if [ ! -f "$NGINX_CONFIG" ]; then
    echo "Nginx configuration is missing: $NGINX_CONFIG" >&2
    exit 1
fi

BACKUP_PATH="$(mktemp /tmp/ariane-nginx.XXXXXX)"
cp "$NGINX_CONFIG" "$BACKUP_PATH"

restore_config() {
    cp "$BACKUP_PATH" "$NGINX_CONFIG"
    rm -f "$BACKUP_PATH"
}

# Earlier installations limited the browser through the old public API path.
# The browser now uses a separate session-protected transport path.
sed -i 's|/api/classify $binary_remote_addr;|/ui-api/classify $binary_remote_addr;|' "$NGINX_CONFIG"

# Older configurations keyed this zone directly by client IP and applied it in
# location /. That also rate-limited the HTML and every static asset. A page
# with enough parallel CSS/JavaScript requests could consequently start only
# partially. Restrict the key to API paths; nginx ignores empty limit keys.
python3 - "$NGINX_CONFIG" <<'PY'
from pathlib import Path

path = Path(__import__("sys").argv[1])
text = path.read_text(encoding="utf-8")
old = "limit_req_zone $binary_remote_addr zone=ariane_api:10m rate=5r/s;"
old_scoped = "limit_req_zone $ariane_api_ip_key zone=ariane_api:10m rate=5r/s;"
new_zone = "limit_req_zone $ariane_api_ip_key zone=ariane_public_api:10m rate=5r/s;"
new = """map $uri $ariane_api_ip_key {
    default "";
    /api/health "";
    /api/v1/capabilities "";
    ~^/api/ $binary_remote_addr;
}
limit_req_zone $ariane_api_ip_key zone=ariane_public_api:10m rate=5r/s;"""
if old in text:
    text = text.replace(old, new, 1)
elif old_scoped in text:
    text = text.replace(old_scoped, new_zone, 1)
text = text.replace(
    "limit_req zone=ariane_api burst=20 nodelay;",
    "limit_req zone=ariane_public_api burst=20 nodelay;",
)
path.write_text(text, encoding="utf-8")
PY

if ! grep -q '^[[:space:]]*limit_req_status 429;' "$NGINX_CONFIG"; then
    sed -i '/^[[:space:]]*server_tokens off;[[:space:]]*$/a\    limit_req_status 429;' "$NGINX_CONFIG"
fi
if ! grep -q '^[[:space:]]*limit_conn_status 429;' "$NGINX_CONFIG"; then
    sed -i '/^[[:space:]]*limit_req_status 429;[[:space:]]*$/a\    limit_conn_status 429;' "$NGINX_CONFIG"
fi

if ! grep -q 'zone=ariane_web_classify:' "$NGINX_CONFIG"; then
    sed -i '/^limit_req_zone .* zone=ariane_public_api:/a\
map $uri $ariane_web_classify_key {\
    default "";\
    /ui-api/classify $binary_remote_addr;\
}\
limit_req_zone $ariane_web_classify_key zone=ariane_web_classify:10m rate=30r/m;' "$NGINX_CONFIG"
fi
if ! grep -q 'zone=ariane_keyed_api:' "$NGINX_CONFIG"; then
    sed -i '/^limit_req_zone .* zone=ariane_public_api:/a\limit_req_zone $http_x_ariane_api_key zone=ariane_keyed_api:10m rate=30r/m;' "$NGINX_CONFIG"
fi
if ! grep -q 'zone=ariane_classification_ip_conn:' "$NGINX_CONFIG"; then
    sed -i '/^limit_req_zone \$ariane_web_classify_key zone=ariane_web_classify:/a\
map $uri $ariane_classification_ip_conn_key {\
    default "";\
    ~^/(?:ui-api/classify|api/(?:v1/)?classify(?:/batch)?)$ $binary_remote_addr;\
}\
map $uri $ariane_classification_api_key_conn_key {\
    default "";\
    ~^/api/(?:v1/)?classify(?:/batch)?$ $http_x_ariane_api_key;\
}\
limit_conn_zone $ariane_classification_ip_conn_key zone=ariane_classification_ip_conn:10m;\
limit_conn_zone $ariane_classification_api_key_conn_key zone=ariane_classification_key_conn:10m;' "$NGINX_CONFIG"
fi
if ! grep -q 'limit_req zone=ariane_keyed_api burst=3 nodelay;' "$NGINX_CONFIG"; then
    sed -i '/limit_req zone=ariane_public_api burst=20 nodelay;/a\        limit_req zone=ariane_keyed_api burst=3 nodelay;' "$NGINX_CONFIG"
fi
if ! grep -q 'limit_req zone=ariane_web_classify burst=3 nodelay;' "$NGINX_CONFIG"; then
    sed -i '/limit_req zone=ariane_public_api burst=20 nodelay;/a\        limit_req zone=ariane_web_classify burst=3 nodelay;' "$NGINX_CONFIG"
fi
if ! grep -q 'limit_conn ariane_classification_ip_conn 4;' "$NGINX_CONFIG"; then
    sed -i '/limit_req zone=ariane_web_classify burst=3 nodelay;/a\        limit_conn ariane_classification_ip_conn 4;' "$NGINX_CONFIG"
fi
if ! grep -q 'limit_conn ariane_classification_key_conn 2;' "$NGINX_CONFIG"; then
    sed -i '/limit_conn ariane_classification_ip_conn 4;/a\        limit_conn ariane_classification_key_conn 2;' "$NGINX_CONFIG"
fi

sed -i 's/^[[:space:]]*proxy_read_timeout 60s;[[:space:]]*$/    proxy_read_timeout 180s;/' "$NGINX_CONFIG"

if ! grep -q '^[[:space:]]*limit_req_status 429;' "$NGINX_CONFIG"; then
    echo "Could not configure HTTP 429 for rate limiting" >&2
    restore_config
    exit 1
fi
if ! grep -q '^[[:space:]]*limit_conn_status 429;' "$NGINX_CONFIG"; then
    echo "Could not configure HTTP 429 for connection limiting" >&2
    restore_config
    exit 1
fi
if ! grep -q '^[[:space:]]*proxy_read_timeout 180s;' "$NGINX_CONFIG"; then
    echo "Could not configure the API response timeout" >&2
    restore_config
    exit 1
fi
if ! grep -q '^limit_req_zone \$ariane_api_ip_key zone=ariane_public_api:10m rate=5r/s;' "$NGINX_CONFIG"; then
    echo "Could not restrict the general rate limit to API paths" >&2
    restore_config
    exit 1
fi
if grep -q '^limit_req_zone \$binary_remote_addr zone=ariane_api:' "$NGINX_CONFIG"; then
    echo "The obsolete page-wide request limit is still present" >&2
    restore_config
    exit 1
fi
if grep -q 'limit_req zone=ariane_api ' "$NGINX_CONFIG"; then
    echo "The obsolete page-wide request limit is still applied" >&2
    restore_config
    exit 1
fi
if ! grep -q 'zone=ariane_web_classify:10m rate=30r/m;' "$NGINX_CONFIG"; then
    echo "Could not configure the interactive classification rate zone" >&2
    restore_config
    exit 1
fi
if ! grep -q 'zone=ariane_keyed_api:10m rate=30r/m;' "$NGINX_CONFIG"; then
    echo "Could not configure the per-key API rate zone" >&2
    restore_config
    exit 1
fi
if ! grep -q 'limit_req zone=ariane_keyed_api burst=3 nodelay;' "$NGINX_CONFIG"; then
    echo "Could not apply the per-key API rate limit" >&2
    restore_config
    exit 1
fi
if ! grep -q 'limit_req zone=ariane_web_classify burst=3 nodelay;' "$NGINX_CONFIG"; then
    echo "Could not apply the interactive classification rate limit" >&2
    restore_config
    exit 1
fi
if ! grep -q 'zone=ariane_classification_ip_conn:10m;' "$NGINX_CONFIG"; then
    echo "Could not configure the classification connection zone" >&2
    restore_config
    exit 1
fi
if ! grep -q 'limit_conn ariane_classification_ip_conn 4;' "$NGINX_CONFIG"; then
    echo "Could not apply the per-IP classification connection limit" >&2
    restore_config
    exit 1
fi
if ! grep -q 'limit_conn ariane_classification_key_conn 2;' "$NGINX_CONFIG"; then
    echo "Could not apply the per-key classification connection limit" >&2
    restore_config
    exit 1
fi
if ! nginx -t; then
    echo "Nginx validation failed. Restoring the previous configuration." >&2
    restore_config
    nginx -t
    exit 1
fi

systemctl reload nginx
rm -f "$BACKUP_PATH"
echo "ARIANE API proxy settings updated"
