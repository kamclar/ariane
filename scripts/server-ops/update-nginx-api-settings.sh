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

if ! grep -q '^[[:space:]]*limit_req_status 429;' "$NGINX_CONFIG"; then
    sed -i '/^[[:space:]]*server_tokens off;[[:space:]]*$/a\    limit_req_status 429;' "$NGINX_CONFIG"
fi

if ! grep -q 'zone=ariane_web_classify:' "$NGINX_CONFIG"; then
    sed -i '/^limit_req_zone \$binary_remote_addr zone=ariane_api:/a\
map $uri $ariane_web_classify_key {\
    default "";\
    /api/classify $binary_remote_addr;\
}\
limit_req_zone $ariane_web_classify_key zone=ariane_web_classify:10m rate=30r/m;' "$NGINX_CONFIG"
fi
if ! grep -q 'zone=ariane_keyed_api:' "$NGINX_CONFIG"; then
    sed -i '/^limit_req_zone \$binary_remote_addr zone=ariane_api:/a\limit_req_zone $http_x_ariane_api_key zone=ariane_keyed_api:10m rate=30r/m;' "$NGINX_CONFIG"
fi
if ! grep -q 'limit_req zone=ariane_keyed_api burst=3 nodelay;' "$NGINX_CONFIG"; then
    sed -i '/limit_req zone=ariane_api burst=20 nodelay;/a\        limit_req zone=ariane_keyed_api burst=3 nodelay;' "$NGINX_CONFIG"
fi
if ! grep -q 'limit_req zone=ariane_web_classify burst=3 nodelay;' "$NGINX_CONFIG"; then
    sed -i '/limit_req zone=ariane_api burst=20 nodelay;/a\        limit_req zone=ariane_web_classify burst=3 nodelay;' "$NGINX_CONFIG"
fi

sed -i 's/^[[:space:]]*proxy_read_timeout 60s;[[:space:]]*$/    proxy_read_timeout 180s;/' "$NGINX_CONFIG"

if ! grep -q '^[[:space:]]*limit_req_status 429;' "$NGINX_CONFIG"; then
    echo "Could not configure HTTP 429 for rate limiting" >&2
    restore_config
    exit 1
fi
if ! grep -q '^[[:space:]]*proxy_read_timeout 180s;' "$NGINX_CONFIG"; then
    echo "Could not configure the API response timeout" >&2
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
if ! nginx -t; then
    echo "Nginx validation failed. Restoring the previous configuration." >&2
    restore_config
    nginx -t
    exit 1
fi

systemctl reload nginx
rm -f "$BACKUP_PATH"
echo "ARIANE API proxy settings updated"
