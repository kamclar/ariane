#!/bin/bash

# Restart ARIANE without stopping Nginx.

set -euo pipefail

ARIANE_PORT="${ARIANE_PORT:-8000}"
ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
ARIANE_USER="${ARIANE_USER:-ubuntu}"
ARIANE_ENV_FILE="${ARIANE_ENV_FILE:-/etc/ariane/ariane.env}"
ARIANE_RUNTIME_CACHE_DIR="${ARIANE_RUNTIME_CACHE_DIR:-/var/lib/ariane/runtime-cache}"
ARIANE_RUNTIME_DATA_DIR="${ARIANE_RUNTIME_DATA_DIR:-/var/lib/ariane/runtime-data}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi

if [ ! -x "$ARIANE_HOME/venv/bin/python" ]; then
    echo "ARIANE virtual environment is missing: $ARIANE_HOME/venv" >&2
    exit 1
fi
if [ ! -f "$ARIANE_HOME/requirements.txt" ]; then
    echo "ARIANE requirements file is missing: $ARIANE_HOME/requirements.txt" >&2
    exit 1
fi
if ! command -v pg_config > /dev/null 2>&1; then
    echo "Missing pg_config required to install the pinned hgvs dependency." >&2
    echo "Install it with: apt-get update && apt-get install -y libpq-dev python3-dev build-essential" >&2
    exit 1
fi
if [ ! -f "$ARIANE_ENV_FILE" ]; then
    echo "ARIANE environment file is missing: $ARIANE_ENV_FILE" >&2
    exit 1
fi
install -d -m 0750 -o "$ARIANE_USER" -g "$ARIANE_USER" "$ARIANE_RUNTIME_CACHE_DIR"
install -d -m 0750 -o "$ARIANE_USER" -g "$ARIANE_USER" "$ARIANE_RUNTIME_DATA_DIR"

set_value() {
    local name="$1"
    local value="$2"
    if grep -q "^${name}=" "$ARIANE_ENV_FILE"; then
        sed -i "s|^${name}=.*$|${name}=${value}|" "$ARIANE_ENV_FILE"
    else
        printf '%s=%s\n' "$name" "$value" >> "$ARIANE_ENV_FILE"
    fi
}

set_value ARIANE_RUNTIME_CACHE_DIR "$ARIANE_RUNTIME_CACHE_DIR"
set_value ARIANE_RUNTIME_DATA_DIR "$ARIANE_RUNTIME_DATA_DIR"
chown root:"$ARIANE_USER" "$ARIANE_ENV_FILE"
chmod 0640 "$ARIANE_ENV_FILE"

SERVICE_DEFINITION="$(systemctl cat ariane.service 2>/dev/null || true)"
if ! grep -Fq "ReadWritePaths=/var/lib/ariane/runtime-cache /var/lib/ariane/runtime-data" <<<"$SERVICE_DEFINITION"; then
    echo "ARIANE service does not expose both persistent runtime directories" >&2
    echo "Run: sudo bash $ARIANE_HOME/scripts/server-ops/install-ariane-service.sh" >&2
    exit 1
fi
if grep -Eq '^ReadWritePaths=.*backend/data' <<<"$SERVICE_DEFINITION"; then
    echo "ARIANE service must keep versioned backend/data reference files read-only" >&2
    echo "Run: sudo bash $ARIANE_HOME/scripts/server-ops/install-ariane-service.sh" >&2
    exit 1
fi
if ! grep -Eq 'ExecStart=.*uvicorn .*--workers 1([[:space:]]|$)' <<<"$SERVICE_DEFINITION"; then
    echo "ARIANE service must use one application worker with the shared runtime caches" >&2
    echo "Run: sudo bash $ARIANE_HOME/scripts/server-ops/install-ariane-service.sh" >&2
    exit 1
fi
if ! grep -q '^ARIANE_UI_SESSION_SECRET=' "$ARIANE_ENV_FILE"; then
    printf 'ARIANE_UI_SESSION_SECRET=%s\n' "$(openssl rand -hex 32)" >> "$ARIANE_ENV_FILE"
    chown root:"$ARIANE_USER" "$ARIANE_ENV_FILE"
    chmod 0640 "$ARIANE_ENV_FILE"
    echo "Created the ARIANE browser-session signing secret"
fi
if ! grep -q '^ARIANE_API_DAILY_CLASSIFICATION_LIMIT=' "$ARIANE_ENV_FILE"; then
    printf 'ARIANE_API_DAILY_CLASSIFICATION_LIMIT=5000\n' >> "$ARIANE_ENV_FILE"
    chown root:"$ARIANE_USER" "$ARIANE_ENV_FILE"
    chmod 0640 "$ARIANE_ENV_FILE"
    echo "Configured the public API daily classification limit"
fi
if ! systemctl is-active --quiet ariane-spliceai.service; then
    echo "ARIANE local SpliceAI service is not running" >&2
    echo "Run: sudo bash $ARIANE_HOME/scripts/server-ops/install-spliceai-service.sh" >&2
    exit 1
fi
SPLICEAI_SERVICE_URL="$(sed -n 's/^SPLICEAI_API_URL=//p' "$ARIANE_ENV_FILE" | tail -n 1)"
SPLICEAI_SERVICE_URL="${SPLICEAI_SERVICE_URL:-http://127.0.0.1:8082/spliceai/}"
if ! [[ "$SPLICEAI_SERVICE_URL" =~ ^http://127\.0\.0\.1:[0-9]+/spliceai/$ ]]; then
    echo "ARIANE has an invalid private SpliceAI URL: $SPLICEAI_SERVICE_URL" >&2
    exit 1
fi
if ! python3 "$ARIANE_HOME/scripts/validate_spliceai_service.py" \
    --url "$SPLICEAI_SERVICE_URL" \
    --timeout 180 \
    --concurrent-copies 3; then
    echo "The running SpliceAI service does not match the active ARIANE profile" >&2
    echo "Run: sudo bash $ARIANE_HOME/scripts/server-ops/install-spliceai-service.sh" >&2
    exit 1
fi

# A code update may add or change a pinned runtime dependency. Install through
# the exact interpreter used by systemd, under the service account, before
# restarting. If dependency synchronization or the import preflight fails, the
# currently running service is left untouched.
echo "Synchronizing ARIANE runtime dependencies"
runuser -u "$ARIANE_USER" -- \
    "$ARIANE_HOME/venv/bin/python" -m pip install \
    --disable-pip-version-check \
    -r "$ARIANE_HOME/requirements.txt"

runuser -u "$ARIANE_USER" -- \
    env PYTHONPATH="$ARIANE_HOME" \
    "$ARIANE_HOME/venv/bin/python" -c \
    'import cdot, fastapi, hgvs, pydantic, uvicorn; from backend.variant_processing import hgvs_engine, hgvs_provider'

nginx -t
systemctl restart ariane
systemctl reload nginx

for attempt in $(seq 1 60); do
    if curl --fail --silent --max-time 5 "http://127.0.0.1:${ARIANE_PORT}/api/health" | grep -q '"status":"ok"'; then
        echo "ARIANE is healthy"
        exit 0
    fi
    if ! systemctl is-active --quiet ariane; then
        echo "ARIANE stopped while waiting for the health endpoint" >&2
        break
    fi
    echo "Waiting for ARIANE startup: $attempt/60"
    sleep 1
done

echo "ARIANE health check failed after waiting up to 60 seconds" >&2
systemctl status ariane --no-pager -l >&2 || true
journalctl -u ariane -n 50 --no-pager >&2 || true
exit 1
