#!/bin/bash

# Restart ARIANE without stopping Nginx.

set -euo pipefail

ARIANE_PORT="${ARIANE_PORT:-8000}"
ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
ARIANE_USER="${ARIANE_USER:-ubuntu}"
ARIANE_ENV_FILE="${ARIANE_ENV_FILE:-/etc/ariane/ariane.env}"

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
if ! grep -q '^ARIANE_UI_SESSION_SECRET=' "$ARIANE_ENV_FILE"; then
    printf 'ARIANE_UI_SESSION_SECRET=%s\n' "$(openssl rand -hex 32)" >> "$ARIANE_ENV_FILE"
    chown root:"$ARIANE_USER" "$ARIANE_ENV_FILE"
    chmod 0640 "$ARIANE_ENV_FILE"
    echo "Created the ARIANE browser-session signing secret"
fi
if ! systemctl is-active --quiet ariane-spliceai.service; then
    echo "ARIANE local SpliceAI service is not running" >&2
    echo "Run: sudo bash $ARIANE_HOME/scripts/server-ops/install-spliceai-service.sh" >&2
    exit 1
fi
if ! python3 "$ARIANE_HOME/scripts/validate_spliceai_service.py" \
    --url http://127.0.0.1:8081/spliceai/ \
    --timeout 180; then
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
    'import cdot, fastapi, hgvs, pydantic, uvicorn; from backend.modules import hgvs_engine, hgvs_provider'

nginx -t
systemctl restart ariane
systemctl reload nginx

for attempt in $(seq 1 15); do
    if curl --fail --silent --max-time 5 "http://127.0.0.1:${ARIANE_PORT}/api/health" | grep -q '"status":"ok"'; then
        echo "ARIANE is healthy"
        exit 0
    fi
    echo "Waiting for ARIANE: $attempt/15"
    sleep 1
done

echo "ARIANE health check failed" >&2
journalctl -u ariane -n 30 --no-pager >&2
exit 1
