#!/bin/bash

# Apply the versioned worker setting required by the shared SpliceAI gate.

set -euo pipefail

SERVICE_FILE="${ARIANE_SERVICE_FILE:-/etc/systemd/system/ariane.service}"
ENV_FILE="${ARIANE_ENV_FILE:-/etc/ariane/ariane.env}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if [ ! -f "$SERVICE_FILE" ]; then
    echo "ARIANE service definition is missing: $SERVICE_FILE" >&2
    exit 1
fi
if [ ! -f "$ENV_FILE" ]; then
    echo "ARIANE environment file is missing: $ENV_FILE" >&2
    exit 1
fi
if ! grep -q '^ExecStart=.*uvicorn backend\.main:app.*--workers [0-9][0-9]*' "$SERVICE_FILE"; then
    echo "The expected ARIANE uvicorn command was not found" >&2
    exit 1
fi

BACKUP_PATH="$(mktemp /tmp/ariane-service.XXXXXX)"
ENV_BACKUP_PATH="$(mktemp /tmp/ariane-env.XXXXXX)"
cp "$SERVICE_FILE" "$BACKUP_PATH"
cp "$ENV_FILE" "$ENV_BACKUP_PATH"

restore_settings() {
    cp "$BACKUP_PATH" "$SERVICE_FILE"
    cp "$ENV_BACKUP_PATH" "$ENV_FILE"
    rm -f "$BACKUP_PATH" "$ENV_BACKUP_PATH"
}

trap restore_settings ERR
sed -i 's/--workers [0-9][0-9]*/--workers 1/' "$SERVICE_FILE"

set_value() {
    local name="$1"
    local value="$2"
    if grep -q "^${name}=" "$ENV_FILE"; then
        sed -i "s|^${name}=.*$|${name}=${value}|" "$ENV_FILE"
    else
        printf '%s=%s\n' "$name" "$value" >> "$ENV_FILE"
    fi
}

set_value SPLICEAI_API_TIMEOUT 20
set_value SPLICEAI_LOOKUP_TIMEOUT 55
set_value SPLICEAI_API_ATTEMPTS 2
set_value SPLICEAI_API_RETRY_DELAY 2
set_value SPLICEAI_API_MAX_CONCURRENT 2
set_value SPLICEAI_API_RATE_SLEEP 1.5

if ! grep -q '^ExecStart=.*--workers 1' "$SERVICE_FILE"; then
    echo "Could not configure the ARIANE worker count" >&2
    false
fi

systemctl daemon-reload
trap - ERR
rm -f "$BACKUP_PATH" "$ENV_BACKUP_PATH"
echo "ARIANE SpliceAI runtime settings updated"
