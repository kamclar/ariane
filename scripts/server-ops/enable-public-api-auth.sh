#!/bin/bash

# Configure checksum-only API-key storage for an existing ARIANE service.

set -euo pipefail

ARIANE_USER="${ARIANE_USER:-ubuntu}"
ENV_FILE="${ARIANE_ENV_FILE:-/etc/ariane/ariane.env}"
KEY_FILE="${ARIANE_API_KEYS_FILE:-/etc/ariane/api-keys.json}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if ! id "$ARIANE_USER" >/dev/null 2>&1; then
    echo "ARIANE service user does not exist: $ARIANE_USER" >&2
    exit 1
fi
if [ ! -f "$ENV_FILE" ]; then
    echo "ARIANE environment file is missing: $ENV_FILE" >&2
    exit 1
fi

install -d -m 0750 -o root -g "$ARIANE_USER" "$(dirname "$KEY_FILE")"
if [ ! -f "$KEY_FILE" ]; then
    printf '{"schema_version":1,"keys":[]}\n' > "$KEY_FILE"
fi
chown root:"$ARIANE_USER" "$KEY_FILE"
chmod 0640 "$KEY_FILE"

if grep -q '^ARIANE_API_KEYS_FILE=' "$ENV_FILE"; then
    sed -i "s|^ARIANE_API_KEYS_FILE=.*$|ARIANE_API_KEYS_FILE=$KEY_FILE|" "$ENV_FILE"
else
    printf 'ARIANE_API_KEYS_FILE=%s\n' "$KEY_FILE" >> "$ENV_FILE"
fi
if ! grep -q '^ARIANE_UI_SESSION_SECRET=' "$ENV_FILE"; then
    printf 'ARIANE_UI_SESSION_SECRET=%s\n' "$(openssl rand -hex 32)" >> "$ENV_FILE"
fi
if ! grep -q '^ARIANE_API_DAILY_CLASSIFICATION_LIMIT=' "$ENV_FILE"; then
    printf 'ARIANE_API_DAILY_CLASSIFICATION_LIMIT=5000\n' >> "$ENV_FILE"
fi
chown root:"$ARIANE_USER" "$ENV_FILE"
chmod 0640 "$ENV_FILE"

echo "ARIANE public API authentication configured"
echo "Registry: $KEY_FILE"
echo "Create at least one key before restarting the service"
