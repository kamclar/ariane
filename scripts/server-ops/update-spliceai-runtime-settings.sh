#!/bin/bash

# Apply the versioned worker setting required by the shared SpliceAI gate.

set -euo pipefail

SERVICE_FILE="${ARIANE_SERVICE_FILE:-/etc/systemd/system/ariane.service}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if [ ! -f "$SERVICE_FILE" ]; then
    echo "ARIANE service definition is missing: $SERVICE_FILE" >&2
    exit 1
fi
if ! grep -q '^ExecStart=.*uvicorn backend\.main:app.*--workers [0-9][0-9]*' "$SERVICE_FILE"; then
    echo "The expected ARIANE uvicorn command was not found" >&2
    exit 1
fi

BACKUP_PATH="$(mktemp /tmp/ariane-service.XXXXXX)"
cp "$SERVICE_FILE" "$BACKUP_PATH"
sed -i 's/--workers [0-9][0-9]*/--workers 1/' "$SERVICE_FILE"

if ! grep -q '^ExecStart=.*--workers 1' "$SERVICE_FILE"; then
    cp "$BACKUP_PATH" "$SERVICE_FILE"
    rm -f "$BACKUP_PATH"
    echo "Could not configure the ARIANE worker count" >&2
    exit 1
fi

systemctl daemon-reload
rm -f "$BACKUP_PATH"
echo "ARIANE runtime worker setting updated"
