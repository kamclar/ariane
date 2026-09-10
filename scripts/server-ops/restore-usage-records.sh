#!/bin/bash

# Restore classification usage statistics after checksum and SQLite validation.

set -euo pipefail

ARIANE_USER="${ARIANE_USER:-ubuntu}"
ARIANE_RUNTIME_DATA_DIR="${ARIANE_RUNTIME_DATA_DIR:-/var/lib/ariane/runtime-data}"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <ariane-usage-*.sqlite3.gz>" >&2
    exit 2
fi

BACKUP_FILE="$1"
CHECKSUM_FILE="${BACKUP_FILE}.sha256"
if [ ! -f "$BACKUP_FILE" ] || [ ! -f "$CHECKSUM_FILE" ]; then
    echo "Backup or checksum file not found" >&2
    exit 1
fi

(cd "$(dirname "$BACKUP_FILE")" && sha256sum --check "$(basename "$CHECKSUM_FILE")")
gzip -t "$BACKUP_FILE"
install -d -m 0750 -o "$ARIANE_USER" -g "$ARIANE_USER" "$ARIANE_RUNTIME_DATA_DIR"
STAGING_FILE="$(mktemp "$ARIANE_RUNTIME_DATA_DIR/.classification-usage.restore.XXXXXX")"
trap 'rm -f -- "$STAGING_FILE"' EXIT
gzip -cd "$BACKUP_FILE" > "$STAGING_FILE"
python3 - "$STAGING_FILE" <<'PY'
import sqlite3
import sys

with sqlite3.connect(sys.argv[1]) as connection:
    result = connection.execute("PRAGMA integrity_check").fetchone()[0]
    if result != "ok":
        raise SystemExit(f"SQLite integrity check failed: {result}")
    connection.execute("SELECT COUNT(*) FROM classification_usage_events").fetchone()
PY

read -r -p "Restore classification usage statistics to $ARIANE_RUNTIME_DATA_DIR? [y/N] " CONFIRM
if [ "$CONFIRM" != "y" ] && [ "$CONFIRM" != "Y" ]; then
    echo "Restore cancelled"
    exit 0
fi

TARGET="$ARIANE_RUNTIME_DATA_DIR/classification_usage.sqlite3"
PREVIOUS="$ARIANE_RUNTIME_DATA_DIR/classification_usage.sqlite3.previous"
systemctl stop ariane
rm -f -- "$PREVIOUS"
if [ -f "$TARGET" ]; then
    mv "$TARGET" "$PREVIOUS"
fi
mv "$STAGING_FILE" "$TARGET"
chown "$ARIANE_USER:$ARIANE_USER" "$TARGET"
chmod 0640 "$TARGET"
systemctl start ariane

for _ in $(seq 1 15); do
    if curl --fail --silent --max-time 5 http://127.0.0.1:8000/api/health | grep -q '"status":"ok"'; then
        echo "Classification usage restore completed"
        exit 0
    fi
    sleep 1
done

systemctl stop ariane || true
rm -f -- "$TARGET"
if [ -f "$PREVIOUS" ]; then
    mv "$PREVIOUS" "$TARGET"
fi
systemctl start ariane
echo "Restore failed and the previous usage database was restored" >&2
exit 1
