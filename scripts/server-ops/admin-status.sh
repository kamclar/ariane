#!/bin/bash

# Write a read-only system status snapshot for the admin dashboard.

set -euo pipefail

SERVER_NAME="${ARIANE_SERVER_NAME:-ariane-app.duckdns.org}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"
STATUS_DIR="/run/ariane"
STATUS_FILE="$STATUS_DIR/admin-status.json"

install -d -m 0755 -o root -g root "$STATUS_DIR"

ARIANE_STATUS="$(systemctl is-active ariane 2>/dev/null || true)"
NGINX_STATUS="$(systemctl is-active nginx 2>/dev/null || true)"
DISK_STATUS="$(df -P / | awk 'NR==2 {print $5 " used, " $4 " KB free"}')"

CERT_FILE="/etc/letsencrypt/live/$SERVER_NAME/cert.pem"
if [ -f "$CERT_FILE" ]; then
    CERT_END="$(openssl x509 -in "$CERT_FILE" -noout -enddate 2>/dev/null | cut -d= -f2-)"
    if [ -n "$CERT_END" ]; then
        CERT_DAYS="$(( ($(date -d "$CERT_END" +%s) - $(date +%s)) / 86400 ))"
        CERT_STATUS="valid, $CERT_DAYS days"
    else
        CERT_STATUS="unreadable"
    fi
else
    CERT_STATUS="missing"
fi

LAST_BACKUP="$(find "$BACKUP_DIR" -maxdepth 1 -type f -name 'ariane-full-*.tar.gz' -printf '%T@ %f\n' 2>/dev/null | sort -nr | head -1 | cut -d' ' -f2- || true)"
LAST_BACKUP="${LAST_BACKUP:-none}"
GENERATED_AT="$(date --iso-8601=seconds)"
TEMP_FILE="$(mktemp "$STATUS_DIR/admin-status.XXXXXX")"

python3 - "$TEMP_FILE" "$ARIANE_STATUS" "$NGINX_STATUS" "$CERT_STATUS" "$DISK_STATUS" "$LAST_BACKUP" "$GENERATED_AT" <<'PY'
import json
import sys

path, ariane, nginx, certificate, disk, last_backup, generated_at = sys.argv[1:]
with open(path, "w", encoding="utf-8") as stream:
    json.dump(
        {
            "ariane": ariane,
            "nginx": nginx,
            "certificate": certificate,
            "disk": disk,
            "last_backup": last_backup,
            "generated_at": generated_at,
        },
        stream,
        separators=(",", ":"),
    )
PY

chmod 0644 "$TEMP_FILE"
mv "$TEMP_FILE" "$STATUS_FILE"
