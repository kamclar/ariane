#!/bin/bash

# Restart ARIANE without stopping Nginx.

set -euo pipefail

ARIANE_PORT="${ARIANE_PORT:-8000}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi

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
