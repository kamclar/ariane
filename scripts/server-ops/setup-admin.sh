#!/bin/bash

# Configure credentials and storage for the audit administration page.

set -euo pipefail

ARIANE_USER="${ARIANE_USER:-ubuntu}"
ADMIN_USER="${1:-admin}"
ENV_FILE="/etc/ariane/ariane.env"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if ! [[ "$ADMIN_USER" =~ ^[A-Za-z0-9_.-]{3,64}$ ]]; then
    echo "Invalid administrator username" >&2
    exit 2
fi

ADMIN_PASSWORD="$(openssl rand -base64 30 | tr -d '\n')"
install -d -m 0750 -o root -g "$ARIANE_USER" /etc/ariane
touch "$ENV_FILE"

TEMP_FILE="$(mktemp)"
trap 'rm -f "$TEMP_FILE"' EXIT
grep -vE '^(ARIANE_ADMIN_USER|ARIANE_ADMIN_PASSWORD|ARIANE_AUDIT_LOG)=' "$ENV_FILE" > "$TEMP_FILE" || true
printf 'ARIANE_ADMIN_USER=%s\n' "$ADMIN_USER" >> "$TEMP_FILE"
printf 'ARIANE_ADMIN_PASSWORD=%s\n' "$ADMIN_PASSWORD" >> "$TEMP_FILE"
printf 'ARIANE_AUDIT_LOG=/var/log/ariane/audit.jsonl\n' >> "$TEMP_FILE"
install -m 0640 -o root -g "$ARIANE_USER" "$TEMP_FILE" "$ENV_FILE"

install -d -m 0750 -o "$ARIANE_USER" -g "$ARIANE_USER" /var/log/ariane
touch /var/log/ariane/audit.jsonl
chown "$ARIANE_USER:$ARIANE_USER" /var/log/ariane/audit.jsonl
chmod 0640 /var/log/ariane/audit.jsonl

cat > /etc/logrotate.d/ariane-audit <<EOF
/var/log/ariane/audit.jsonl {
    daily
    rotate 30
    compress
    delaycompress
    missingok
    notifempty
    copytruncate
    create 0640 $ARIANE_USER $ARIANE_USER
}
EOF

install -d -m 0755 /etc/systemd/system/ariane.service.d
cat > /etc/systemd/system/ariane.service.d/audit.conf <<'EOF'
[Service]
ReadWritePaths=/var/log/ariane
EOF

systemctl daemon-reload
systemctl restart ariane

echo "Admin URL: https://ariane-app.duckdns.org/admin/audit"
echo "Username: $ADMIN_USER"
echo "Password: $ADMIN_PASSWORD"
echo "Store the password in a password manager. It will not be shown again."
