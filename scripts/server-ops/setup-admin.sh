#!/bin/bash

# Configure credentials and storage for the audit administration page.

set -euo pipefail

ARIANE_USER="${ARIANE_USER:-ubuntu}"
ADMIN_USER="${1:-admin}"
ROTATE_MODE="${2:-}"
ENV_FILE="/etc/ariane/ariane.env"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if ! [[ "$ADMIN_USER" =~ ^[A-Za-z0-9_.-]{3,64}$ ]]; then
    echo "Invalid administrator username" >&2
    exit 2
fi

install -d -m 0750 -o root -g "$ARIANE_USER" /etc/ariane
touch "$ENV_FILE"

EXISTING_PASSWORD="$(grep '^ARIANE_ADMIN_PASSWORD=' "$ENV_FILE" 2>/dev/null | tail -1 | cut -d= -f2- || true)"
if [ -n "$EXISTING_PASSWORD" ] && [ "$ROTATE_MODE" != "--rotate" ]; then
    ADMIN_PASSWORD="$EXISTING_PASSWORD"
    PASSWORD_CHANGED=false
else
    ADMIN_PASSWORD="$(openssl rand -base64 30 | tr -d '\n')"
    PASSWORD_CHANGED=true
fi

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
if [ ! -s /var/log/ariane/audit.jsonl ]; then
    journalctl -u ariane -o cat --no-pager 2>/dev/null | \
        grep -F '"log_type":"ariane_audit"' > /var/log/ariane/audit.jsonl || true
    chown "$ARIANE_USER:$ARIANE_USER" /var/log/ariane/audit.jsonl
    chmod 0640 /var/log/ariane/audit.jsonl
fi

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

install -m 0750 -o root -g root "$SCRIPT_DIR/admin-status.sh" /usr/local/sbin/ariane-admin-status
cat > /etc/systemd/system/ariane-admin-status.service <<'EOF'
[Unit]
Description=ARIANE admin status snapshot

[Service]
Type=oneshot
User=root
ExecStart=/usr/local/sbin/ariane-admin-status
RuntimeDirectory=ariane
RuntimeDirectoryMode=0755
RuntimeDirectoryPreserve=yes
NoNewPrivileges=true
PrivateTmp=true
ProtectSystem=strict
ProtectHome=true
ReadWritePaths=/run/ariane
EOF

cat > /etc/systemd/system/ariane-admin-status.timer <<'EOF'
[Unit]
Description=Update ARIANE admin status snapshot

[Timer]
OnBootSec=30s
OnUnitActiveSec=5m
Unit=ariane-admin-status.service

[Install]
WantedBy=timers.target
EOF

systemctl daemon-reload
systemctl enable --now ariane-admin-status.timer
systemctl start ariane-admin-status.service
systemctl restart ariane

echo "Admin URL: https://ariane-app.duckdns.org/admin/audit"
echo "Username: $ADMIN_USER"
if [ "$PASSWORD_CHANGED" = true ]; then
    echo "Password: $ADMIN_PASSWORD"
    echo "Store the password in a password manager. It will not be shown again."
else
    echo "The existing password was kept. Use --rotate to generate a new password."
fi
