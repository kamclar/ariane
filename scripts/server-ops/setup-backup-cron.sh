#!/bin/bash

# Install one systemd timer for ARIANE backups.

set -euo pipefail

ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"
BACKUP_RETENTION_DAYS="${BACKUP_RETENTION_DAYS:-30}"
BACKUP_TIME="${1:-03:00}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if ! [[ "$BACKUP_TIME" =~ ^([01][0-9]|2[0-3]):[0-5][0-9]$ ]]; then
    echo "Invalid time. Use HH:MM" >&2
    exit 2
fi

install -d -m 0700 -o root -g root "$BACKUP_DIR"
install -m 0750 -o root -g root "$SCRIPT_DIR/backup-ariane.sh" /usr/local/sbin/ariane-backup
install -m 0750 -o root -g root "$SCRIPT_DIR/restore-ariane.sh" /usr/local/sbin/ariane-restore
install -d -m 0750 -o root -g root /etc/ariane

cat > /etc/ariane/backup.env <<EOF
ARIANE_HOME=$ARIANE_HOME
BACKUP_DIR=$BACKUP_DIR
BACKUP_RETENTION_DAYS=$BACKUP_RETENTION_DAYS
EOF
chmod 0600 /etc/ariane/backup.env

cat > /etc/systemd/system/ariane-backup.service <<'EOF'
[Unit]
Description=ARIANE backup

[Service]
Type=oneshot
User=root
EnvironmentFile=/etc/ariane/backup.env
ExecStart=/usr/local/sbin/ariane-backup
UMask=0077
Nice=10
IOSchedulingClass=idle
NoNewPrivileges=true
PrivateTmp=true
ProtectSystem=strict
ProtectHome=read-only
ReadWritePaths=/backup
EOF

cat > /etc/systemd/system/ariane-backup.timer <<EOF
[Unit]
Description=Daily ARIANE backup

[Timer]
OnCalendar=*-*-* $BACKUP_TIME:00
Persistent=true
RandomizedDelaySec=5m
Unit=ariane-backup.service

[Install]
WantedBy=timers.target
EOF

# Remove the obsolete cron entry to prevent duplicate backups.
if command -v crontab >/dev/null && crontab -l >/tmp/ariane-crontab 2>/dev/null; then
    grep -v 'ariane-backup' /tmp/ariane-crontab >/tmp/ariane-crontab.new || true
    crontab /tmp/ariane-crontab.new
    rm -f /tmp/ariane-crontab /tmp/ariane-crontab.new
fi

systemctl daemon-reload
systemctl enable --now ariane-backup.timer
systemctl start ariane-backup.service
systemctl list-timers ariane-backup.timer --no-pager
