#!/bin/bash

# ARIANE Backup Cron Setup Script
# Configures automatic daily backups

set -e

# Configuration
ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"
BACKUP_TIME="${1:-03:00}"  # Default: 3:00 AM
ARIANE_USER="${ARIANE_USER:-ubuntu}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${BLUE}   Automatic Backup Setup${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root (use sudo)${NC}"
    exit 1
fi

# 1. Create backup directory
echo -e "${YELLOW}[1] Setting up backup directory...${NC}"
if [ ! -d "$BACKUP_DIR" ]; then
    mkdir -p "$BACKUP_DIR"
    chmod 755 "$BACKUP_DIR"
    echo -e "${GREEN}✓ Created: $BACKUP_DIR${NC}"
else
    echo -e "${GREEN}✓ Using existing: $BACKUP_DIR${NC}"
fi

# 2. Create backup script
echo -e "\n${YELLOW}[2] Creating backup script...${NC}"
BACKUP_SCRIPT="/usr/local/bin/ariane-backup"
BACKUP_LOG="$BACKUP_DIR/backup.log"
cat > "$BACKUP_SCRIPT" << 'SCRIPT'
#!/bin/bash

# ARIANE Automated Backup Script
ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$BACKUP_DIR/backup.log"

mkdir -p "$BACKUP_DIR"

echo "[$(date)] Starting backup..." >> "$LOG_FILE"

# Backup data directory
if [ -d "$ARIANE_HOME/backend/data" ]; then
    tar -czf "$BACKUP_DIR/ariane-data-$TIMESTAMP.tar.gz" \
        -C "$ARIANE_HOME/backend" data/ 2>/dev/null
    echo "[$(date)] Data backup completed" >> "$LOG_FILE"
fi

# Backup full application
tar -czf "$BACKUP_DIR/ariane-full-$TIMESTAMP.tar.gz" \
    --exclude='venv' \
    --exclude='__pycache__' \
    --exclude='.git' \
    --exclude='*.pyc' \
    -C "$(dirname $ARIANE_HOME)" "$(basename $ARIANE_HOME)" 2>/dev/null
echo "[$(date)] Full backup completed" >> "$LOG_FILE"

# Cleanup old backups (keep 30 days)
find "$BACKUP_DIR" -name "ariane-*.tar.gz" -mtime +30 -delete
echo "[$(date)] Old backups cleaned up" >> "$LOG_FILE"
echo "[$(date)] Backup finished successfully" >> "$LOG_FILE"
SCRIPT

chmod +x "$BACKUP_SCRIPT"
echo -e "${GREEN}✓ Backup script created: $BACKUP_SCRIPT${NC}"

# 2.5 Create restore helper
echo -e "\n${YELLOW}[2.5] Installing restore helper...${NC}"
RESTORE_SCRIPT="/usr/local/bin/restore-ariane"
cat > "$RESTORE_SCRIPT" << 'RESTORE'
#!/bin/bash

# ARIANE Restore Script
# Restores a selected ARIANE backup archive.

set -euo pipefail

ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"

usage() {
    cat <<EOF
Usage: $0 [OPTIONS] <backup-file>

Options:
  --list               List available backups in $BACKUP_DIR
  --data               Restore only the backend/data directory from a full backup
  --full               Restore the full application from a full backup
  -h, --help           Show this help message

Examples:
  $0 --list
  $0 /backup/ariane-full-20260710_030000.tar.gz
  $0 --data /backup/ariane-full-20260710_030000.tar.gz
EOF
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

RESTORE_MODE="full"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --list)
            echo "Available backups in $BACKUP_DIR:"
            ls -1 "$BACKUP_DIR"/ariane-*.tar.gz 2>/dev/null || echo "  (no backups found)"
            exit 0
            ;;
        --data)
            RESTORE_MODE="data"
            shift
            ;;
        --full)
            RESTORE_MODE="full"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            BACKUP_FILE="$1"
            shift
            ;;
    esac
done

if [ -z "${BACKUP_FILE:-}" ]; then
    echo "Error: backup file is required." >&2
    usage
    exit 1
fi

if [ ! -f "$BACKUP_FILE" ]; then
    echo "Error: backup file not found: $BACKUP_FILE" >&2
    exit 1
fi

echo "Restoring backup: $BACKUP_FILE"
echo "Restore mode: $RESTORE_MODE"

read -r -p "This will overwrite files in $ARIANE_HOME. Continue? [y/N] " CONFIRM
if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
    echo "Aborted."
    exit 0
fi

if [ "$RESTORE_MODE" = "data" ]; then
    echo "Extracting backend/data from backup..."
    tar -xzf "$BACKUP_FILE" -C "$(dirname "$ARIANE_HOME")" --strip-components=1 "$(basename "$ARIANE_HOME")/backend/data"
    echo "Data restore complete."
else
    echo "Extracting full backup to $(dirname "$ARIANE_HOME")..."
    tar -xzf "$BACKUP_FILE" -C "$(dirname "$ARIANE_HOME")"
    echo "Full restore complete."
fi

echo "Restore finished."
RESTORE
chmod +x "$RESTORE_SCRIPT"
echo -e "${GREEN}✓ Restore helper installed: $RESTORE_SCRIPT${NC}"

# 3. Create cron job
echo -e "\n${YELLOW}[3] Creating cron job...${NC}"

# Parse time
HOUR=$(echo "$BACKUP_TIME" | cut -d: -f1)
MINUTE=$(echo "$BACKUP_TIME" | cut -d: -f2)

# Create cron entry
CRON_JOB="$MINUTE $HOUR * * * $BACKUP_SCRIPT"

# Add to root crontab
if crontab -l 2>/dev/null | grep -q "$BACKUP_SCRIPT"; then
    echo -e "${YELLOW}⚠ Cron job already exists${NC}"
else
    (crontab -l 2>/dev/null; echo "$CRON_JOB") | crontab -
    echo -e "${GREEN}✓ Cron job added${NC}"
fi

# 4. Verify cron job
echo -e "\n${YELLOW}[4] Verifying cron job...${NC}"
if crontab -l 2>/dev/null | grep -q "$BACKUP_SCRIPT"; then
    echo -e "${GREEN}✓ Cron job configured${NC}"
    echo -e "   Schedule: Daily at $BACKUP_TIME"
else
    echo -e "${RED}✗ Failed to configure cron job${NC}"
    exit 1
fi

# 5. Create systemd timer (alternative to cron)
echo -e "\n${YELLOW}[5] Setting up systemd timer (alternative)...${NC}"

cat > /etc/systemd/system/ariane-backup.service << 'SYSD_SERVICE'
[Unit]
Description=ARIANE Backup Service
After=network.target

[Service]
Type=oneshot
User=root
ExecStart=/usr/local/bin/ariane-backup
StandardOutput=journal
StandardError=journal
SYSD_SERVICE

cat > /etc/systemd/system/ariane-backup.timer << SYSD_TIMER
[Unit]
Description=Daily ARIANE Backup Timer
Requires=ariane-backup.service

[Timer]
OnCalendar=daily
OnCalendar=*-*-* $HOUR:$MINUTE:00
Persistent=true

[Install]
WantedBy=timers.target
SYSD_TIMER

systemctl daemon-reload
systemctl enable ariane-backup.timer
systemctl start ariane-backup.timer

echo -e "${GREEN}✓ Systemd timer configured${NC}"

# 6. Run initial backup
echo -e "\n${YELLOW}[6] Running initial backup...${NC}"
$BACKUP_SCRIPT
echo -e "${GREEN}✓ Initial backup completed${NC}"

# 7. Show backup log
echo -e "\n${YELLOW}[7] Backup log:${NC}"
LOG_FILE="$BACKUP_DIR/backup.log"
tail -10 "$LOG_FILE" | sed 's/^/   /'

# Summary
echo -e "\n${BLUE}════════════════════════════════════════${NC}"
echo -e "${GREEN}✓ Backup setup completed!${NC}\n"
echo -e "Configuration:"
echo -e "  Backup time: ${BACKUP_TIME} daily"
echo -e "  Backup script: ${BLUE}$BACKUP_SCRIPT${NC}"
echo -e "  Backup directory: ${BLUE}$BACKUP_DIR${NC}"
echo -e "  Retention: 30 days"
echo -e "\nUseful commands:"
echo -e "  Run backup manually: ${BLUE}sudo $BACKUP_SCRIPT${NC}"
echo -e "  View backup log: ${BLUE}tail -f $BACKUP_DIR/backup.log${NC}"
echo -e "  List backups: ${BLUE}ls -lh $BACKUP_DIR/ariane-*.tar.gz${NC}"
echo -e "  Check cron: ${BLUE}crontab -l${NC}"
echo -e "  Check timer: ${BLUE}systemctl status ariane-backup.timer${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"
