#!/bin/bash

# ARIANE Backup Script
# Backs up data and cache files

set -euo pipefail

# Configuration
ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"
BACKUP_RETENTION_DAYS="${BACKUP_RETENTION_DAYS:-30}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOCK_FILE="${BACKUP_DIR}/.ariane-backup.lock"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   ARIANE Backup Script${NC}"
echo -e "${BLUE}========================================${NC}\n"

# Create backup directory
if [ ! -d "$BACKUP_DIR" ]; then
    echo -e "${YELLOW}Creating backup directory: $BACKUP_DIR${NC}"
    install -d -m 0700 "$BACKUP_DIR"
fi

# Check permissions
if [ ! -w "$BACKUP_DIR" ]; then
    echo -e "${RED}Error: No write permission to $BACKUP_DIR${NC}"
    exit 1
fi

exec 9>"$LOCK_FILE"
if ! flock -n 9; then
    echo "Error: another backup is already running" >&2
    exit 1
fi

# Backup 1: Data directory
echo -e "${YELLOW}[1] Backing up data directory...${NC}"
DATA_DIR="$ARIANE_HOME/backend/data"
if [ -d "$DATA_DIR" ]; then
    BACKUP_FILE="$BACKUP_DIR/ariane-data-$TIMESTAMP.tar.gz"
    TEMP_FILE="${BACKUP_FILE}.tmp"
    tar -czf "$TEMP_FILE" -C "$ARIANE_HOME/backend" data/
    tar -tzf "$TEMP_FILE" >/dev/null
    mv "$TEMP_FILE" "$BACKUP_FILE"
    sha256sum "$BACKUP_FILE" > "${BACKUP_FILE}.sha256"
    BACKUP_SIZE=$(du -h "$BACKUP_FILE" | awk '{print $1}')
    echo -e "   ${GREEN}OK Created: $(basename $BACKUP_FILE) (${BACKUP_SIZE})${NC}"
else
    echo -e "   ${RED}ERROR Data directory not found: $DATA_DIR${NC}"
fi

# Backup 2: Full application directory
echo -e "${YELLOW}[2] Backing up full application...${NC}"
if [ -d "$ARIANE_HOME" ]; then
    BACKUP_FILE="$BACKUP_DIR/ariane-full-$TIMESTAMP.tar.gz"
    TEMP_FILE="${BACKUP_FILE}.tmp"
    tar -czf "$TEMP_FILE" \
        --exclude='venv' \
        --exclude='__pycache__' \
        --exclude='.git' \
        --exclude='*.pyc' \
        -C "$(dirname "$ARIANE_HOME")" "$(basename "$ARIANE_HOME")"
    tar -tzf "$TEMP_FILE" >/dev/null
    mv "$TEMP_FILE" "$BACKUP_FILE"
    sha256sum "$BACKUP_FILE" > "${BACKUP_FILE}.sha256"
    BACKUP_SIZE=$(du -h "$BACKUP_FILE" | awk '{print $1}')
    echo -e "   ${GREEN}OK Created: $(basename $BACKUP_FILE) (${BACKUP_SIZE})${NC}"
else
    echo -e "   ${RED}ERROR Application directory not found: $ARIANE_HOME${NC}"
fi

# List recent backups
echo -e "\n${YELLOW}[3] Recent Backups:${NC}"
ls -lh "$BACKUP_DIR"/ariane-*.tar.gz 2>/dev/null | tail -10 | awk '{print "   " $9 " (" $5 ")"}'

# Cleanup old backups
echo -e "\n${YELLOW}[4] Cleaning up old backups (retention: ${BACKUP_RETENTION_DAYS} days)...${NC}"
DELETED_COUNT=0
while IFS= read -r file; do
    rm -f "$file"
    echo -e "   ${YELLOW}Deleted: $(basename $file)${NC}"
    rm -f "${file}.sha256"
    DELETED_COUNT=$((DELETED_COUNT + 1))
done < <(find "$BACKUP_DIR" -type f -name "ariane-*.tar.gz" -mtime "+${BACKUP_RETENTION_DAYS}" 2>/dev/null)

if [ $DELETED_COUNT -eq 0 ]; then
    echo -e "   ${GREEN}No old backups to delete${NC}"
fi

# Summary
echo -e "\n${BLUE}========================================${NC}"
TOTAL_SIZE=$(du -sh "$BACKUP_DIR" | awk '{print $1}')
BACKUP_COUNT=$(ls -1 "$BACKUP_DIR"/ariane-*.tar.gz 2>/dev/null | wc -l)
echo -e "${GREEN}Backup completed at $(date)${NC}"
echo -e "Total backups: ${BACKUP_COUNT}"
echo -e "Total size: ${TOTAL_SIZE}"
echo -e "${BLUE}========================================${NC}\n"
