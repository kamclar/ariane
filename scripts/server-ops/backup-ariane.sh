#!/bin/bash

# ARIANE Backup Script
# Backs up data and cache files

set -e

# Configuration
ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"
BACKUP_RETENTION_DAYS="${BACKUP_RETENTION_DAYS:-30}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${BLUE}   ARIANE Backup Script${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"

# Create backup directory
if [ ! -d "$BACKUP_DIR" ]; then
    echo -e "${YELLOW}Creating backup directory: $BACKUP_DIR${NC}"
    mkdir -p "$BACKUP_DIR"
fi

# Check permissions
if [ ! -w "$BACKUP_DIR" ]; then
    echo -e "${RED}Error: No write permission to $BACKUP_DIR${NC}"
    exit 1
fi

# Backup 1: Data directory
echo -e "${YELLOW}[1] Backing up data directory...${NC}"
DATA_DIR="$ARIANE_HOME/backend/data"
if [ -d "$DATA_DIR" ]; then
    BACKUP_FILE="$BACKUP_DIR/ariane-data-$TIMESTAMP.tar.gz"
    tar -czf "$BACKUP_FILE" -C "$ARIANE_HOME/backend" data/ 2>/dev/null
    BACKUP_SIZE=$(du -h "$BACKUP_FILE" | awk '{print $1}')
    echo -e "   ${GREEN}✓ Created: $(basename $BACKUP_FILE) (${BACKUP_SIZE})${NC}"
else
    echo -e "   ${RED}✗ Data directory not found: $DATA_DIR${NC}"
fi

# Backup 2: Full application directory
echo -e "${YELLOW}[2] Backing up full application...${NC}"
if [ -d "$ARIANE_HOME" ]; then
    BACKUP_FILE="$BACKUP_DIR/ariane-full-$TIMESTAMP.tar.gz"
    tar -czf "$BACKUP_FILE" \
        --exclude='venv' \
        --exclude='__pycache__' \
        --exclude='.git' \
        --exclude='*.pyc' \
        -C "$(dirname $ARIANE_HOME)" "$(basename $ARIANE_HOME)" 2>/dev/null
    BACKUP_SIZE=$(du -h "$BACKUP_FILE" | awk '{print $1}')
    echo -e "   ${GREEN}✓ Created: $(basename $BACKUP_FILE) (${BACKUP_SIZE})${NC}"
else
    echo -e "   ${RED}✗ Application directory not found: $ARIANE_HOME${NC}"
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
    ((DELETED_COUNT++))
done < <(find "$BACKUP_DIR" -name "ariane-*.tar.gz" -mtime +${BACKUP_RETENTION_DAYS} 2>/dev/null)

if [ $DELETED_COUNT -eq 0 ]; then
    echo -e "   ${GREEN}No old backups to delete${NC}"
fi

# Summary
echo -e "\n${BLUE}════════════════════════════════════════${NC}"
TOTAL_SIZE=$(du -sh "$BACKUP_DIR" | awk '{print $1}')
BACKUP_COUNT=$(ls -1 "$BACKUP_DIR"/ariane-*.tar.gz 2>/dev/null | wc -l)
echo -e "${GREEN}Backup completed at $(date)${NC}"
echo -e "Total backups: ${BACKUP_COUNT}"
echo -e "Total size: ${TOTAL_SIZE}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"
