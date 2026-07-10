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
    tar -xzf "$BACKUP_FILE" -C "$ARIANE_HOME" --wildcards --no-anchored 'backend/data/*'
    echo "Data restore complete."
else
    echo "Extracting full backup to $(dirname "$ARIANE_HOME")..."
    tar -xzf "$BACKUP_FILE" -C "$(dirname "$ARIANE_HOME")"
    echo "Full restore complete."
fi

echo "Restore finished."
