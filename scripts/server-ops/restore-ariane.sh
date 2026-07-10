#!/bin/bash

# Restore an ARIANE backup after validating its checksum and paths.

set -euo pipefail

ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"
RESTORE_MODE="full"

usage() {
    echo "Usage: $0 [--list] [--data|--full] <backup-file>"
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --list)
            find "$BACKUP_DIR" -maxdepth 1 -type f -name 'ariane-*.tar.gz' -print | sort
            exit 0
            ;;
        --data) RESTORE_MODE="data" ;;
        --full) RESTORE_MODE="full" ;;
        -h|--help) usage; exit 0 ;;
        -*) echo "Unknown option: $1" >&2; exit 2 ;;
        *)
            if [ -n "${BACKUP_FILE:-}" ]; then
                echo "Only one backup file is allowed" >&2
                exit 2
            fi
            BACKUP_FILE="$1"
            ;;
    esac
    shift
done

if [ -z "${BACKUP_FILE:-}" ] || [ ! -f "$BACKUP_FILE" ]; then
    echo "Backup file not found" >&2
    exit 1
fi

CHECKSUM_FILE="${BACKUP_FILE}.sha256"
if [ ! -f "$CHECKSUM_FILE" ]; then
    echo "Checksum file not found: $CHECKSUM_FILE" >&2
    exit 1
fi
(cd "$(dirname "$BACKUP_FILE")" && sha256sum --check "$(basename "$CHECKSUM_FILE")")
tar -tzf "$BACKUP_FILE" >/dev/null
if tar -tvzf "$BACKUP_FILE" | grep -qE '^[lh]'; then
    echo "Archive contains links" >&2
    exit 1
fi

ARCHIVE_ROOT="$(basename "$ARIANE_HOME")"
while IFS= read -r entry; do
    case "$entry" in
        /*|../*|*/../*|*/..) echo "Unsafe archive path: $entry" >&2; exit 1 ;;
        "$ARCHIVE_ROOT"|"$ARCHIVE_ROOT"/*|data|data/*) ;;
        *) echo "Unexpected archive path: $entry" >&2; exit 1 ;;
    esac
done < <(tar -tzf "$BACKUP_FILE")

read -r -p "Restore $RESTORE_MODE backup to $ARIANE_HOME? [y/N] " CONFIRM
if [ "$CONFIRM" != "y" ] && [ "$CONFIRM" != "Y" ]; then
    echo "Restore cancelled"
    exit 0
fi

PARENT_DIR="$(dirname "$ARIANE_HOME")"
STAGING_DIR="$(mktemp -d "$PARENT_DIR/.ariane-restore.XXXXXX")"
trap 'rm -rf -- "$STAGING_DIR"' EXIT

tar -xzf "$BACKUP_FILE" -C "$STAGING_DIR" --no-same-owner --no-same-permissions

systemctl stop ariane
if [ "$RESTORE_MODE" = "data" ]; then
    if [ -d "$STAGING_DIR/$ARCHIVE_ROOT/backend/data" ]; then
        SOURCE_DATA="$STAGING_DIR/$ARCHIVE_ROOT/backend/data"
    elif [ -d "$STAGING_DIR/data" ]; then
        SOURCE_DATA="$STAGING_DIR/data"
    else
        echo "Data directory not found in archive" >&2
        systemctl start ariane
        exit 1
    fi
    PREVIOUS_DATA="$ARIANE_HOME/backend/data.previous"
    rm -rf -- "$PREVIOUS_DATA"
    mv "$ARIANE_HOME/backend/data" "$PREVIOUS_DATA"
    mv "$SOURCE_DATA" "$ARIANE_HOME/backend/data"
else
    if [ ! -d "$STAGING_DIR/$ARCHIVE_ROOT" ]; then
        echo "Application root not found in archive" >&2
        systemctl start ariane
        exit 1
    fi
    PREVIOUS_HOME="${ARIANE_HOME}.previous"
    rm -rf -- "$PREVIOUS_HOME"
    mv "$ARIANE_HOME" "$PREVIOUS_HOME"
    mv "$STAGING_DIR/$ARCHIVE_ROOT" "$ARIANE_HOME"
    if [ -d "$PREVIOUS_HOME/venv" ] && [ ! -d "$ARIANE_HOME/venv" ]; then
        mv "$PREVIOUS_HOME/venv" "$ARIANE_HOME/venv"
    fi
fi
systemctl start ariane

for _ in $(seq 1 15); do
    if curl --fail --silent --max-time 5 http://127.0.0.1:8000/api/health | grep -q '"status":"ok"'; then
        echo "Restore completed"
        exit 0
    fi
    sleep 1
done

# Roll back when the restored service does not become healthy.
systemctl stop ariane || true
if [ "$RESTORE_MODE" = "data" ] && [ -d "${PREVIOUS_DATA:-}" ]; then
    rm -rf -- "$ARIANE_HOME/backend/data"
    mv "$PREVIOUS_DATA" "$ARIANE_HOME/backend/data"
elif [ "$RESTORE_MODE" = "full" ] && [ -d "${PREVIOUS_HOME:-}" ]; then
    if [ -d "$ARIANE_HOME/venv" ] && [ ! -d "$PREVIOUS_HOME/venv" ]; then
        mv "$ARIANE_HOME/venv" "$PREVIOUS_HOME/venv"
    fi
    rm -rf -- "$ARIANE_HOME"
    mv "$PREVIOUS_HOME" "$ARIANE_HOME"
fi
systemctl start ariane
echo "Restore failed and the previous version was restored" >&2
exit 1
