#!/bin/bash

# Show structured ARIANE audit events from the systemd journal.

set -euo pipefail

MODE="${1:-all}"
SINCE="${2:-today}"

case "$MODE" in
    all) FILTER='.' ;;
    errors) FILTER='select(.event | test("error|exception"))' ;;
    requests) FILTER='select(.event == "classification_completed" or .event == "batch_item_completed" or .event == "manual_evidence_completed")' ;;
    *)
        echo "Usage: $0 [all|errors|requests] [journalctl-since]" >&2
        exit 2
        ;;
esac

if ! command -v jq >/dev/null; then
    echo "jq is required" >&2
    exit 1
fi

printf 'REQUEST_ID\tSOURCE_IP\tEVENT\tINPUT\tRESULT\tERROR\n'
journalctl -u ariane --since "$SINCE" -o cat --no-pager | \
    jq -Rr "fromjson? | select(.log_type == \"ariane_audit\") | $FILTER | [(.request_id // \"\"), (.source_ip // \"\"), (.event // \"\"), ((.input // {}) | tojson), ((.result // {}) | tojson), (.error // \"\")] | @tsv"
