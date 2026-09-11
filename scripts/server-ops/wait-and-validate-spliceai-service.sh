#!/bin/bash

# Keep the SpliceAI systemd unit in "starting" state until every model worker
# has completed a profile validation request. This also initializes TensorFlow
# before ARIANE accepts user traffic.

set -euo pipefail

ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
SPLICEAI_PORT="${SPLICEAI_PORT:-8082}"
SPLICEAI_WORKERS="${SPLICEAI_WORKERS:-3}"
VALIDATOR="$ARIANE_HOME/scripts/validate_spliceai_service.py"
PROFILE_FILE="$ARIANE_HOME/data/spliceai/enigma_v1_2_spliceai_profile.json"
REFERENCE_ANNOTATION_GZIP="${SPLICEAI_REFERENCE_ANNOTATION_GZIP:-/var/lib/ariane/spliceai/gencode.v49.basic.reference-transcripts.annotation.txt.gz}"

if ! [[ "$SPLICEAI_PORT" =~ ^[0-9]+$ ]] || [ "$SPLICEAI_PORT" -lt 1024 ] || [ "$SPLICEAI_PORT" -gt 65535 ]; then
    echo "Invalid SPLICEAI_PORT: $SPLICEAI_PORT" >&2
    exit 1
fi
if ! [[ "$SPLICEAI_WORKERS" =~ ^[0-9]+$ ]] || [ "$SPLICEAI_WORKERS" -lt 1 ] || [ "$SPLICEAI_WORKERS" -gt 8 ]; then
    echo "Invalid SPLICEAI_WORKERS: $SPLICEAI_WORKERS" >&2
    exit 1
fi
if [ ! -f "$VALIDATOR" ] || [ ! -f "$PROFILE_FILE" ] || [ ! -f "$REFERENCE_ANNOTATION_GZIP" ]; then
    echo "SpliceAI validator, profile, or reference-transcript annotation is missing" >&2
    exit 1
fi

EXPECTED_ANNOTATION_SHA256="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["approved_engine"]["reference_transcript_annotation_sha256"])' "$PROFILE_FILE")"
ACTUAL_ANNOTATION_SHA256="$(gzip -cd "$REFERENCE_ANNOTATION_GZIP" | sha256sum | awk '{print $1}')"
if [ "$ACTUAL_ANNOTATION_SHA256" != "$EXPECTED_ANNOTATION_SHA256" ]; then
    echo "Installed reference-transcript SpliceAI annotation checksum does not match the active profile" >&2
    exit 1
fi

READY=0
for _ in $(seq 1 180); do
    if curl --silent --output /dev/null --max-time 2 \
        "http://127.0.0.1:${SPLICEAI_PORT}/"; then
        READY=1
        break
    fi
    sleep 2
done
if [ "$READY" -ne 1 ]; then
    echo "ARIANE SpliceAI service did not become ready" >&2
    exit 1
fi

python3 "$VALIDATOR" \
    --url "http://127.0.0.1:${SPLICEAI_PORT}/spliceai/" \
    --timeout 180 \
    --concurrent-copies "$SPLICEAI_WORKERS"
