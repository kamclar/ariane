#!/bin/bash

# Install the digest-pinned SpliceAI service used by ARIANE.

set -euo pipefail

ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
ARIANE_ENV_FILE="${ARIANE_ENV_FILE:-/etc/ariane/ariane.env}"
SPLICEAI_PORT="${SPLICEAI_PORT:-8082}"
SPLICEAI_WORKERS="${SPLICEAI_WORKERS:-3}"
ARIANE_PORT="${ARIANE_PORT:-8000}"
PROFILE_FILE="$ARIANE_HOME/data/spliceai/enigma_v1_2_spliceai_profile.json"
VALIDATOR="$ARIANE_HOME/scripts/validate_spliceai_service.py"
STARTUP_VALIDATOR="$ARIANE_HOME/scripts/server-ops/wait-and-validate-spliceai-service.sh"
REFERENCE_ANNOTATION_SOURCE="$ARIANE_HOME/data/spliceai/gencode_v49_basic_reference_transcripts.tsv"
REFERENCE_ANNOTATION_DIR="/var/lib/ariane/spliceai"
REFERENCE_ANNOTATION_GZIP="$REFERENCE_ANNOTATION_DIR/gencode.v49.basic.reference-transcripts.annotation.txt.gz"
SERVICE_FILE="/etc/systemd/system/ariane-spliceai.service"
VALIDATION_CONTAINER="ariane-spliceai-validation-$$"
VALIDATION_PORT="${SPLICEAI_VALIDATION_PORT:-18081}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
for command_name in docker python3 curl systemctl ss gzip sha256sum; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "Required command is missing: $command_name" >&2
        exit 1
    fi
done
if [ ! -f "$PROFILE_FILE" ] || [ ! -f "$VALIDATOR" ] || [ ! -f "$STARTUP_VALIDATOR" ] || [ ! -f "$REFERENCE_ANNOTATION_SOURCE" ]; then
    echo "ARIANE SpliceAI profile or validator is missing" >&2
    exit 1
fi

if [ ! -f "$ARIANE_ENV_FILE" ]; then
    echo "ARIANE environment file is missing: $ARIANE_ENV_FILE" >&2
    exit 1
fi
if ! [[ "$SPLICEAI_PORT" =~ ^[0-9]+$ ]] || [ "$SPLICEAI_PORT" -lt 1024 ] || [ "$SPLICEAI_PORT" -gt 65535 ]; then
    echo "Invalid SPLICEAI_PORT: $SPLICEAI_PORT" >&2
    exit 1
fi
if ! [[ "$VALIDATION_PORT" =~ ^[0-9]+$ ]] || [ "$VALIDATION_PORT" -lt 1024 ] || [ "$VALIDATION_PORT" -gt 65535 ]; then
    echo "Invalid SPLICEAI_VALIDATION_PORT: $VALIDATION_PORT" >&2
    exit 1
fi
if ! [[ "$ARIANE_PORT" =~ ^[0-9]+$ ]] || [ "$ARIANE_PORT" -lt 1024 ] || [ "$ARIANE_PORT" -gt 65535 ]; then
    echo "Invalid ARIANE_PORT: $ARIANE_PORT" >&2
    exit 1
fi
if ! [[ "$SPLICEAI_WORKERS" =~ ^[0-9]+$ ]] || [ "$SPLICEAI_WORKERS" -lt 1 ] || [ "$SPLICEAI_WORKERS" -gt 8 ]; then
    echo "Invalid SPLICEAI_WORKERS: $SPLICEAI_WORKERS" >&2
    exit 1
fi
if [ "$SPLICEAI_PORT" = "$ARIANE_PORT" ] || [ "$SPLICEAI_PORT" = "$VALIDATION_PORT" ]; then
    echo "SpliceAI port conflicts with another ARIANE port: $SPLICEAI_PORT" >&2
    exit 1
fi

SPLICEAI_IMAGE="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["approved_engine"]["docker_image"])' "$PROFILE_FILE")"
if ! [[ "$SPLICEAI_IMAGE" =~ ^docker\.io/weisburd/spliceai-38@sha256:[0-9a-f]{64}$ ]]; then
    echo "The active profile does not contain an immutable SpliceAI GRCh38 image: $SPLICEAI_IMAGE" >&2
    exit 1
fi

EXPECTED_ANNOTATION_SHA256="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["approved_engine"]["reference_transcript_annotation_sha256"])' "$PROFILE_FILE")"
ACTUAL_ANNOTATION_SHA256="$(sha256sum "$REFERENCE_ANNOTATION_SOURCE" | awk '{print $1}')"
if [ "$ACTUAL_ANNOTATION_SHA256" != "$EXPECTED_ANNOTATION_SHA256" ]; then
    echo "Reference-transcript SpliceAI annotation checksum does not match the active profile" >&2
    exit 1
fi
install -d -m 0755 -o root -g root "$REFERENCE_ANNOTATION_DIR"
ANNOTATION_TEMP="$(mktemp /tmp/ariane-spliceai-annotation.XXXXXX)"
gzip -n -9 -c "$REFERENCE_ANNOTATION_SOURCE" > "$ANNOTATION_TEMP"
install -m 0644 -o root -g root "$ANNOTATION_TEMP" "${REFERENCE_ANNOTATION_GZIP}.new"
mv "${REFERENCE_ANNOTATION_GZIP}.new" "$REFERENCE_ANNOTATION_GZIP"
rm -f "$ANNOTATION_TEMP"

cleanup_validation_container() {
    docker rm -f "$VALIDATION_CONTAINER" >/dev/null 2>&1 || true
}
trap cleanup_validation_container EXIT

echo "Pulling $SPLICEAI_IMAGE"
docker pull "$SPLICEAI_IMAGE"
docker image inspect "$SPLICEAI_IMAGE" >/dev/null

echo "Validating the pinned image on 127.0.0.1:$VALIDATION_PORT"
docker run --detach \
    --name "$VALIDATION_CONTAINER" \
    --init \
    --cap-drop ALL \
    --security-opt no-new-privileges:true \
    --pids-limit 512 \
    --publish "127.0.0.1:${VALIDATION_PORT}:8080" \
    --mount "type=bind,src=${REFERENCE_ANNOTATION_GZIP},dst=/gencode.v49.basic.annotation.txt.gz,readonly" \
    --env GENE_SET=basic \
    --env WORKERS=1 \
    --env DATABASE_ENABLED=0 \
    --env DISABLE_RATE_LIMIT=1 \
    "$SPLICEAI_IMAGE" >/dev/null

READY=0
for _ in $(seq 1 180); do
    if curl --silent --output /dev/null --max-time 2 \
        "http://127.0.0.1:${VALIDATION_PORT}/"; then
        READY=1
        break
    fi
    if ! docker inspect --format '{{.State.Running}}' "$VALIDATION_CONTAINER" 2>/dev/null | grep -q '^true$'; then
        echo "SpliceAI validation container stopped during startup" >&2
        docker logs "$VALIDATION_CONTAINER" >&2 || true
        exit 1
    fi
    sleep 2
done
if [ "$READY" -ne 1 ]; then
    echo "SpliceAI validation container did not become ready" >&2
    docker logs "$VALIDATION_CONTAINER" >&2 || true
    exit 1
fi
python3 "$VALIDATOR" \
    --url "http://127.0.0.1:${VALIDATION_PORT}/spliceai/" \
    --timeout 180
cleanup_validation_container

if ss -H -ltn "sport = :${SPLICEAI_PORT}" | grep -q .; then
    CURRENT_SPLICEAI_MAPPING=""
    if systemctl is-active --quiet ariane-spliceai.service; then
        CURRENT_SPLICEAI_MAPPING="$(docker port ariane-spliceai 8080/tcp 2>/dev/null || true)"
    fi
    if [ "$CURRENT_SPLICEAI_MAPPING" != "127.0.0.1:${SPLICEAI_PORT}" ]; then
        echo "SpliceAI port is already used by another service: 127.0.0.1:${SPLICEAI_PORT}" >&2
        ss -ltnp "sport = :${SPLICEAI_PORT}" >&2 || true
        exit 1
    fi
fi

SERVICE_BACKUP=""
SERVICE_FILE_EXISTED=0
SERVICE_WAS_ACTIVE=0
SERVICE_WAS_ENABLED=0
if [ -f "$SERVICE_FILE" ]; then
    SERVICE_BACKUP="$(mktemp /tmp/ariane-spliceai-service-backup.XXXXXX)"
    cp "$SERVICE_FILE" "$SERVICE_BACKUP"
    SERVICE_FILE_EXISTED=1
fi
if systemctl is-active --quiet ariane-spliceai.service; then
    SERVICE_WAS_ACTIVE=1
fi
if systemctl is-enabled --quiet ariane-spliceai.service; then
    SERVICE_WAS_ENABLED=1
fi

restore_service() {
    local exit_code="$?"
    trap - ERR
    set +e
    systemctl stop ariane-spliceai.service >/dev/null 2>&1
    if [ "$SERVICE_FILE_EXISTED" -eq 1 ]; then
        install -m 0644 -o root -g root "$SERVICE_BACKUP" "$SERVICE_FILE"
    else
        systemctl disable ariane-spliceai.service >/dev/null 2>&1
        rm -f "$SERVICE_FILE"
    fi
    systemctl daemon-reload
    if [ "$SERVICE_WAS_ENABLED" -eq 1 ]; then
        systemctl enable ariane-spliceai.service >/dev/null 2>&1
    else
        systemctl disable ariane-spliceai.service >/dev/null 2>&1
    fi
    if [ "$SERVICE_WAS_ACTIVE" -eq 1 ]; then
        systemctl start ariane-spliceai.service
    fi
    systemctl reset-failed ariane-spliceai.service >/dev/null 2>&1
    if [ -n "$SERVICE_BACKUP" ]; then
        rm -f "$SERVICE_BACKUP"
    fi
    exit "$exit_code"
}
trap restore_service ERR

TEMP_SERVICE="$(mktemp /tmp/ariane-spliceai-service.XXXXXX)"
cat > "$TEMP_SERVICE" <<EOF
[Unit]
Description=ARIANE local SpliceAI service
After=docker.service network.target
Requires=docker.service

[Service]
Type=simple
ExecStartPre=-/usr/bin/docker rm -f ariane-spliceai
ExecStart=/usr/bin/docker run --name ariane-spliceai --init --cap-drop ALL --security-opt no-new-privileges:true --pids-limit 512 --publish 127.0.0.1:${SPLICEAI_PORT}:8080 --mount type=bind,src=${REFERENCE_ANNOTATION_GZIP},dst=/gencode.v49.basic.annotation.txt.gz,readonly --env GENE_SET=basic --env WORKERS=${SPLICEAI_WORKERS} --env DATABASE_ENABLED=0 --env DISABLE_RATE_LIMIT=1 ${SPLICEAI_IMAGE}
ExecStartPost=/usr/bin/env ARIANE_HOME=${ARIANE_HOME} SPLICEAI_PORT=${SPLICEAI_PORT} SPLICEAI_WORKERS=${SPLICEAI_WORKERS} SPLICEAI_REFERENCE_ANNOTATION_GZIP=${REFERENCE_ANNOTATION_GZIP} /bin/bash ${STARTUP_VALIDATOR}
ExecStop=/usr/bin/docker stop --time 60 ariane-spliceai
ExecStopPost=-/usr/bin/docker rm -f ariane-spliceai
Restart=on-failure
RestartSec=10
TimeoutStartSec=900
TimeoutStopSec=90

[Install]
WantedBy=multi-user.target
EOF
install -m 0644 -o root -g root "$TEMP_SERVICE" "$SERVICE_FILE"
rm -f "$TEMP_SERVICE"

systemctl daemon-reload
systemctl enable ariane-spliceai.service
systemctl restart ariane-spliceai.service

trap - ERR
if [ -n "$SERVICE_BACKUP" ]; then
    rm -f "$SERVICE_BACKUP"
fi

ENV_BACKUP="$(mktemp /tmp/ariane-env.XXXXXX)"
cp "$ARIANE_ENV_FILE" "$ENV_BACKUP"
restore_environment() {
    trap - ERR
    cp "$ENV_BACKUP" "$ARIANE_ENV_FILE"
    rm -f "$ENV_BACKUP"
    if systemctl cat ariane.service >/dev/null 2>&1; then
        systemctl restart ariane.service || true
    fi
}
trap 'restore_environment; cleanup_validation_container' ERR

set_value() {
    local name="$1"
    local value="$2"
    if grep -q "^${name}=" "$ARIANE_ENV_FILE"; then
        sed -i "s|^${name}=.*$|${name}=${value}|" "$ARIANE_ENV_FILE"
    else
        printf '%s=%s\n' "$name" "$value" >> "$ARIANE_ENV_FILE"
    fi
}

set_value SPLICEAI_API_URL "http://127.0.0.1:${SPLICEAI_PORT}/spliceai/"
set_value SPLICEAI_API_SOURCE "ARIANE local SpliceAI service"
set_value SPLICEAI_API_TIMEOUT 120
set_value SPLICEAI_LOOKUP_TIMEOUT 135
set_value SPLICEAI_API_ATTEMPTS 1
set_value SPLICEAI_API_RETRY_DELAY 0
set_value SPLICEAI_API_MAX_CONCURRENT 2
set_value SPLICEAI_API_RATE_SLEEP 0

if [ "${ARIANE_RESTART_AFTER_SPLICEAI:-1}" = "1" ] && systemctl cat ariane.service >/dev/null 2>&1; then
    systemctl restart ariane.service
    HEALTHY=0
    for _ in $(seq 1 60); do
        if curl --fail --silent --max-time 5 \
            "http://127.0.0.1:${ARIANE_PORT}/api/health" | grep -q '"status":"ok"'; then
            HEALTHY=1
            break
        fi
        sleep 2
    done
    if [ "$HEALTHY" -ne 1 ]; then
        echo "ARIANE did not become healthy with the local SpliceAI configuration" >&2
        false
    fi
fi

trap - ERR
trap cleanup_validation_container EXIT
rm -f "$ENV_BACKUP"
echo "ARIANE now uses the validated local SpliceAI service on 127.0.0.1:$SPLICEAI_PORT"
