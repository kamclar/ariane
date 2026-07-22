#!/bin/bash

# ARIANE Service Installation Script
# Installs ARIANE as a systemd service

set -e

# Configuration
ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
ARIANE_USER="${ARIANE_USER:-ubuntu}"
ARIANE_PORT="${ARIANE_PORT:-8000}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   ARIANE Service Installation${NC}"
echo -e "${BLUE}========================================${NC}\n"

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root (use sudo)${NC}"
    exit 1
fi

# 1. Check if ARIANE directory exists
echo -e "${YELLOW}[1] Checking ARIANE installation...${NC}"
if [ ! -d "$ARIANE_HOME" ]; then
    echo -e "${RED}Error: ARIANE directory not found: $ARIANE_HOME${NC}"
    exit 1
fi
echo -e "${GREEN}OK Found: $ARIANE_HOME${NC}"

# 2. Check if venv exists
if [ ! -d "$ARIANE_HOME/venv" ]; then
    echo -e "${YELLOW}WARN Virtual environment not found, creating...${NC}"
    cd "$ARIANE_HOME"
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    deactivate
fi
echo -e "${GREEN}OK Virtual environment ready${NC}"

# 3. Get Python path
PYTHON_PATH="$ARIANE_HOME/venv/bin/python"
if [ ! -f "$PYTHON_PATH" ]; then
    echo -e "${RED}Error: Python executable not found at $PYTHON_PATH${NC}"
    exit 1
fi
echo -e "${GREEN}OK Python: $PYTHON_PATH${NC}"

# 4. Create the environment file and systemd service.
echo -e "\n${YELLOW}[2] Creating systemd service...${NC}"
install -d -m 0750 -o root -g "$ARIANE_USER" /etc/ariane
install -d -m 0750 -o "$ARIANE_USER" -g "$ARIANE_USER" /var/log/ariane
install -d -m 0750 -o "$ARIANE_USER" -g "$ARIANE_USER" /var/lib/ariane/runtime-cache
touch /var/log/ariane/audit.jsonl
chown "$ARIANE_USER:$ARIANE_USER" /var/log/ariane/audit.jsonl
chmod 0640 /var/log/ariane/audit.jsonl
if [ ! -f /etc/ariane/ariane.env ]; then
    printf 'ARIANE_ADMIN_TOKEN=%s\n' "$(openssl rand -hex 32)" > /etc/ariane/ariane.env
fi
if ! grep -q '^ARIANE_RUNTIME_CACHE_DIR=' /etc/ariane/ariane.env; then
    printf 'ARIANE_RUNTIME_CACHE_DIR=/var/lib/ariane/runtime-cache\n' >> /etc/ariane/ariane.env
fi
chown root:"$ARIANE_USER" /etc/ariane/ariane.env
chmod 0640 /etc/ariane/ariane.env
cat > /etc/systemd/system/ariane.service << EOF
[Unit]
Description=ARIANE FastAPI Application
After=network.target

[Service]
Type=simple
User=$ARIANE_USER
WorkingDirectory=$ARIANE_HOME
Environment="PATH=$ARIANE_HOME/venv/bin"
EnvironmentFile=-/etc/ariane/ariane.env
ExecStart=$ARIANE_HOME/venv/bin/uvicorn backend.main:app --host 127.0.0.1 --port $ARIANE_PORT --workers 4 --proxy-headers --forwarded-allow-ips=127.0.0.1
Restart=always
RestartSec=10
TimeoutStopSec=30
UMask=0027
NoNewPrivileges=true
PrivateTmp=true
PrivateDevices=true
ProtectSystem=strict
ProtectHome=read-only
ReadWritePaths=/var/lib/ariane/runtime-cache
ProtectKernelTunables=true
ProtectKernelModules=true
ProtectControlGroups=true
RestrictSUIDSGID=true
LockPersonality=true
RestrictRealtime=true
CapabilityBoundingSet=
AmbientCapabilities=
ReadWritePaths=$ARIANE_HOME/backend/data /var/log/ariane
StandardOutput=journal
StandardError=journal

[Install]
WantedBy=multi-user.target
EOF

echo -e "${GREEN}OK Service file created: /etc/systemd/system/ariane.service${NC}"

# 5. Reload systemd
echo -e "\n${YELLOW}[3] Reloading systemd daemon...${NC}"
systemctl daemon-reload
echo -e "${GREEN}OK Systemd reloaded${NC}"

# 6. Enable service
echo -e "\n${YELLOW}[4] Enabling service on boot...${NC}"
systemctl enable ariane
echo -e "${GREEN}OK Service enabled${NC}"

# 7. Start service
echo -e "\n${YELLOW}[5] Starting ARIANE service...${NC}"
systemctl start ariane
sleep 2

# 8. Check status
echo -e "\n${YELLOW}[6] Checking service status...${NC}"
if systemctl is-active --quiet ariane; then
    echo -e "${GREEN}OK ARIANE service is running${NC}"
else
    echo -e "${RED}ERROR Service failed to start${NC}"
    echo -e "${YELLOW}Check logs with: journalctl -u ariane -f${NC}"
    exit 1
fi

# 9. Health check
echo -e "\n${YELLOW}[7] Running health check...${NC}"
sleep 2
if command -v curl &> /dev/null; then
    for i in {1..5}; do
        if curl --fail --silent --max-time 5 "http://127.0.0.1:$ARIANE_PORT/api/health" | grep -q '"status":"ok"'; then
            echo -e "${GREEN}OK API is responding${NC}"
            break
        fi
        if [ $i -lt 5 ]; then
            echo -e "   Attempt $i/5..."
            sleep 1
        fi
    done
    if ! curl --fail --silent --max-time 5 "http://127.0.0.1:$ARIANE_PORT/api/health" | grep -q '"status":"ok"'; then
        echo "Health check failed" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}WARN curl not available, skipping health check${NC}"
fi

# Summary
echo -e "\n${BLUE}========================================${NC}"
echo -e "${GREEN}Installation completed successfully!${NC}\n"
echo -e "Service commands:"
echo -e "  Status:  ${BLUE}systemctl status ariane${NC}"
echo -e "  Start:   ${BLUE}systemctl start ariane${NC}"
echo -e "  Stop:    ${BLUE}systemctl stop ariane${NC}"
echo -e "  Restart: ${BLUE}systemctl restart ariane${NC}"
echo -e "  Logs:    ${BLUE}journalctl -u ariane -f${NC}"
echo -e "\nApplication URL: ${BLUE}http://localhost:$ARIANE_PORT${NC}"
echo -e "${BLUE}========================================${NC}\n"
