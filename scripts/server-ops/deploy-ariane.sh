#!/bin/bash

# ARIANE Full Deployment Script
# Complete setup from scratch

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   ARIANE Full Deployment${NC}"
echo -e "${BLUE}========================================${NC}\n"

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root (use sudo)${NC}"
    exit 1
fi

# Get repository URL
REPO_URL="${1:-}"
if [ -z "$REPO_URL" ]; then
    echo -e "${RED}Usage: sudo bash deploy-ariane.sh <git-repository-url>${NC}"
    echo -e "Example: sudo bash deploy-ariane.sh https://github.com/user/ariane.git"
    exit 1
fi

# Configuration
ARIANE_HOME="/home/ubuntu/ariane"
ARIANE_USER="ubuntu"
ARIANE_SERVER_NAME="${ARIANE_SERVER_NAME:-ariane-app.duckdns.org}"
DEPLOY_LOG="/var/log/ariane-deploy.log"

if ! [[ "$ARIANE_SERVER_NAME" =~ ^[A-Za-z0-9.-]+$ ]] || [[ "$ARIANE_SERVER_NAME" == *..* ]]; then
    echo "Invalid ARIANE_SERVER_NAME" >&2
    exit 2
fi

# Initialize log
echo "ARIANE Deployment started at $(date)" > "$DEPLOY_LOG"

# 1. System Updates
echo -e "${YELLOW}[1] Updating system packages...${NC}"
apt-get update -qq >> "$DEPLOY_LOG" 2>&1
echo -e "${GREEN}OK System updated${NC}"

# 2. Install dependencies
echo -e "\n${YELLOW}[2] Installing Python and tools...${NC}"
apt-get install -y -qq \
    python3-pip \
    python3-venv \
    python3-dev \
    git \
    curl \
    wget \
    nginx \
    ufw \
    htop \
    build-essential \
    >> "$DEPLOY_LOG" 2>&1
echo -e "${GREEN}OK Dependencies installed${NC}"

# 3. Clone repository
echo -e "\n${YELLOW}[3] Cloning ARIANE repository...${NC}"
if [ -d "$ARIANE_HOME" ]; then
    echo -e "   ${YELLOW}Directory exists, pulling latest...${NC}"
    cd "$ARIANE_HOME"
    git pull >> "$DEPLOY_LOG" 2>&1
else
    git clone "$REPO_URL" "$ARIANE_HOME" >> "$DEPLOY_LOG" 2>&1
fi
chown -R "$ARIANE_USER:$ARIANE_USER" "$ARIANE_HOME"
echo -e "${GREEN}OK Repository ready at $ARIANE_HOME${NC}"

# 4. Create virtual environment
echo -e "\n${YELLOW}[4] Creating Python virtual environment...${NC}"
cd "$ARIANE_HOME"
python3 -m venv venv >> "$DEPLOY_LOG" 2>&1
chown -R "$ARIANE_USER:$ARIANE_USER" "$ARIANE_HOME/venv"
echo -e "${GREEN}OK Virtual environment created${NC}"

# 5. Install Python requirements
echo -e "\n${YELLOW}[5] Installing Python packages...${NC}"
source venv/bin/activate >> "$DEPLOY_LOG" 2>&1
pip install --upgrade pip setuptools wheel >> "$DEPLOY_LOG" 2>&1
pip install -r requirements.txt >> "$DEPLOY_LOG" 2>&1
deactivate
echo -e "${GREEN}OK Python packages installed${NC}"

# 6. Create the environment file and systemd service.
echo -e "\n${YELLOW}[6] Setting up systemd service...${NC}"
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
ExecStart=$ARIANE_HOME/venv/bin/uvicorn backend.main:app --host 127.0.0.1 --port 8000 --workers 4 --proxy-headers --forwarded-allow-ips=127.0.0.1
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

systemctl daemon-reload
systemctl enable ariane >> "$DEPLOY_LOG" 2>&1
echo -e "${GREEN}OK Service configured${NC}"

# 7. Configure Nginx
echo -e "\n${YELLOW}[7] Configuring Nginx...${NC}"
cat > /etc/nginx/sites-available/ariane << 'EOF'
limit_req_zone $binary_remote_addr zone=ariane_api:10m rate=5r/s;

# Reject direct IP access and unknown hostnames.
server {
    listen 80 default_server;
    listen [::]:80 default_server;
    server_name _;
    return 444;
}

server {
    listen 80;
    listen [::]:80;

    server_name ARIANE_SERVER_NAME_PLACEHOLDER;
    server_tokens off;

    # Limit request size
    client_max_body_size 50M;

    # Timeouts
    proxy_connect_timeout 60s;
    proxy_send_timeout 60s;
    proxy_read_timeout 60s;

    location / {
        limit_req zone=ariane_api burst=20 nodelay;
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_hide_header X-Powered-By;
        add_header X-Content-Type-Options nosniff always;
        add_header X-Frame-Options DENY always;
        add_header Referrer-Policy no-referrer always;
        
        # Gzip compression
        gzip on;
        gzip_types text/plain text/css application/json application/javascript;
    }

    # Health check endpoint
    location /health {
        access_log off;
        proxy_pass http://127.0.0.1:8000/api/health;
    }

    location = /api/clear-cache {
        deny all;
    }

    # Deny access to sensitive files
    location ~ /\.env {
        deny all;
    }
    
    location ~ /\.git {
        deny all;
    }
}
EOF

sed -i "s/ARIANE_SERVER_NAME_PLACEHOLDER/$ARIANE_SERVER_NAME/" /etc/nginx/sites-available/ariane

rm -f /etc/nginx/sites-enabled/default
ln -sf /etc/nginx/sites-available/ariane /etc/nginx/sites-enabled/ariane
nginx -t >> "$DEPLOY_LOG" 2>&1
systemctl enable nginx >> "$DEPLOY_LOG" 2>&1
echo -e "${GREEN}OK Nginx configured${NC}"

# Allow SSH and web traffic. Keep Uvicorn private.
ufw allow OpenSSH >> "$DEPLOY_LOG" 2>&1
ufw allow 'Nginx Full' >> "$DEPLOY_LOG" 2>&1
ufw deny 8000/tcp >> "$DEPLOY_LOG" 2>&1
ufw --force enable >> "$DEPLOY_LOG" 2>&1

# 8. Create backup directory
echo -e "\n${YELLOW}[8] Setting up backup directory...${NC}"
install -d -m 0700 -o root -g root /backup
echo -e "${GREEN}OK Backup directory ready${NC}"

# 9. Start services
echo -e "\n${YELLOW}[9] Starting services...${NC}"
systemctl start ariane >> "$DEPLOY_LOG" 2>&1
sleep 3
systemctl start nginx >> "$DEPLOY_LOG" 2>&1
echo -e "${GREEN}OK Services started${NC}"

# 10. Health check
echo -e "\n${YELLOW}[10] Running health check...${NC}"
sleep 3
for i in {1..10}; do
    if curl --fail --silent --max-time 5 http://127.0.0.1:8000/api/health | grep -q '"status":"ok"'; then
        echo -e "${GREEN}OK API is responding${NC}"
        break
    fi
    if [ $i -lt 10 ]; then
        echo -e "   Waiting... attempt $i/10"
        sleep 1
    fi
done
if ! curl --fail --silent --max-time 5 http://127.0.0.1:8000/api/health | grep -q '"status":"ok"'; then
    echo "Health check failed" >&2
    exit 1
fi

# 11. Check status
echo -e "\n${YELLOW}[11] Verifying deployment...${NC}"
if systemctl is-active --quiet ariane; then
    echo -e "   ${GREEN}OK ARIANE service: RUNNING${NC}"
else
    echo -e "   ${RED}ERROR ARIANE service: FAILED${NC}"
fi

if systemctl is-active --quiet nginx; then
    echo -e "   ${GREEN}OK Nginx: RUNNING${NC}"
else
    echo -e "   ${RED}ERROR Nginx: FAILED${NC}"
fi

# Summary
echo -e "\n${BLUE}========================================${NC}"
echo -e "${GREEN}OK Deployment completed successfully!${NC}\n"
echo -e "Application URL: ${BLUE}http://$ARIANE_SERVER_NAME${NC}"
echo -e "API URL: ${BLUE}http://$ARIANE_SERVER_NAME/api${NC}"
echo -e "Installation directory: ${BLUE}$ARIANE_HOME${NC}"
echo -e "Deployment log: ${BLUE}$DEPLOY_LOG${NC}\n"
echo -e "Useful commands:"
echo -e "  Status: ${BLUE}systemctl status ariane${NC}"
echo -e "  Logs: ${BLUE}journalctl -u ariane -f${NC}"
echo -e "  Restart: ${BLUE}systemctl restart ariane${NC}"
echo -e "\nNext steps:"
echo -e "  1. Configure TLS: sudo bash scripts/server-ops/setup-tls.sh $ARIANE_SERVER_NAME <email>"
echo -e "  2. Set up automated backups"
echo -e "  3. Install monitoring tools"
echo -e "  4. Configure security settings"
echo -e "\nFor more information, see:"
echo -e "  ${BLUE}./README.md${NC} (in server-ops directory)"
echo -e "${BLUE}========================================${NC}\n"

echo "Deployment finished at $(date)" >> "$DEPLOY_LOG"
