#!/bin/bash

# ARIANE Restart Script
# Safely restarts ARIANE and Nginx

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${BLUE}   ARIANE Restart Script${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root (use sudo)${NC}"
    exit 1
fi

# 1. Stop services
echo -e "${YELLOW}[1] Stopping services...${NC}"
systemctl stop ariane || true
sleep 1
systemctl stop nginx || true
sleep 1
echo -e "${GREEN}✓ Services stopped${NC}"

# 2. Wait a moment
echo -e "${YELLOW}[2] Waiting 3 seconds...${NC}"
sleep 3

# 3. Start services
echo -e "${YELLOW}[3] Starting services...${NC}"
systemctl start ariane
sleep 2
systemctl start nginx
sleep 2
echo -e "${GREEN}✓ Services started${NC}"

# 4. Check status
echo -e "\n${YELLOW}[4] Checking status...${NC}"
if systemctl is-active --quiet ariane; then
    echo -e "   ${GREEN}✓ ARIANE: RUNNING${NC}"
else
    echo -e "   ${RED}✗ ARIANE: FAILED TO START${NC}"
fi

if systemctl is-active --quiet nginx; then
    echo -e "   ${GREEN}✓ Nginx: RUNNING${NC}"
else
    echo -e "   ${RED}✗ Nginx: FAILED TO START${NC}"
fi

# 5. Health check
echo -e "\n${YELLOW}[5] Health check...${NC}"
sleep 2
if command -v curl &> /dev/null; then
    for i in {1..10}; do
        if curl -s http://localhost:8000/api/health > /dev/null 2>&1; then
            echo -e "   ${GREEN}✓ API responding${NC}"
            break
        fi
        if [ $i -lt 10 ]; then
            echo -e "   Waiting... attempt $i/10"
            sleep 1
        fi
    done
else
    echo -e "   ${YELLOW}⚠ curl not available, skipping health check${NC}"
fi

echo -e "\n${BLUE}════════════════════════════════════════${NC}"
echo -e "${GREEN}Restart completed at $(date)${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"
