#!/bin/bash

# ARIANE Status Check Script
# Quick overview of application health and server status

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${BLUE}   ARIANE Status Check${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"

# 1. Service Status
echo -e "${YELLOW}[1] Service Status${NC}"
if systemctl is-active --quiet ariane; then
    echo -e "   ${GREEN}✓ ARIANE service: RUNNING${NC}"
else
    echo -e "   ${RED}✗ ARIANE service: STOPPED${NC}"
fi

if systemctl is-active --quiet nginx; then
    echo -e "   ${GREEN}✓ Nginx: RUNNING${NC}"
else
    echo -e "   ${RED}✗ Nginx: STOPPED${NC}"
fi

# 2. Health Check
echo -e "\n${YELLOW}[2] API Health Check${NC}"
if command -v curl &> /dev/null; then
    HEALTH=$(curl -s http://localhost:8000/api/health 2>/dev/null | grep -o '"status":"[^"]*"' | cut -d'"' -f4 || echo "UNKNOWN")
    if [ "$HEALTH" = "ok" ]; then
        echo -e "   ${GREEN}✓ API Health: OK${NC}"
    else
        echo -e "   ${RED}✗ API Health: $HEALTH${NC}"
    fi
else
    echo -e "   ${YELLOW}⚠ curl not available, skipping${NC}"
fi

# 3. Process Status
echo -e "\n${YELLOW}[3] Process Status${NC}"
if pgrep -f "uvicorn" > /dev/null; then
    PID=$(pgrep -f "uvicorn" | head -1)
    echo -e "   ${GREEN}✓ FastAPI process running (PID: $PID)${NC}"
else
    echo -e "   ${RED}✗ FastAPI process not found${NC}"
fi

# 4. Memory & CPU
echo -e "\n${YELLOW}[4] System Resources${NC}"
TOTAL_MEM=$(free -h | awk 'NR==2 {print $2}')
USED_MEM=$(free -h | awk 'NR==2 {print $3}')
MEM_PERCENT=$(free | awk 'NR==2 {printf "%.1f", ($3/$2)*100}')
echo -e "   Memory: ${USED_MEM} / ${TOTAL_MEM} (${MEM_PERCENT}%)"

# CPU
CPU_COUNT=$(nproc)
LOAD=$(uptime | awk -F'load average:' '{print $2}' | awk '{print $1}')
echo -e "   CPU Cores: ${CPU_COUNT}"
echo -e "   Load Average: ${LOAD}"

# 5. Disk Space
echo -e "\n${YELLOW}[5] Disk Space${NC}"
DISK_OUTPUT=$(df -h / | tail -1)
DISK_USED=$(echo $DISK_OUTPUT | awk '{print $3}')
DISK_TOTAL=$(echo $DISK_OUTPUT | awk '{print $2}')
DISK_PERCENT=$(echo $DISK_OUTPUT | awk '{print $5}')
echo -e "   Root: ${DISK_USED} / ${DISK_TOTAL} (${DISK_PERCENT})"

# Check data directory
if [ -d "/home/ubuntu/ariane/backend/data" ]; then
    DATA_SIZE=$(du -sh /home/ubuntu/ariane/backend/data 2>/dev/null | awk '{print $1}')
    echo -e "   Data dir: ${DATA_SIZE}"
fi

# 6. Recent Errors
echo -e "\n${YELLOW}[6] Recent Errors (last 10)${NC}"
ERROR_COUNT=$(journalctl -u ariane -n 100 2>/dev/null | grep -i "error\|exception\|failed\|critical" | wc -l || echo "0")
if [ "$ERROR_COUNT" -gt 0 ]; then
    echo -e "   ${RED}Found ${ERROR_COUNT} error(s)${NC}"
    journalctl -u ariane -n 100 2>/dev/null | grep -i "error\|exception\|failed\|critical" | tail -5 | sed 's/^/   /'
else
    echo -e "   ${GREEN}No errors found${NC}"
fi

# 7. Port Status
echo -e "\n${YELLOW}[7] Port Status${NC}"
if command -v ss &> /dev/null; then
    ss -tlnp 2>/dev/null | grep -E ":(8000|80|443)" | sed 's/^/   /' || echo "   (no services found on expected ports)"
fi

# 8. Last Restart
echo -e "\n${YELLOW}[8] Service Uptime${NC}"
if systemctl is-active --quiet ariane; then
    UPTIME=$(systemctl show ariane --property=ActiveEnterTimestamp --value)
    echo -e "   Started: ${UPTIME}"
fi

echo -e "\n${BLUE}════════════════════════════════════════${NC}"
echo -e "${GREEN}Check completed at $(date)${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"
