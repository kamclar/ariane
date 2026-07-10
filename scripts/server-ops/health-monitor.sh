#!/bin/bash

# ARIANE Health Monitor Script
# Continuously monitors application health

# Configuration
INTERVAL="${1:-30}"  # Check interval in seconds
CHECK_COUNT=0
ERROR_COUNT=0

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

clear

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${BLUE}   ARIANE Health Monitor${NC}"
echo -e "${BLUE}   Interval: ${INTERVAL}s (Press Ctrl+C to exit)${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"

# Trap Ctrl+C
trap 'echo -e "\n${YELLOW}Monitoring stopped${NC}"; exit 0' INT

# Main monitoring loop
while true; do
    ((CHECK_COUNT++))
    TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
    
    # Clear previous output
    clear
    
    echo -e "${BLUE}════════════════════════════════════════${NC}"
    echo -e "${BLUE}   ARIANE Health Monitor${NC}"
    echo -e "${BLUE}   Check #${CHECK_COUNT} | Errors: ${ERROR_COUNT}${NC}"
    echo -e "${BLUE}════════════════════════════════════════${NC}\n"
    
    # 1. Service Status
    echo -e "${YELLOW}[Service Status]${NC}"
    if systemctl is-active --quiet ariane; then
        echo -e "   ${GREEN}✓ ARIANE service: RUNNING${NC}"
    else
        echo -e "   ${RED}✗ ARIANE service: STOPPED${NC}"
        ((ERROR_COUNT++))
    fi
    
    if systemctl is-active --quiet nginx; then
        echo -e "   ${GREEN}✓ Nginx: RUNNING${NC}"
    else
        echo -e "   ${RED}✗ Nginx: STOPPED${NC}"
        ((ERROR_COUNT++))
    fi
    
    # 2. API Health
    echo -e "\n${YELLOW}[API Health]${NC}"
    if command -v curl &> /dev/null; then
        HEALTH_RESPONSE=$(curl -s -m 5 http://localhost:8000/api/health 2>/dev/null)
        
        if echo "$HEALTH_RESPONSE" | grep -q '"status":"ok"'; then
            echo -e "   ${GREEN}✓ API Health: OK${NC}"
            
            # Extract additional info
            TABLE4=$(echo "$HEALTH_RESPONSE" | grep -o '"table4":[^,}]*' | cut -d':' -f2)
            TABLE9=$(echo "$HEALTH_RESPONSE" | grep -o '"table9":[^,}]*' | cut -d':' -f2)
            ST7=$(echo "$HEALTH_RESPONSE" | grep -o '"st7":[^,}]*' | cut -d':' -f2)
            echo -e "   Data files: Table4=$TABLE4, Table9=$TABLE9, ST7=$ST7"
        else
            echo -e "   ${RED}✗ API Health: FAILED${NC}"
            ((ERROR_COUNT++))
        fi
    else
        echo -e "   ${YELLOW}⚠ curl not available${NC}"
    fi
    
    # 3. System Resources
    echo -e "\n${YELLOW}[System Resources]${NC}"
    MEM_USED=$(free | awk 'NR==2 {printf "%.1f", ($3/$2)*100}')
    LOAD=$(uptime | awk -F'load average:' '{print $2}' | awk '{print $1}')
    DISK=$(df / | tail -1 | awk '{print $5}')
    
    echo -e "   Memory: ${MEM_USED}% | Load: ${LOAD} | Disk: ${DISK}"
    
    # 4. Process Check
    echo -e "\n${YELLOW}[Process Check]${NC}"
    if pgrep -f "uvicorn" > /dev/null; then
        PID=$(pgrep -f "uvicorn" | head -1)
        # Get process memory usage
        PROC_MEM=$(ps -p $PID -o %mem= 2>/dev/null | awk '{print $1}')
        PROC_CPU=$(ps -p $PID -o %cpu= 2>/dev/null | awk '{print $1}')
        echo -e "   ${GREEN}✓ FastAPI running (PID: $PID)${NC}"
        echo -e "   CPU: ${PROC_CPU}% | Memory: ${PROC_MEM}%"
    else
        echo -e "   ${RED}✗ FastAPI process not found${NC}"
        ((ERROR_COUNT++))
    fi
    
    # 5. Nginx Check
    echo -e "\n${YELLOW}[Nginx Check]${NC}"
    NGINX_WORKERS=$(ps aux | grep -c "[n]ginx: worker" || echo "0")
    if [ "$NGINX_WORKERS" -gt 0 ]; then
        echo -e "   ${GREEN}✓ Nginx workers running: $NGINX_WORKERS${NC}"
    else
        echo -e "   ${YELLOW}⚠ No Nginx workers found${NC}"
    fi
    
    # 6. Recent Errors
    echo -e "\n${YELLOW}[Recent Issues]${NC}"
    ERROR_LINES=$(journalctl -u ariane -n 50 2>/dev/null | grep -i "error\|exception\|failed" | head -3 || echo "")
    if [ -z "$ERROR_LINES" ]; then
        echo -e "   ${GREEN}✓ No recent errors${NC}"
    else
        echo -e "   ${RED}Recent errors:${NC}"
        echo "$ERROR_LINES" | sed 's/^/   /'
    fi
    
    # Footer
    echo -e "\n${BLUE}════════════════════════════════════════${NC}"
    echo -e "Last update: ${TIMESTAMP}"
    echo -e "Next check in: ${INTERVAL}s"
    echo -e "${BLUE}════════════════════════════════════════${NC}"
    
    # Wait for next check
    sleep "$INTERVAL"
done
