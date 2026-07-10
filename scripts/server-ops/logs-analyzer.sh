#!/bin/bash

# ARIANE Logs Analyzer Script
# Analyzes recent logs for errors and warnings

LINES="${1:-100}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
NC='\033[0m'

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${BLUE}   ARIANE Logs Analyzer${NC}"
echo -e "${BLUE}   Analysis of last $LINES lines${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"

# 1. Get raw logs
echo -e "${YELLOW}[1] Fetching logs from journalctl...${NC}\n"
LOGS=$(journalctl -u ariane -n $LINES 2>/dev/null)

if [ -z "$LOGS" ]; then
    echo -e "${RED}No logs found for ariane service${NC}"
    exit 1
fi

# 2. Count different severity levels
echo -e "${YELLOW}[2] Severity Analysis${NC}"
ERROR_COUNT=$(echo "$LOGS" | grep -ic "error")
WARNING_COUNT=$(echo "$LOGS" | grep -ic "warning\|warn")
EXCEPTION_COUNT=$(echo "$LOGS" | grep -ic "exception")
CRITICAL_COUNT=$(echo "$LOGS" | grep -ic "critical\|fatal")

echo -e "   ${RED}Errors: $ERROR_COUNT${NC}"
echo -e "   ${YELLOW}Warnings: $WARNING_COUNT${NC}"
echo -e "   ${MAGENTA}Exceptions: $EXCEPTION_COUNT${NC}"
echo -e "   ${RED}Critical: $CRITICAL_COUNT${NC}"

# 3. Show errors
if [ $ERROR_COUNT -gt 0 ]; then
    echo -e "\n${YELLOW}[3] Recent Errors${NC}"
    echo "$LOGS" | grep -i "error" | head -5 | nl | sed 's/^/   /'
fi

# 4. Show exceptions
if [ $EXCEPTION_COUNT -gt 0 ]; then
    echo -e "\n${YELLOW}[4] Recent Exceptions${NC}"
    echo "$LOGS" | grep -i "exception" | head -5 | nl | sed 's/^/   /'
fi

# 5. Show warnings
if [ $WARNING_COUNT -gt 0 ]; then
    echo -e "\n${YELLOW}[5] Recent Warnings${NC}"
    echo "$LOGS" | grep -i "warning\|warn" | head -5 | nl | sed 's/^/   /'
fi

# 6. API endpoints analysis
echo -e "\n${YELLOW}[6] API Endpoints (last 10 requests)${NC}"
ENDPOINTS=$(echo "$LOGS" | grep -o 'GET\|POST\|PUT\|DELETE' | head -10 || echo "No endpoints found")
if [ ! -z "$ENDPOINTS" ]; then
    echo "$ENDPOINTS" | sort | uniq -c | sed 's/^/   /'
else
    echo -e "   No endpoint logs found"
fi

# 7. Failed requests
echo -e "\n${YELLOW}[7] Failed Requests (4xx, 5xx)${NC}"
FAILED=$(echo "$LOGS" | grep -E ' [45][0-9]{2} ' | wc -l)
echo -e "   Total failed: $FAILED"
if [ $FAILED -gt 0 ]; then
    echo "$LOGS" | grep -E ' [45][0-9]{2} ' | head -5 | sed 's/^/   /'
fi

# 8. Performance analysis
echo -e "\n${YELLOW}[8] Performance Indicators${NC}"
SLOW_REQUESTS=$(echo "$LOGS" | grep -i "slow\|timeout\|took" | wc -l)
if [ $SLOW_REQUESTS -gt 0 ]; then
    echo -e "   ${YELLOW}Slow requests detected: $SLOW_REQUESTS${NC}"
fi

# 9. Disk/Memory issues
echo -e "\n${YELLOW}[9] Resource Issues${NC}"
DISK_ISSUES=$(echo "$LOGS" | grep -ic "disk\|space\|storage")
MEM_ISSUES=$(echo "$LOGS" | grep -ic "memory\|ram\|out of")
if [ $DISK_ISSUES -gt 0 ]; then
    echo -e "   ${RED}Disk issues: $DISK_ISSUES${NC}"
fi
if [ $MEM_ISSUES -gt 0 ]; then
    echo -e "   ${RED}Memory issues: $MEM_ISSUES${NC}"
fi

# 10. Database/External service issues
echo -e "\n${YELLOW}[10] External Service Issues${NC}"
DB_ISSUES=$(echo "$LOGS" | grep -ic "database\|db\|connection refused")
API_ISSUES=$(echo "$LOGS" | grep -ic "api error\|external\|unreachable")
if [ $DB_ISSUES -gt 0 ]; then
    echo -e "   ${RED}Database issues: $DB_ISSUES${NC}"
fi
if [ $API_ISSUES -gt 0 ]; then
    echo -e "   ${RED}API/External issues: $API_ISSUES${NC}"
fi

# 11. Raw logs (optional)
echo -e "\n${YELLOW}[11] Raw Logs (last $LINES entries)${NC}"
echo "$LOGS" | tail -20 | sed 's/^/   /'

# Summary
echo -e "\n${BLUE}════════════════════════════════════════${NC}"
TOTAL_ISSUES=$((ERROR_COUNT + WARNING_COUNT + EXCEPTION_COUNT + CRITICAL_COUNT))
if [ $TOTAL_ISSUES -eq 0 ]; then
    echo -e "${GREEN}✓ No issues detected${NC}"
else
    echo -e "${YELLOW}⚠ Found $TOTAL_ISSUES total issues${NC}"
fi
echo -e "Analysis completed at $(date)"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"
