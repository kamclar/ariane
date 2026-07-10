#!/bin/bash

# ARIANE Security Check Script
# Verifies security settings and best practices

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

WARNINGS=0
PASSED=0

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${BLUE}   ARIANE Security Audit${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"

# 1. SSH Configuration
echo -e "${YELLOW}[1] SSH Security${NC}"

if grep -q "^PermitRootLogin no" /etc/ssh/sshd_config; then
    echo -e "   ${GREEN}✓ Root login disabled${NC}"
    ((PASSED++))
else
    echo -e "   ${RED}✗ Root login not disabled${NC}"
    ((WARNINGS++))
fi

if grep -q "^PasswordAuthentication no" /etc/ssh/sshd_config; then
    echo -e "   ${GREEN}✓ Password authentication disabled${NC}"
    ((PASSED++))
else
    echo -e "   ${YELLOW}⚠ Password authentication enabled (use key-based auth)${NC}"
    ((WARNINGS++))
fi

# 2. Firewall Status
echo -e "\n${YELLOW}[2] Firewall Configuration${NC}"

if command -v ufw &> /dev/null; then
    UFW_STATUS=$(ufw status | head -1)
    if [[ "$UFW_STATUS" == *"active"* ]]; then
        echo -e "   ${GREEN}✓ UFW firewall is active${NC}"
        ((PASSED++))
    else
        echo -e "   ${RED}✗ UFW firewall is inactive${NC}"
        ((WARNINGS++))
    fi
    
    # Check open ports
    echo -e "   Open ports:"
    ufw status | grep ALLOW | sed 's/^/      /'
else
    echo -e "   ${YELLOW}⚠ UFW not installed${NC}"
fi

# 3. File Permissions
echo -e "\n${YELLOW}[3] File Permissions${NC}"

# Check SSH keys
if [ -d ~/.ssh ]; then
    SSH_PERMS=$(stat -c %a ~/.ssh)
    if [ "$SSH_PERMS" = "700" ]; then
        echo -e "   ${GREEN}✓ SSH directory permissions correct (700)${NC}"
        ((PASSED++))
    else
        echo -e "   ${RED}✗ SSH directory permissions incorrect ($SSH_PERMS, should be 700)${NC}"
        ((WARNINGS++))
    fi
fi

# Check /etc/shadow
SHADOW_PERMS=$(stat -c %a /etc/shadow)
if [ "$SHADOW_PERMS" = "640" ] || [ "$SHADOW_PERMS" = "600" ]; then
    echo -e "   ${GREEN}✓ /etc/shadow permissions correct${NC}"
    ((PASSED++))
else
    echo -e "   ${RED}✗ /etc/shadow permissions might be too open ($SHADOW_PERMS)${NC}"
    ((WARNINGS++))
fi

# 4. Sudo Configuration
echo -e "\n${YELLOW}[4] Sudo Access${NC}"

if sudo -l -U ubuntu | grep -q "NOPASSWD"; then
    echo -e "   ${RED}✗ NOPASSWD sudo detected (security risk)${NC}"
    ((WARNINGS++))
else
    echo -e "   ${GREEN}✓ Sudo password required${NC}"
    ((PASSED++))
fi

# 5. System Updates
echo -e "\n${YELLOW}[5] System Updates${NC}"

apt-get update -qq > /dev/null 2>&1
UPDATES=$(apt list --upgradable 2>/dev/null | wc -l)
if [ $UPDATES -le 1 ]; then
    echo -e "   ${GREEN}✓ System is up to date${NC}"
    ((PASSED++))
else
    echo -e "   ${YELLOW}⚠ $((UPDATES-1)) package(s) available for update${NC}"
    echo -e "   Run: ${BLUE}sudo apt update && sudo apt upgrade${NC}"
    ((WARNINGS++))
fi

# 6. Port Security
echo -e "\n${YELLOW}[6] Open Ports${NC}"

OPEN_PORTS=$(ss -tlnp 2>/dev/null | grep LISTEN | grep -v "127.0.0.1" | awk '{print $4}' | cut -d: -f2 | sort -u)
if [ -z "$OPEN_PORTS" ]; then
    echo -e "   ${GREEN}✓ No unexpected ports open${NC}"
    ((PASSED++))
else
    echo -e "   Open ports: $(echo $OPEN_PORTS | tr '\n' ' ')"
fi

# 7. SSL/TLS Certificates
echo -e "\n${YELLOW}[7] SSL/TLS Configuration${NC}"

if [ -f /etc/nginx/ssl/ariane.crt ]; then
    EXPIRY=$(openssl x509 -in /etc/nginx/ssl/ariane.crt -noout -enddate 2>/dev/null | cut -d= -f2)
    DAYS_LEFT=$(( ($(date -d "$EXPIRY" +%s) - $(date +%s)) / 86400 ))
    if [ $DAYS_LEFT -gt 30 ]; then
        echo -e "   ${GREEN}✓ SSL certificate valid (expires in $DAYS_LEFT days)${NC}"
        ((PASSED++))
    else
        echo -e "   ${RED}✗ SSL certificate expiring soon ($DAYS_LEFT days)${NC}"
        ((WARNINGS++))
    fi
else
    echo -e "   ${YELLOW}⚠ No SSL certificate found${NC}"
fi

# 8. Process Security
echo -e "\n${YELLOW}[8] Process Isolation${NC}"

if pgrep -f "uvicorn" > /dev/null; then
    ARIANE_USER=$(ps aux | grep "[u]vicorn" | awk '{print $1}' | head -1)
    if [ "$ARIANE_USER" != "root" ]; then
        echo -e "   ${GREEN}✓ ARIANE running as non-root user ($ARIANE_USER)${NC}"
        ((PASSED++))
    else
        echo -e "   ${RED}✗ ARIANE running as root (security risk)${NC}"
        ((WARNINGS++))
    fi
fi

# 9. Log Monitoring
echo -e "\n${YELLOW}[9] Log Monitoring${NC}"

if [ -f /var/log/auth.log ]; then
    FAILED_LOGINS=$(grep "Failed password" /var/log/auth.log 2>/dev/null | wc -l)
    if [ $FAILED_LOGINS -gt 10 ]; then
        echo -e "   ${RED}⚠ Multiple failed login attempts detected ($FAILED_LOGINS)${NC}"
        ((WARNINGS++))
    else
        echo -e "   ${GREEN}✓ No unusual login activity${NC}"
        ((PASSED++))
    fi
fi

# 10. Fail2ban Status
echo -e "\n${YELLOW}[10] Intrusion Protection${NC}"

if systemctl is-active --quiet fail2ban; then
    echo -e "   ${GREEN}✓ Fail2ban is running${NC}"
    ((PASSED++))
else
    echo -e "   ${YELLOW}⚠ Fail2ban is not running (consider installing)${NC}"
    ((WARNINGS++))
fi

# 11. Database Security (if applicable)
echo -e "\n${YELLOW}[11] Sensitive Data${NC}"

if grep -q "password\|API_KEY\|SECRET" /home/ubuntu/ariane/backend/config.py 2>/dev/null; then
    echo -e "   ${RED}✗ Potential sensitive data in config files${NC}"
    echo -e "   Consider using environment variables"
    ((WARNINGS++))
else
    echo -e "   ${GREEN}✓ No hardcoded secrets in accessible files${NC}"
    ((PASSED++))
fi

# Summary
echo -e "\n${BLUE}════════════════════════════════════════${NC}"
echo -e "Security Audit Results:"
echo -e "  ${GREEN}Passed: $PASSED${NC}"
echo -e "  ${YELLOW}Warnings: $WARNINGS${NC}"

if [ $WARNINGS -eq 0 ]; then
    echo -e "  ${GREEN}Overall: EXCELLENT${NC}"
elif [ $WARNINGS -le 2 ]; then
    echo -e "  ${YELLOW}Overall: GOOD${NC}"
else
    echo -e "  ${RED}Overall: NEEDS ATTENTION${NC}"
fi

echo -e "\n${YELLOW}Recommendations:${NC}"
echo -e "  1. Enable SSL/TLS certificates (Let's Encrypt)"
echo -e "  2. Configure automated updates (unattended-upgrades)"
echo -e "  3. Install and configure Fail2ban"
echo -e "  4. Use environment variables for secrets"
echo -e "  5. Regular security audits and patches"
echo -e "  6. Keep SSH key secure and backed up"
echo -e "  7. Monitor logs regularly"
echo -e "${BLUE}════════════════════════════════════════${NC}\n"
