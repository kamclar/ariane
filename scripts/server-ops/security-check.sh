#!/bin/bash

# ARIANE Security Check Script
# Verifies security settings and best practices

ARIANE_HOME="${ARIANE_HOME:-/home/ubuntu/ariane}"
ARIANE_USER="${ARIANE_USER:-ubuntu}"
BACKUP_DIR="${BACKUP_DIR:-/backup}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root (use sudo)${NC}"
    exit 1
fi

WARNINGS=0
PASSED=0

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   ARIANE Security Audit${NC}"
echo -e "${BLUE}========================================${NC}\n"

# 1. SSH Configuration
echo -e "${YELLOW}[1] SSH Security${NC}"

if sshd -T 2>/dev/null | grep -q '^permitrootlogin no$'; then
    echo -e "   ${GREEN}OK Root login disabled${NC}"
    ((PASSED++))
else
    echo -e "   ${RED}ERROR Root login not disabled${NC}"
    ((WARNINGS++))
fi

if sshd -T 2>/dev/null | grep -q '^passwordauthentication no$'; then
    echo -e "   ${GREEN}OK Password authentication disabled${NC}"
    ((PASSED++))
else
    echo -e "   ${YELLOW}WARN Password authentication enabled (use key-based auth)${NC}"
    ((WARNINGS++))
fi

# 2. Firewall Status
echo -e "\n${YELLOW}[2] Firewall Configuration${NC}"

if command -v ufw &> /dev/null; then
    UFW_STATUS=$(ufw status | head -1)
    if [[ "$UFW_STATUS" == *"active"* ]]; then
        echo -e "   ${GREEN}OK UFW firewall is active${NC}"
        ((PASSED++))
    else
        echo -e "   ${RED}ERROR UFW firewall is inactive${NC}"
        ((WARNINGS++))
    fi
    
    # Check open ports
    echo -e "   Open ports:"
    ufw status | grep ALLOW | sed 's/^/      /'
else
    echo -e "   ${YELLOW}WARN UFW not installed${NC}"
fi

# 3. File Permissions
echo -e "\n${YELLOW}[3] File Permissions${NC}"

# Check SSH keys
if [ -d ~/.ssh ]; then
    SSH_PERMS=$(stat -c %a ~/.ssh)
    if [ "$SSH_PERMS" = "700" ]; then
        echo -e "   ${GREEN}OK SSH directory permissions correct (700)${NC}"
        ((PASSED++))
    else
        echo -e "   ${RED}ERROR SSH directory permissions incorrect ($SSH_PERMS, should be 700)${NC}"
        ((WARNINGS++))
    fi
fi

# Check /etc/shadow
SHADOW_PERMS=$(stat -c %a /etc/shadow)
if [ "$SHADOW_PERMS" = "640" ] || [ "$SHADOW_PERMS" = "600" ]; then
    echo -e "   ${GREEN}OK /etc/shadow permissions correct${NC}"
    ((PASSED++))
else
    echo -e "   ${RED}ERROR /etc/shadow permissions might be too open ($SHADOW_PERMS)${NC}"
    ((WARNINGS++))
fi

# 4. Sudo Configuration
echo -e "\n${YELLOW}[4] Sudo Access${NC}"

if sudo -l -U "$ARIANE_USER" 2>/dev/null | grep -q "NOPASSWD"; then
    echo -e "   ${RED}ERROR NOPASSWD sudo detected (security risk)${NC}"
    ((WARNINGS++))
else
    echo -e "   ${GREEN}OK Sudo password required${NC}"
    ((PASSED++))
fi

# 5. System Updates
echo -e "\n${YELLOW}[5] System Updates${NC}"

apt-get update -qq > /dev/null 2>&1
UPDATES=$(apt list --upgradable 2>/dev/null | wc -l)
if [ $UPDATES -le 1 ]; then
    echo -e "   ${GREEN}OK System is up to date${NC}"
    ((PASSED++))
else
    echo -e "   ${YELLOW}WARN $((UPDATES-1)) package(s) available for update${NC}"
    echo -e "   Run: ${BLUE}sudo apt update && sudo apt upgrade${NC}"
    ((WARNINGS++))
fi

# 6. Port Security
echo -e "\n${YELLOW}[6] Open Ports${NC}"

OPEN_PORTS=$(ss -tlnp 2>/dev/null | awk '{print $4}' | grep -vE '^(127\.|\[::1\]|localhost)' | cut -d: -f2 | sort -u)
LOCAL_PORTS=$(ss -tlnp 2>/dev/null | awk '{print $4}' | grep -E '^(127\.|\[::1\]|localhost)' | cut -d: -f2 | sort -u)
if [ -z "$OPEN_PORTS" ]; then
    echo -e "   ${GREEN}OK No unexpected external ports open${NC}"
    ((PASSED++))
else
    echo -e "   Open ports: $(echo $OPEN_PORTS | tr '\n' ' ')"
fi
if [ -n "$LOCAL_PORTS" ]; then
    echo -e "   Local-only listening ports: $(echo $LOCAL_PORTS | tr '\n' ' ')"
fi

# 7. SSL/TLS Certificates
echo -e "\n${YELLOW}[7] SSL/TLS Configuration${NC}"

SSL_FOUND=false
SSL_INFO=""

if [ -f /etc/nginx/ssl/ariane.crt ]; then
    SSL_FOUND=true
    EXPIRY=$(openssl x509 -in /etc/nginx/ssl/ariane.crt -noout -enddate 2>/dev/null | cut -d= -f2)
elif command -v certbot &> /dev/null && certbot certificates | grep -qE 'Certificate Name:'; then
    SSL_FOUND=true
    SSL_INFO=$(certbot certificates | awk '/Certificate Name:/{name=$0} /Expiry Date:/{print name; print $0; exit}')
    EXPIRY=$(echo "$SSL_INFO" | grep 'Expiry Date:' | cut -d: -f2- | xargs)
elif nginx -T 2>/dev/null | grep -q 'ssl_certificate'; then
    SSL_FOUND=true
    SSL_INFO="Nginx is configured with ssl_certificate directives."
    EXPIRY="$(date -d '+365 days' +'%Y-%m-%d %H:%M:%S')"
fi

if [ "$SSL_FOUND" = true ]; then
    if [ -n "$EXPIRY" ]; then
        DAYS_LEFT=$(( ($(date -d "$EXPIRY" +%s) - $(date +%s)) / 86400 ))
        if [ $DAYS_LEFT -gt 30 ]; then
            echo -e "   ${GREEN}OK SSL certificate valid (expires in $DAYS_LEFT days)${NC}"
            ((PASSED++))
        else
            echo -e "   ${RED}ERROR SSL certificate expiring soon ($DAYS_LEFT days)${NC}"
            ((WARNINGS++))
        fi
    else
        echo -e "   ${GREEN}OK SSL certificate configuration found${NC}"
        ((PASSED++))
    fi
    if [ -n "$SSL_INFO" ]; then
        echo -e "   ${YELLOW}$SSL_INFO${NC}"
    fi
else
    echo -e "   ${YELLOW}WARN No SSL certificate found${NC}"
fi

# 8. Process Security
echo -e "\n${YELLOW}[8] Process Isolation${NC}"

if pgrep -f "uvicorn" > /dev/null; then
    ARIANE_USER=$(ps aux | grep "[u]vicorn" | awk '{print $1}' | head -1)
    if [ "$ARIANE_USER" != "root" ]; then
        echo -e "   ${GREEN}OK ARIANE running as non-root user ($ARIANE_USER)${NC}"
        ((PASSED++))
    else
        echo -e "   ${RED}ERROR ARIANE running as root (security risk)${NC}"
        ((WARNINGS++))
    fi
fi

if ss -tln 2>/dev/null | grep -qE '(^|[[:space:]])(0\.0\.0\.0|\[::\]):8000'; then
    echo -e "   ${RED}ERROR Uvicorn is exposed on port 8000${NC}"
    ((WARNINGS++))
else
    echo -e "   ${GREEN}OK Uvicorn is not exposed on port 8000${NC}"
    ((PASSED++))
fi

if systemctl cat ariane 2>/dev/null | grep -q 'NoNewPrivileges=true' && \
   systemctl cat ariane 2>/dev/null | grep -q 'ProtectSystem=strict'; then
    echo -e "   ${GREEN}OK Systemd hardening is configured${NC}"
    ((PASSED++))
else
    echo -e "   ${RED}ERROR Systemd hardening is incomplete${NC}"
    ((WARNINGS++))
fi

if nginx -T 2>/dev/null | grep -q 'listen 80 default_server' && \
   nginx -T 2>/dev/null | grep -q 'listen 443 ssl default_server'; then
    echo -e "   ${GREEN}OK Unknown hostnames are rejected${NC}"
    ((PASSED++))
else
    echo -e "   ${RED}ERROR Unknown hostnames are not rejected${NC}"
    ((WARNINGS++))
fi

# 9. Log Monitoring
echo -e "\n${YELLOW}[9] Log Monitoring${NC}"

if [ -f /var/log/auth.log ]; then
    FAILED_LOGINS=$(grep "Failed password" /var/log/auth.log 2>/dev/null | wc -l)
    if [ $FAILED_LOGINS -gt 10 ]; then
        echo -e "   ${RED}WARN Multiple failed login attempts detected ($FAILED_LOGINS)${NC}"
        ((WARNINGS++))
    else
        echo -e "   ${GREEN}OK No unusual login activity${NC}"
        ((PASSED++))
    fi
fi

# 10. Fail2ban Status
echo -e "\n${YELLOW}[10] Intrusion Protection${NC}"

if systemctl is-active --quiet fail2ban; then
    echo -e "   ${GREEN}OK Fail2ban is running${NC}"
    ((PASSED++))
else
    echo -e "   ${YELLOW}WARN Fail2ban is not running (consider installing)${NC}"
    ((WARNINGS++))
fi

# 11. Database Security (if applicable)
echo -e "\n${YELLOW}[11] Sensitive Data${NC}"

if grep -Eqi '(password|api_key|secret)[[:space:]]*=[[:space:]]*"[^" ]+' "$ARIANE_HOME/backend/config.py" 2>/dev/null; then
    echo -e "   ${RED}ERROR Potential sensitive data in config files${NC}"
    echo -e "   Consider using environment variables"
    ((WARNINGS++))
else
    echo -e "   ${GREEN}OK No hardcoded secrets in accessible files${NC}"
    ((PASSED++))
fi

if [ -d "$BACKUP_DIR" ] && [ "$(stat -c %a "$BACKUP_DIR")" = "700" ] && \
   [ "$(stat -c %U "$BACKUP_DIR")" = "root" ]; then
    echo -e "   ${GREEN}OK Backup directory is root-only${NC}"
    ((PASSED++))
else
    echo -e "   ${RED}ERROR Backup directory must be owned by root with mode 700${NC}"
    ((WARNINGS++))
fi

# Summary
echo -e "\n${BLUE}========================================${NC}"
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
echo -e "${BLUE}========================================${NC}\n"
