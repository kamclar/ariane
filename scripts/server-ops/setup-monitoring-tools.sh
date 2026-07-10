#!/bin/bash

# ARIANE Monitoring Tools Setup Script
# Installs useful monitoring and debugging tools

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   Monitoring Tools Setup${NC}"
echo -e "${BLUE}========================================${NC}\n"

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Error: This script must be run as root (use sudo)${NC}"
    exit 1
fi

# Update package lists
echo -e "${YELLOW}[1] Updating package lists...${NC}"
apt-get update -qq
echo -e "${GREEN}OK Package lists updated${NC}"

# Tools to install
TOOLS=(
    "htop"
    "iotop"
    "net-tools"
    "nethogs"
    "lnav"
    "curl"
    "wget"
    "vim"
    "git"
    "jq"
    "unzip"
    "gzip"
    "tar"
)

# Install tools
echo -e "\n${YELLOW}[2] Installing tools...${NC}"
INSTALLED=0
FAILED=0

for tool in "${TOOLS[@]}"; do
    echo -ne "   Installing $tool... "
    if apt-get install -y -qq "$tool" 2>/dev/null; then
        echo -e "${GREEN}OK${NC}"
        ((INSTALLED++))
    else
        echo -e "${RED}ERROR${NC}"
        ((FAILED++))
    fi
done

echo -e "\n${GREEN}Installed: $INSTALLED tools${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED tools${NC}"
fi

# Install additional useful packages
echo -e "\n${YELLOW}[3] Installing additional utilities...${NC}"

# tmux (for persistent sessions)
echo -ne "   Installing tmux... "
if apt-get install -y -qq tmux 2>/dev/null; then
    echo -e "${GREEN}OK${NC}"
fi

# ranger (file manager)
echo -ne "   Installing ranger... "
if apt-get install -y -qq ranger 2>/dev/null; then
    echo -e "${GREEN}OK${NC}"
fi

# Install Docker (optional)
echo -e "\n${YELLOW}[4] Docker (optional)${NC}"
echo "   Install Docker CLI? (y/n)"
read -r response
if [[ "$response" == "y" || "$response" == "Y" ]]; then
    apt-get install -y -qq docker.io 2>/dev/null || true
    echo -e "${GREEN}OK Docker installed${NC}"
fi

# Install ctop (Docker container monitoring)
echo -e "\n${YELLOW}[5] Container monitoring tools${NC}"
if command -v docker &> /dev/null; then
    echo "   Installing ctop..."
    if command -v ctop &> /dev/null; then
        echo -e "   ${GREEN}OK ctop already installed${NC}"
    else
        echo -e "   ${YELLOW}Note: Install ctop manually with: docker run --rm -ti --name=ctop --volume /var/run/docker.sock:/var/run/docker.sock:ro quay.io/vektorlab/ctop:latest${NC}"
    fi
fi

# Create helpful alias file
echo -e "\n${YELLOW}[6] Creating helpful aliases...${NC}"
cat > /tmp/ariane-aliases.sh << 'EOF'
# ARIANE Monitoring Aliases
alias ariane-status='systemctl status ariane'
alias ariane-logs='journalctl -u ariane -f'
alias ariane-restart='sudo systemctl restart ariane'
alias ariane-stop='sudo systemctl stop ariane'
alias ariane-start='sudo systemctl start ariane'

# System monitoring
alias mem-top='free -h && echo "---" && top -b -n 1 | head -20'
alias disk-usage='du -sh /* 2>/dev/null | sort -h'
alias port-check='sudo netstat -tlnp | grep LISTEN'
alias process-top='ps aux --sort=-%cpu | head -20'

# Logs
alias logs-errors='journalctl -u ariane -p err'
alias logs-warnings='journalctl -u ariane -p warning'
alias logs-today='journalctl -u ariane --since=today'

# Network
alias net-active='ss -tlnp'
alias net-connections='ss -tpn | grep ESTAB'
EOF

if [ -f /home/ubuntu/.bashrc ]; then
    echo "   Adding aliases to /home/ubuntu/.bashrc..."
    if ! grep -q "ARIANE Monitoring Aliases" /home/ubuntu/.bashrc; then
        cat /tmp/ariane-aliases.sh >> /home/ubuntu/.bashrc
        chown ubuntu:ubuntu /home/ubuntu/.bashrc
        echo -e "   ${GREEN}OK Aliases added${NC}"
    fi
fi

# Create monitoring dashboard script
echo -e "\n${YELLOW}[7] Creating dashboard script...${NC}"
cat > /usr/local/bin/ariane-dashboard << 'EOF'
#!/bin/bash
# ARIANE Dashboard

while true; do
    clear
    echo "==========================================="
    echo "     ARIANE Monitoring Dashboard"
    echo "==========================================="
    echo ""
    echo "SERVICE STATUS:"
    systemctl status ariane --no-pager | grep "Active" | awk '{print "  " $0}'
    echo ""
    echo "RESOURCE USAGE:"
    echo "  $(free -h | awk 'NR==2 {print "Memory: " $3 " / " $2}')"
    echo "  $(df -h / | awk 'NR==2 {print "Disk: " $3 " / " $2}')"
    echo ""
    echo "PROCESSES:"
    ps aux | grep "[u]vicorn" | awk '{print "  FastAPI: PID " $2 " (" $3 "% CPU, " $4 "% MEM)"}'
    ps aux | grep "[n]ginx: worker" | wc -l | awk '{print "  Nginx workers: " $1}'
    echo ""
    echo "API STATUS:"
    curl --fail --silent --max-time 5 http://127.0.0.1:8000/api/health | jq -r '"  Status: " + .status' 2>/dev/null || echo "  (API not responding)"
    echo ""
    echo "==========================================="
    echo "Updated: $(date)"
    echo "Press Ctrl+C to exit, refreshes every 5s"
    echo "==========================================="
    sleep 5
done
EOF

chmod +x /usr/local/bin/ariane-dashboard
echo -e "   ${GREEN}OK Dashboard script created${NC}"
echo -e "   Run with: ${BLUE}ariane-dashboard${NC}"

# Summary
echo -e "\n${BLUE}========================================${NC}"
echo -e "${GREEN}OK Monitoring tools setup completed!${NC}\n"
echo -e "Useful commands:"
echo -e "  ${BLUE}htop${NC} - Interactive process monitor"
echo -e "  ${BLUE}iotop${NC} - I/O monitoring (requires sudo)"
echo -e "  ${BLUE}nethogs${NC} - Network monitoring by process"
echo -e "  ${BLUE}lnav${NC} - Log file navigator"
echo -e "  ${BLUE}ariane-dashboard${NC} - Custom dashboard"
echo -e "\nUseful aliases:"
echo -e "  ${BLUE}ariane-logs${NC} - Follow ARIANE logs"
echo -e "  ${BLUE}ariane-status${NC} - Check service status"
echo -e "  ${BLUE}mem-top${NC} - Show memory usage"
echo -e "  ${BLUE}disk-usage${NC} - Show disk usage"
echo -e "\nSource updated aliases with: ${BLUE}source ~/.bashrc${NC}"
echo -e "${BLUE}========================================${NC}\n"
