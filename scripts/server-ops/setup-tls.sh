#!/bin/bash

# Configure a Let's Encrypt certificate for the ARIANE Nginx site.

set -euo pipefail

SERVER_NAME="${1:-}"
EMAIL="${2:-}"

if [ "$EUID" -ne 0 ]; then
    echo "Run this script as root" >&2
    exit 1
fi
if [ -z "$SERVER_NAME" ] || [ -z "$EMAIL" ]; then
    echo "Usage: $0 <server-name> <email>" >&2
    exit 2
fi
if ! [[ "$SERVER_NAME" =~ ^[A-Za-z0-9.-]+$ ]] || [[ "$SERVER_NAME" == *..* ]]; then
    echo "Invalid server name" >&2
    exit 2
fi

apt-get update -qq
apt-get install -y -qq certbot python3-certbot-nginx
if ! grep -Fq "server_name $SERVER_NAME;" /etc/nginx/sites-available/ariane; then
    echo "Nginx is not configured for $SERVER_NAME" >&2
    echo "Set ARIANE_SERVER_NAME during deployment" >&2
    exit 1
fi
nginx -t
systemctl reload nginx
certbot --nginx --non-interactive --agree-tos --redirect --hsts --staple-ocsp \
    --email "$EMAIL" -d "$SERVER_NAME"

cat > /etc/nginx/conf.d/ariane-default-tls.conf <<EOF
# Reject direct IP access and unknown TLS hostnames.
server {
    listen 443 ssl default_server;
    listen [::]:443 ssl default_server;
    server_name _;
    ssl_certificate /etc/letsencrypt/live/$SERVER_NAME/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/$SERVER_NAME/privkey.pem;
    return 444;
}
EOF

nginx -t
systemctl reload nginx
systemctl enable --now certbot.timer
certbot renew --dry-run
