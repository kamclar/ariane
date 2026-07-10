# ARIANE server operations

Production scripts for an Ubuntu server with systemd and Nginx.

## Initial deployment

```bash
sudo bash deploy-ariane.sh <git-repository-url>
sudo bash setup-tls.sh ariane.example.org admin@example.org
sudo bash setup-backup-cron.sh 03:00
sudo bash security-check.sh
```

`deploy-ariane.sh` binds Uvicorn to `127.0.0.1`, configures Nginx, enables UFW, creates a protected admin token, and installs a hardened systemd service. TLS requires a DNS name that resolves to the server.

## Service operations

```bash
bash check-ariane.sh
sudo bash restart-ariane.sh
bash logs-analyzer.sh 100
bash health-monitor.sh 30
```

The application service is `ariane`. Application logs are stored in the systemd journal. Nginx logs are stored under `/var/log/nginx`.

## Backups

The backup installer creates one systemd timer. It removes the old cron entry to avoid duplicate runs. Backups are locked with `flock`, verified after creation, and accompanied by SHA-256 checksum files.

```bash
sudo bash setup-backup-cron.sh 03:00
sudo systemctl start ariane-backup.service
systemctl list-timers ariane-backup.timer
sudo /usr/local/sbin/ariane-restore --list
sudo /usr/local/sbin/ariane-restore --data /backup/ariane-full-YYYYMMDD_HHMMSS.tar.gz
```

The local `/backup` directory is root-only. Copy backups and checksum files to encrypted storage outside the VM. Test restore procedures regularly.

## Configuration

Service secrets are stored in `/etc/ariane/ariane.env` with restricted permissions. Backup settings are stored in `/etc/ariane/backup.env`.

The `/api/clear-cache` endpoint requires the `X-ARIANE-Admin-Token` header and is blocked by the public Nginx site. Use it only through a trusted local administrative channel.

## Updates

Deployment does not run a full operating system upgrade. Apply OS updates during a planned maintenance window. Use a pinned release or commit for production deployments and keep a tested rollback copy.
