# ARIANE Server Operations Scripts

Sbírka skriptů pro provozování a monitoring ARIANE aplikace na produkčním serveru.

## 📋 Skripty

### 1. **check-ariane.sh** - Status Check
Rychlý přehled stavu aplikace.
```bash
bash check-ariane.sh
```
Zobrazí: Service status, health check, CPU/RAM, disk space, Nginx status, poslední chyby.

### 2. **backup-ariane.sh** - Backup Dat
Zálohuje `backend/data` a celý aplikaci adresář bez `venv`, `.git` a běhových souborů.
```bash
bash backup-ariane.sh
```
Vytvoří timestampovaný `.tar.gz` soubor v `/backup/` adresáři.

### 3. **restore-ariane.sh** - Restore Backup
Obnoví zálohu ARIANE z archivu.
```bash
bash restore-ariane.sh --list
bash restore-ariane.sh /backup/ariane-full-20260710_030000.tar.gz
```

### 4. **restart-ariane.sh** - Restart Aplikace
Bezpečně restartuje ARIANE a Nginx.
```bash
bash restart-ariane.sh
```

### 4. **health-monitor.sh** - Continuous Monitoring
Sleduje zdraví aplikace v reálném čase (každých 30 sekund).
```bash
bash health-monitor.sh
```
Zastavení: Ctrl+C

### 5. **logs-analyzer.sh** - Log Analysis
Analyzuje poslední logy a hledá chyby.
```bash
bash logs-analyzer.sh [počet řádků=100]
```
Příklad: `bash logs-analyzer.sh 50`

### 6. **install-ariane-service.sh** - Systemd Service
Instaluje ARIANE jako systemd service (autostart).
```bash
sudo bash install-ariane-service.sh
```

### 7. **setup-monitoring-tools.sh** - Install Tools
Nainstaluje monitoring nástroje (htop, iotop, nethogs, lnav).
```bash
sudo bash setup-monitoring-tools.sh
```

### 8. **setup-backup-cron.sh** - Automatic Backups
Nastaví denní automatické backupy pomocí cron a/nebo systemd timeru.
```bash
sudo bash setup-backup-cron.sh
```

### 9. **restore-ariane.sh** - Restore Backup
Obnoví zálohu z `/backup`.
```bash
bash restore-ariane.sh --list
bash restore-ariane.sh /backup/ariane-full-20260710_030000.tar.gz
```

### 10. **deploy-ariane.sh** - Full Deployment
Kompletní nasazení aplikace (klonování, setup, spuštění).
```bash
sudo bash deploy-ariane.sh <git-repository-url>
```

### 10. **security-check.sh** - Security Audit
Ověřuje bezpečnostní nastavení.
```bash
bash security-check.sh
```

---

## 🚀 Doporučený setup

### Při prvním nasazení:
```bash
# 1. Instaluj service
sudo bash install-ariane-service.sh

# 2. Instaluj monitoring nástroje
sudo bash setup-monitoring-tools.sh

# 3. Nastav automatické backupy
sudo bash setup-backup-cron.sh

# 4. Zkontroluj status
bash check-ariane.sh
```

### Denní rutina:
```bash
# Ráno
bash check-ariane.sh

# Při potřebě restartovat
bash restart-ariane.sh

# V případě problémů
bash logs-analyzer.sh
```

---

## 📌 Instalace na VM

1. Nakopíruj adresář `scripts/server-ops/` na server:
```bash
scp -r scripts/server-ops/ ubuntu@<vm-ip>:~/
```

2. Udělej skripty spustitelné:
```bash
chmod +x ~/server-ops/*.sh
```

3. Spusť je podle potřeby (viz výše).

---

## ⚙️ Konfigurace

Některé skripty si vyžadují konfiguraci v souboru `.env.ops`:

```bash
# Vytvoř .env.ops v domovském adresáři
cat > ~/.env.ops << 'EOF'
# ARIANE Configuration
ARIANE_HOME="/home/ubuntu/ariane"
ARIANE_USER="ubuntu"
ARIANE_PORT="8000"
BACKUP_DIR="/backup"
BACKUP_RETENTION_DAYS="30"
EOF
```

---

## 📊 Logování

Logy najdeš v:
- **Aplikace**: `journalctl -u ariane -f`
- **Nginx**: `/var/log/nginx/access.log` a `/var/log/nginx/error.log`
- **System**: `/var/log/syslog`

---

## 🔧 Troubleshooting

### Aplikace nereaguje?
```bash
bash restart-ariane.sh
bash check-ariane.sh
```

### Disk space problém?
```bash
bash backup-ariane.sh
# Pak smaž staré backupy
ls -lah /backup/
```

### Problém s Nginx?
```bash
sudo nginx -t
sudo systemctl restart nginx
```

---

## 📞 Support

Při problémech se podívej na logy:
```bash
bash logs-analyzer.sh 100
```

---

**Created for ARIANE v1.8.0 | e-infra.cz OpenStack**
