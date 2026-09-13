# ARIANE server operations

Production scripts for an Ubuntu server with systemd and Nginx.

## Initial deployment

```bash
sudo bash deploy-ariane.sh <git-repository-url>
sudo bash setup-tls.sh ariane-app.duckdns.org admin@example.org
sudo bash setup-admin.sh admin
sudo bash setup-backup-cron.sh 03:00
sudo bash security-check.sh
```

`deploy-ariane.sh` uses `ariane-app.duckdns.org` by default. Set `ARIANE_SERVER_NAME` to override it. The DNS A record must resolve to `195.113.167.59` before running the TLS setup.

The deploy script binds Uvicorn to `127.0.0.1`, configures Nginx, enables UFW, creates a protected admin token, and installs a hardened systemd service.

Nginx serves the application only for `ariane-app.duckdns.org`. Direct IP access and unknown hostnames are rejected on HTTP and HTTPS.

## Service operations

```bash
bash check-ariane.sh
sudo bash restart-ariane.sh
bash logs-analyzer.sh 100
bash health-monitor.sh 30
```

The application service is `ariane`. Application logs are stored in the systemd journal. Nginx logs are stored under `/var/log/nginx`.

## Local SpliceAI service

ARIANE uses a separate SpliceAI container bound to `127.0.0.1:8082`. Install it
on an existing server after updating the repository:

Port 8081 is used by the ARIANE issue tracker on the current host. The SpliceAI
port can be changed through `SPLICEAI_PORT`. The installer rejects an occupied
port before replacing the systemd service.

```bash
sudo bash /home/ubuntu/ariane/scripts/server-ops/install-spliceai-service.sh
```

The script reads the immutable image digest from the active SpliceAI profile. It
tests the image on a temporary private port before changing `/etc/ariane/ariane.env`.
The test requires the Appendix J parameters, the BRCA1 and BRCA2 reference
transcripts and the exact delta, REF and ALT scores in the versioned validation cases. A failed
validation leaves ARIANE on its previous configuration.

The production container runs three model workers. Installation validates them
concurrently so that every worker loads and traces the model before user traffic
is accepted. The first request handled by an uninitialized worker can otherwise
take much longer than subsequent requests.

Useful checks:

```bash
systemctl status ariane-spliceai
journalctl -u ariane-spliceai -f
python3 /home/ubuntu/ariane/scripts/validate_spliceai_service.py \
  --url http://127.0.0.1:8082/spliceai/
```

The service has no public listener and no server-side database. ARIANE keeps its
own profile-specific runtime cache. Do not replace the digest with a mutable
`latest` tag. A new image is promoted by updating the profile and validation
case together, running the test suite and then running the installer.

`restart-ariane.sh` verifies the running local service against the active profile
before restarting the application. If a repository update contains a new
SpliceAI profile, run `install-spliceai-service.sh`. It installs and verifies the
new image before restarting ARIANE.

ARIANE itself runs as one application worker. Its in-memory request gate and the
JSON SpliceAI runtime cache are process-local, so multiple application workers
would bypass the configured concurrency bound and could overwrite concurrent
cache updates. `restart-ariane.sh` rejects an older service definition with more
workers or missing persistent runtime directories and points to
`install-ariane-service.sh` for the one-time service migration.

The Broad SpliceAI Lookup wrapper is MIT licensed. The pinned SpliceAI commit is
GPLv3 and its model weights are CC BY-NC 4.0. The current deployment is intended
for free academic development and evaluation. Confirm the licence scope with
Illumina before allowing use as part of a paid diagnostic service.

## Public API keys

Configure the protected registry on an existing server, then create one key
for each external integration:

```bash
sudo bash /home/ubuntu/ariane/scripts/server-ops/enable-public-api-auth.sh
sudo /home/ubuntu/ariane/venv/bin/python \
  /home/ubuntu/ariane/scripts/manage_api_keys.py \
  --file /etc/ariane/api-keys.json create \
  --id external-batch-01 \
  --description "External batch testing"
sudo bash /home/ubuntu/ariane/scripts/server-ops/restart-ariane.sh
```

The create command prints the key once. Send it through a suitable private
channel and do not add it to Git, issue reports or application logs. The server
stores only its SHA-256 digest.

List IDs or disable a key without changing other clients:

```bash
sudo /home/ubuntu/ariane/venv/bin/python \
  /home/ubuntu/ariane/scripts/manage_api_keys.py \
  --file /etc/ariane/api-keys.json list
sudo /home/ubuntu/ariane/venv/bin/python \
  /home/ubuntu/ariane/scripts/manage_api_keys.py \
  --file /etc/ariane/api-keys.json disable --id external-batch-01
```

Registry updates are read on each authenticated request. Disabling a key does
not require a service restart. The key ID, never the secret, identifies API
usage in classification statistics.

The same key is required by the legacy `/api/classify` and
`/api/classify/batch` endpoints. The browser uses `/ui-api` routes with a signed
HttpOnly session and never receives a reusable API key. Interactive
classification is limited to 30 requests per minute per IP with a burst of 3.
Authenticated routes are limited to 30 HTTP requests per minute for each key
with a burst of 3. The browser batch tool and reference API client space or
sequence requests accordingly.

The backend also reserves one quota unit per submitted variant before starting
classification. The default allowance is 5,000 classifications per API key and
UTC day. Configure it with `ARIANE_API_DAILY_CLASSIFICATION_LIMIT`. Batch items
count separately, so batching cannot bypass the allowance. If quota storage is
unavailable, public classification fails closed with HTTP 503 while the UI
remains available.

Nginx permits four concurrent classification requests per IP address and two
per API key. The connection limits apply only to classification paths. Apply the
current proxy settings to an existing installation with:

```bash
sudo bash /home/ubuntu/ariane/scripts/server-ops/update-nginx-api-settings.sh
```

The general per-IP request limit applies only to `/api/*` routes. The page,
stylesheets, JavaScript modules and images are not charged against that limit,
because browsers fetch these files concurrently. Interactive classification is
limited separately on `/ui-api/classify`.

`ARIANE_UI_SESSION_SECRET` signs browser sessions. Installation scripts create
a random 32-byte secret in `/etc/ariane/ariane.env`. Keep the value outside Git
and use the same value for every ARIANE worker. ARIANE does not start if this
value is missing or shorter than 32 bytes. The restart script provisions the
value on installations created before this setting was introduced.
An ordinary restart keeps the existing value. Replacing it invalidates active
browser cookies, so users must reload the page to obtain a new session.

Structured audit events include the request ID, source IP, endpoint, submitted values, predicted class, class label, total points, and error details. Tokens and request headers are not logged.

```bash
sudo bash audit-log.sh all today
sudo bash audit-log.sh errors "7 days ago"
sudo bash audit-log.sh requests today
journalctl -u ariane --since today -o cat | grep '"log_type":"ariane_audit"'
```

Audit events can contain variant data, free-text notes, assessor names, and IP addresses. Limit journal access to administrators and define a retention period that matches local privacy requirements.

## Updating Python dependencies on an existing server

`restart-ariane.sh` synchronizes the pinned packages in `requirements.txt`
with the existing virtual environment and performs an import preflight before
restarting the service. If synchronization fails, the running service is not
stopped.

The equivalent manual recovery command is:

```bash
sudo apt-get update
sudo apt-get install -y libpq-dev python3-dev build-essential
sudo -u ubuntu /home/ubuntu/ariane/venv/bin/python -m pip install -r /home/ubuntu/ariane/requirements.txt
sudo bash /home/ubuntu/ariane/scripts/server-ops/restart-ariane.sh
```

The service intentionally refuses to start if the pinned HGVS packages or the
checksum-verified panel reference bundle are missing or inconsistent.

## Audit administration page

Configure a generated password and open the page over HTTPS:

```bash
sudo bash setup-admin.sh admin
```

The script prints a new password once and keeps it on later runs. Store it in a password manager. Run `sudo bash setup-admin.sh admin --rotate` over SSH to generate a replacement. The page is available at `https://ariane-app.duckdns.org/admin/audit` and uses browser Basic Auth. Audit files are stored in `/var/log/ariane/audit.jsonl`, readable only by the service account and administrators, and rotated daily for 30 days.

The dashboard is read-only. It provides time, text, gene, and event filters; persistent search and cache statistics; class counts; common variants; account or browser counts; performance data; request details; pagination; CSV and JSON export; login history; and service, certificate, disk, and backup status. It does not provide restart, deploy, restore, cache, or password actions.

## Backups

The backup installer creates one systemd timer. It removes the old cron entry to avoid duplicate runs. Backups are locked with `flock`, verified after creation, and accompanied by SHA-256 checksum files.

```bash
sudo bash setup-backup-cron.sh 03:00
sudo systemctl start ariane-backup.service
systemctl list-timers ariane-backup.timer
sudo /usr/local/sbin/ariane-restore --list
sudo /usr/local/sbin/ariane-restore --data /backup/ariane-full-YYYYMMDD_HHMMSS.tar.gz
sudo /usr/local/sbin/ariane-review-restore /backup/ariane-reviews-YYYYMMDD_HHMMSS.sqlite3.gz
sudo /usr/local/sbin/ariane-usage-restore /backup/ariane-usage-YYYYMMDD_HHMMSS.sqlite3.gz
```

The local `/backup` directory is root-only. Off-site backup is optional and is not configured by these scripts. Local backups do not protect against loss of the VM or its disk. Test restore procedures regularly.

Manual-review drafts and approvals are stored separately from replaceable API
caches in `/var/lib/ariane/runtime-data/review_records.sqlite3`. The backup job
uses the SQLite online-backup API, verifies database integrity, compresses the
copy and writes a SHA-256 checksum. Review-record backups follow the same
retention period as the application backups.

Classification usage events are stored in
`/var/lib/ariane/runtime-data/classification_usage.sqlite3` and are included in
the same backup job. Replaceable classification results are stored in
`/var/lib/ariane/runtime-cache/classification_results.sqlite3` and are not
backed up. The result cache has no time-based expiration. Its implementation and data
fingerprint invalidates it after a relevant update. An optional maximum age can
be configured with `ARIANE_CLASSIFICATION_CACHE_MAX_AGE_SECONDS`. The default
usage retention is 365 days.

## Configuration

Service secrets are stored in `/etc/ariane/ariane.env` with restricted permissions. Backup settings are stored in `/etc/ariane/backup.env`.

Public API key digests are stored in `/etc/ariane/api-keys.json`, owned by root
and readable by the ARIANE service group. Plaintext keys cannot be recovered
from this file.

`ARIANE_RUNTIME_DATA_DIR` must point to persistent storage writable only by the
ARIANE service account. It must not point to the runtime-cache directory.

Without user authentication, usage records identify only a persistent browser
visitor. To store authenticated account names, the reverse proxy must remove
any client-supplied identity header, set its own header after authentication,
and `ARIANE_TRUSTED_USER_HEADER` must name that header. Do not enable this
setting for a header that can be supplied directly by an unauthenticated
client.

The `/api/clear-cache` endpoint requires the `X-ARIANE-Admin-Token` header and is blocked by the public Nginx site. Use it only through a trusted local administrative channel.

## Updates

Deployment does not run a full operating system upgrade. Apply OS updates during a planned maintenance window. Use a pinned release or commit for production deployments and keep a tested rollback copy.
