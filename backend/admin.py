from __future__ import annotations

import csv
import gzip
import html
import io
import json
import logging
import os
import secrets
from collections import Counter
from datetime import datetime, time, timedelta, timezone
from pathlib import Path
from urllib.parse import urlencode

from fastapi import APIRouter, Depends, HTTPException, Query, Request, status
from fastapi.responses import HTMLResponse, Response
from fastapi.security import HTTPBasic, HTTPBasicCredentials


router = APIRouter()
security = HTTPBasic()
logger = logging.getLogger("ariane.audit")
AUDIT_LOG_PATH = Path(os.getenv("ARIANE_AUDIT_LOG", "/var/log/ariane/audit.jsonl"))
STATUS_PATH = Path(os.getenv("ARIANE_ADMIN_STATUS", "/run/ariane/admin-status.json"))
REQUEST_EVENTS = {
    "classification_completed",
    "batch_item_completed",
    "manual_evidence_completed",
}
ERROR_MARKERS = ("error", "exception")
RETENTION_DAYS = 30


def _context(request: Request) -> dict:
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "log_type": "ariane_audit",
        "request_id": getattr(request.state, "request_id", ""),
        "source_ip": request.client.host if request.client else "unknown",
        "method": request.method,
        "path": request.url.path,
        "user_agent": request.headers.get("user-agent", "")[:300],
    }


def _log_admin(request: Request, event: str, **fields) -> None:
    record = {**_context(request), "event": event, **fields}
    logger.info(json.dumps(record, ensure_ascii=True, separators=(",", ":")))


def _require_admin(
    request: Request,
    credentials: HTTPBasicCredentials = Depends(security),
) -> str:
    expected_user = os.getenv("ARIANE_ADMIN_USER", "")
    expected_password = os.getenv("ARIANE_ADMIN_PASSWORD", "")
    user_matches = secrets.compare_digest(
        credentials.username.encode("utf-8"), expected_user.encode("utf-8")
    )
    password_matches = secrets.compare_digest(
        credentials.password.encode("utf-8"), expected_password.encode("utf-8")
    )
    if not expected_user or not expected_password or not user_matches or not password_matches:
        _log_admin(request, "admin_login_failed", attempted_user=credentials.username[:100])
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid administrator credentials",
            headers={"WWW-Authenticate": 'Basic realm="ARIANE audit"'},
        )
    _log_admin(request, "admin_login_succeeded", admin_user=expected_user)
    return credentials.username


def _iter_log_paths():
    parent = AUDIT_LOG_PATH.parent
    name = AUDIT_LOG_PATH.name
    paths = [AUDIT_LOG_PATH, parent / f"{name}.1"]
    paths.extend(sorted(parent.glob(f"{name}.*.gz"), reverse=True))
    return [path for path in paths if path.is_file()]


def _load_records(max_records: int = 100000) -> list[dict]:
    records: list[dict] = []
    for path in _iter_log_paths():
        opener = gzip.open if path.suffix == ".gz" else open
        try:
            with opener(path, "rt", encoding="utf-8", errors="replace") as stream:
                for line in stream:
                    try:
                        record = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if record.get("log_type") == "ariane_audit":
                        records.append(record)
                        if len(records) >= max_records:
                            break
        except OSError:
            continue
        if len(records) >= max_records:
            break
    records.sort(key=lambda item: item.get("timestamp", ""), reverse=True)
    return records


def _parse_timestamp(value: str) -> datetime | None:
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)
    except (TypeError, ValueError):
        return None


def _time_bounds(period: str, date_from: str, date_to: str) -> tuple[datetime, datetime]:
    now = datetime.now(timezone.utc)
    if period == "today":
        return datetime.combine(now.date(), time.min, tzinfo=timezone.utc), now
    if period == "7d":
        return now - timedelta(days=7), now
    if period == "30d":
        return now - timedelta(days=30), now
    try:
        start = datetime.combine(datetime.strptime(date_from, "%Y-%m-%d").date(), time.min, tzinfo=timezone.utc)
        end = datetime.combine(datetime.strptime(date_to, "%Y-%m-%d").date(), time.max, tzinfo=timezone.utc)
        return start, min(end, now)
    except ValueError:
        return now - timedelta(days=7), now


def _matches(record: dict, *, start: datetime, end: datetime, query: str, gene: str, event: str) -> bool:
    timestamp = _parse_timestamp(record.get("timestamp", ""))
    if not timestamp or timestamp < start or timestamp > end:
        return False
    if event and record.get("event") != event:
        return False
    input_data = record.get("input") or {}
    if gene and input_data.get("gene") != gene:
        return False
    if query:
        searchable = json.dumps(
            {
                "request_id": record.get("request_id"),
                "source_ip": record.get("source_ip"),
                "input": input_data,
                "result": record.get("result"),
                "error": record.get("error"),
            },
            ensure_ascii=False,
        ).lower()
        if query.lower() not in searchable:
            return False
    return True


def _filtered_records(records: list[dict], params: dict) -> tuple[list[dict], datetime, datetime]:
    start, end = _time_bounds(params["period"], params["date_from"], params["date_to"])
    filtered = [
        record
        for record in records
        if _matches(
            record,
            start=start,
            end=end,
            query=params["q"],
            gene=params["gene"],
            event=params["event"],
        )
    ]
    return filtered, start, end


def _status_data() -> dict:
    try:
        return json.loads(STATUS_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"status": "unavailable"}


def _safe_json(value) -> str:
    return html.escape(json.dumps(value or {}, ensure_ascii=False, indent=2))


def _url(params: dict, **updates) -> str:
    values = {**params, **updates}
    return "?" + urlencode({key: value for key, value in values.items() if value not in ("", None)})


def _export_response(records: list[dict], export_format: str) -> Response:
    if export_format == "json":
        return Response(
            json.dumps(records, ensure_ascii=False, indent=2),
            media_type="application/json",
            headers={"Content-Disposition": "attachment; filename=ariane-audit.json"},
        )
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["timestamp", "request_id", "source_ip", "event", "input", "result", "error"])
    for record in records:
        writer.writerow([
            record.get("timestamp", ""),
            record.get("request_id", ""),
            record.get("source_ip", ""),
            record.get("event", ""),
            json.dumps(record.get("input") or {}, ensure_ascii=False),
            json.dumps(record.get("result") or {}, ensure_ascii=False),
            record.get("error", ""),
        ])
    return Response(
        output.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=ariane-audit.csv"},
    )


@router.get("/admin/audit", response_class=HTMLResponse, include_in_schema=False)
async def admin_audit(
    request: Request,
    period: str = "7d",
    date_from: str = "",
    date_to: str = "",
    q: str = Query(default="", max_length=200),
    gene: str = "",
    event: str = "",
    page: int = 1,
    page_size: int = 50,
    detail: str = Query(default="", max_length=100),
    export: str = "",
    admin_user: str = Depends(_require_admin),
):
    period = period if period in {"today", "7d", "30d", "custom"} else "7d"
    gene = gene if gene in {"", "BRCA1", "BRCA2"} else ""
    page = max(1, page)
    page_size = page_size if page_size in {25, 50, 100} else 50
    params = {
        "period": period,
        "date_from": date_from,
        "date_to": date_to,
        "q": q,
        "gene": gene,
        "event": event,
        "page": page,
        "page_size": page_size,
    }
    all_records = _load_records()
    filtered, start, end = _filtered_records(all_records, params)
    if export in {"csv", "json"}:
        _log_admin(request, "admin_audit_exported", admin_user=admin_user, format=export, count=len(filtered))
        return _export_response(filtered, export)

    request_records = [record for record in filtered if record.get("event") in REQUEST_EVENTS]
    error_records = [record for record in filtered if any(value in record.get("event", "") for value in ERROR_MARKERS)]
    durations = {
        record.get("request_id"): float(record.get("duration_ms", 0))
        for record in filtered
        if record.get("event") == "request_completed" and record.get("duration_ms") is not None
    }
    request_durations = [durations.get(record.get("request_id"), 0) for record in request_records]
    class_counts = Counter(
        str((record.get("result") or {}).get("predicted_class"))
        for record in request_records
        if (record.get("result") or {}).get("predicted_class") is not None
    )
    variants = Counter(
        f"{(record.get('input') or {}).get('gene', '')} {(record.get('input') or {}).get('c_notation', '')}".strip()
        for record in request_records
    )
    variants.pop("", None)
    unique_ips = len({record.get("source_ip") for record in request_records if record.get("source_ip")})
    average_duration = sum(request_durations) / len(request_durations) if request_durations else 0
    maximum_duration = max(request_durations, default=0)
    slow_records = sorted(
        request_records,
        key=lambda record: durations.get(record.get("request_id"), 0),
        reverse=True,
    )[:10]
    login_success = sum(record.get("event") == "admin_login_succeeded" for record in filtered)
    login_failed = sum(record.get("event") == "admin_login_failed" for record in filtered)

    total_pages = max(1, (len(filtered) + page_size - 1) // page_size)
    page = min(page, total_pages)
    page_records = filtered[(page - 1) * page_size:page * page_size]
    event_options = sorted({str(record.get("event", "")) for record in all_records if record.get("event")})

    detail_records = [record for record in all_records if detail and record.get("request_id") == detail]
    detail_html = ""
    if detail:
        detail_html = "<section><h2>Request detail</h2>" + "".join(
            f"<details open><summary>{html.escape(str(record.get('timestamp', '')))} | {html.escape(str(record.get('event', '')))}</summary><pre>{_safe_json(record)}</pre></details>"
            for record in detail_records
        ) + "</section>"

    rows = []
    for record in page_records:
        request_id = str(record.get("request_id", ""))
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(record.get('timestamp', '')))}</td>"
            f"<td><a href='{html.escape(_url(params, detail=request_id))}'>{html.escape(request_id[:12])}</a></td>"
            f"<td>{html.escape(str(record.get('source_ip', '')))}</td>"
            f"<td>{html.escape(str(record.get('event', '')))}</td>"
            f"<td><details><summary>Input</summary><pre>{_safe_json(record.get('input'))}</pre></details></td>"
            f"<td><details><summary>Result</summary><pre>{_safe_json(record.get('result'))}</pre></details></td>"
            f"<td>{html.escape(str(record.get('error', '')))}</td>"
            "</tr>"
        )

    status_data = _status_data()
    status_cards = "".join(
        f"<div class='card'><span>{html.escape(label)}</span><strong>{html.escape(str(status_data.get(key, 'unknown')))}</strong></div>"
        for key, label in [
            ("ariane", "ARIANE"), ("nginx", "Nginx"), ("certificate", "Certificate"),
            ("disk", "Disk"), ("last_backup", "Last backup"),
        ]
    )
    class_cards = "".join(
        f"<div class='card'><span>Class {number}</span><strong>{class_counts.get(str(number), 0)}</strong></div>"
        for number in range(1, 6)
    )
    top_variants = "".join(
        f"<li>{html.escape(variant)} <strong>{count}</strong></li>" for variant, count in variants.most_common(10)
    ) or "<li>No data</li>"
    slow_rows = "".join(
        f"<tr><td>{html.escape(str(record.get('request_id', '')))}</td><td>{html.escape(str((record.get('input') or {}).get('gene', '')))} {html.escape(str((record.get('input') or {}).get('c_notation', '')))}</td><td>{durations.get(record.get('request_id'), 0):.1f} ms</td></tr>"
        for record in slow_records
    )
    event_select = "<option value=''>All events</option>" + "".join(
        f"<option value='{html.escape(option)}'{' selected' if event == option else ''}>{html.escape(option)}</option>"
        for option in event_options
    )
    prev_link = _url(params, page=max(1, page - 1))
    next_link = _url(params, page=min(total_pages, page + 1))
    export_csv = _url(params, export="csv", page=None)
    export_json = _url(params, export="json", page=None)

    _log_admin(request, "admin_audit_viewed", admin_user=admin_user, filters=params, count=len(filtered))
    body = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="robots" content="noindex,nofollow"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ARIANE audit</title><style>
body{{font-family:system-ui,sans-serif;margin:20px;background:#f4f6f8;color:#17202a}} h1,h2{{margin-bottom:10px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0}}
.card{{background:white;border:1px solid #d7dde3;border-radius:7px;padding:12px}} .card span{{display:block;color:#5c6773;font-size:12px}} .card strong{{font-size:20px}}
form{{background:white;border:1px solid #d7dde3;padding:12px;display:flex;gap:8px;flex-wrap:wrap;align-items:end}} label{{font-size:12px}} input,select,button{{display:block;padding:7px}}
table{{border-collapse:collapse;width:100%;background:white;margin-top:10px}} th,td{{border:1px solid #d7dde3;padding:7px;text-align:left;vertical-align:top;font-size:12px}} th{{background:#e8edf2;position:sticky;top:0}} pre{{white-space:pre-wrap;max-width:520px;max-height:360px;overflow:auto}}
.columns{{display:grid;grid-template-columns:1fr 2fr;gap:14px}} section{{margin-top:20px}} nav{{display:flex;gap:12px;margin:12px 0}} @media(max-width:850px){{.columns{{grid-template-columns:1fr}}}}
</style></head><body><h1>ARIANE audit dashboard</h1>
<p>UTC range: {start.isoformat()} to {end.isoformat()}. Retention: {RETENTION_DAYS} days.</p>
<h2>System status</h2><div class="grid">{status_cards}</div>
<form method="get"><label>Period<select name="period"><option value="today"{' selected' if period == 'today' else ''}>Today</option><option value="7d"{' selected' if period == '7d' else ''}>7 days</option><option value="30d"{' selected' if period == '30d' else ''}>30 days</option><option value="custom"{' selected' if period == 'custom' else ''}>Custom</option></select></label>
<label>From<input type="date" name="date_from" value="{html.escape(date_from)}"></label><label>To<input type="date" name="date_to" value="{html.escape(date_to)}"></label>
<label>Search<input name="q" value="{html.escape(q)}" placeholder="variant, gene, IP, request ID"></label>
<label>Gene<select name="gene"><option value="">All genes</option><option{' selected' if gene == 'BRCA1' else ''}>BRCA1</option><option{' selected' if gene == 'BRCA2' else ''}>BRCA2</option></select></label>
<label>Event<select name="event">{event_select}</select></label><label>Rows<select name="page_size"><option>25</option><option{' selected' if page_size == 50 else ''}>50</option><option{' selected' if page_size == 100 else ''}>100</option></select></label><button type="submit">Apply</button></form>
<div class="grid"><div class="card"><span>Requests</span><strong>{len(request_records)}</strong></div><div class="card"><span>Errors</span><strong>{len(error_records)}</strong></div><div class="card"><span>Unique IPs</span><strong>{unique_ips}</strong></div><div class="card"><span>Average duration</span><strong>{average_duration:.1f} ms</strong></div><div class="card"><span>Maximum duration</span><strong>{maximum_duration:.1f} ms</strong></div><div class="card"><span>Authenticated admin access</span><strong>{login_success}</strong></div><div class="card"><span>Failed logins</span><strong>{login_failed}</strong></div></div>
<h2>Classification results</h2><div class="grid">{class_cards}</div>
<div class="columns"><section><h2>Top variants</h2><ol>{top_variants}</ol></section><section><h2>Slow requests</h2><table><tr><th>Request ID</th><th>Variant</th><th>Duration</th></tr>{slow_rows}</table></section></div>
{detail_html}<section><h2>Events</h2><nav><a href="{html.escape(export_csv)}">Export CSV</a><a href="{html.escape(export_json)}">Export JSON</a></nav>
<table><thead><tr><th>Time UTC</th><th>Request ID</th><th>Source IP</th><th>Event</th><th>Input</th><th>Result</th><th>Error</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<nav><a href="{html.escape(prev_link)}">Previous</a><span>Page {page} of {total_pages}</span><a href="{html.escape(next_link)}">Next</a></nav></section></body></html>"""
    return HTMLResponse(
        body,
        headers={
            "Cache-Control": "no-store, max-age=0",
            "Pragma": "no-cache",
            "X-Robots-Tag": "noindex, nofollow",
            "Content-Security-Policy": "default-src 'none'; style-src 'unsafe-inline'; form-action 'self'; base-uri 'none'",
        },
    )
