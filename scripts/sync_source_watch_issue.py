"""Create, update or close the GitHub issue for source-watch review notices."""

from __future__ import annotations

import argparse
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


ISSUE_TITLE = "ARIANE source update available for review"
MARKER_PREFIX = "<!-- ariane-source-watch-fingerprint:"


def build_issue_body(
    report: dict[str, Any],
    markdown: str,
    *,
    repository: str,
    run_url: str,
) -> str:
    fingerprint = str(report.get("review_fingerprint") or "")
    prepare_url = (
        f"https://github.com/{repository}/actions/workflows/prepare-source-update.yml"
    )
    return "\n".join(
        [
            markdown.rstrip(),
            "",
            "## GitHub actions",
            "",
            f"- [Source-watch run]({run_url})",
            f"- [Prepare a draft update pull request]({prepare_url})",
            "",
            "The preparation workflow is manual. It does not deploy the update. "
            "The active ARIANE data change only after the resulting pull request "
            "is reviewed, approved and merged.",
            "",
            f"{MARKER_PREFIX}{fingerprint} -->",
            "",
        ]
    )


class GitHubIssues:
    def __init__(self, *, repository: str, token: str) -> None:
        self.repository = repository
        self.token = token
        self.base = f"https://api.github.com/repos/{repository}"

    def request(
        self,
        method: str,
        path: str,
        payload: dict[str, Any] | None = None,
    ) -> Any:
        data = None if payload is None else json.dumps(payload).encode("utf-8")
        request = urllib.request.Request(
            f"{self.base}{path}",
            data=data,
            method=method,
            headers={
                "Accept": "application/vnd.github+json",
                "Authorization": f"Bearer {self.token}",
                "X-GitHub-Api-Version": "2022-11-28",
                "User-Agent": "ARIANE-source-watch/1.0",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                if response.status == 204:
                    return None
                return json.load(response)
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(
                f"GitHub API {method} {path} returned {exc.code}: {detail}"
            ) from exc

    def find_open_notice(self) -> dict[str, Any] | None:
        query = urllib.parse.urlencode({"state": "open", "per_page": 100})
        issues = self.request("GET", f"/issues?{query}")
        return next(
            (
                issue
                for issue in issues
                if issue.get("title") == ISSUE_TITLE and "pull_request" not in issue
            ),
            None,
        )

    def create(self, body: str, *, assignee: str = "") -> dict[str, Any]:
        payload: dict[str, Any] = {"title": ISSUE_TITLE, "body": body}
        if assignee:
            payload["assignees"] = [assignee]
        return self.request("POST", "/issues", payload)

    def edit(self, number: int, **changes: Any) -> dict[str, Any]:
        return self.request("PATCH", f"/issues/{number}", changes)

    def comment(self, number: int, body: str) -> dict[str, Any]:
        return self.request("POST", f"/issues/{number}/comments", {"body": body})


def sync(
    report: dict[str, Any],
    markdown: str,
    *,
    repository: str,
    run_url: str,
    token: str,
    assignee: str = "",
) -> str:
    github = GitHubIssues(repository=repository, token=token)
    existing = github.find_open_notice()
    status = report.get("overall_status")
    if status == "review_required":
        body = build_issue_body(
            report,
            markdown,
            repository=repository,
            run_url=run_url,
        )
        marker = f"{MARKER_PREFIX}{report.get('review_fingerprint', '')} -->"
        if existing is None:
            created = github.create(body, assignee=assignee)
            return f"created issue #{created['number']}"
        if marker in str(existing.get("body") or ""):
            return f"issue #{existing['number']} already describes this change"
        number = int(existing["number"])
        github.edit(number, body=body)
        github.comment(
            number,
            "A further upstream change was detected. The issue body now contains "
            f"the latest diff and links to [workflow run]({run_url}).",
        )
        return f"updated issue #{number}"
    if status == "current" and existing is not None:
        number = int(existing["number"])
        github.comment(
            number,
            "The source watch now matches the approved repository baseline. "
            "Closing this notice.",
        )
        github.edit(number, state="closed")
        return f"closed issue #{number}"
    return "no issue change required"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--repository", required=True)
    parser.add_argument("--run-url", required=True)
    args = parser.parse_args()
    token = os.environ.get("GITHUB_TOKEN", "")
    if not token:
        raise RuntimeError("GITHUB_TOKEN is required")
    result = sync(
        json.loads(args.report.read_text(encoding="utf-8")),
        args.markdown.read_text(encoding="utf-8"),
        repository=args.repository,
        run_url=args.run_url,
        token=token,
        assignee=os.environ.get("SOURCE_WATCH_ASSIGNEE", ""),
    )
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
