"""Command-line entry point for the read-only ARIANE source watch."""
from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.source_watch import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
