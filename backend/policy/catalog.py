"""Lazy read-only views of gene policy owned by the policy layer."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
from typing import Any

from backend.policy.gene import domains_by_gene, transcripts_by_gene


class LazyPolicyMapping(Mapping[str, Any]):
    def __init__(self, loader: Callable[[], dict[str, Any]]) -> None:
        self._loader = loader

    def __getitem__(self, key: str) -> Any:
        return self._loader()[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._loader())

    def __len__(self) -> int:
        return len(self._loader())


FUNCTIONAL_DOMAINS = LazyPolicyMapping(domains_by_gene)
TRANSCRIPTS = LazyPolicyMapping(transcripts_by_gene)
