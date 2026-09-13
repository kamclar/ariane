"""Import-safe access to immutable JSON reference objects."""

from __future__ import annotations

from collections.abc import Iterator, MutableMapping
from typing import Any, Callable


class LazyJsonObject(MutableMapping[str, Any]):
    """Preserve mapping-style APIs while deferring file I/O until first use."""

    def __init__(self, loader: Callable[[], dict[str, Any]]) -> None:
        self._loader = loader

    def _data(self) -> dict[str, Any]:
        return self._loader()

    def __getitem__(self, key: str) -> Any:
        return self._data()[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self._data()[key] = value

    def __delitem__(self, key: str) -> None:
        del self._data()[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data())

    def __len__(self) -> int:
        return len(self._data())

    def __repr__(self) -> str:
        return repr(self._data())
