# Infrastructure

This package contains mutable runtime storage, cache locations, health state,
lookup execution helpers and repository adapters.

`DataHealthRegistry` is created by `backend/bootstrap.py` and injected into the
runtime components that report degraded sources. There is no process-global
health dictionary.
