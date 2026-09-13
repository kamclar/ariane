# Application contracts

These Pydantic models cross the transport, application service and persistence
boundaries. They are grouped by workflow so lower layers do not need to import
the HTTP package.

Contracts contain validation and serialization rules. They do not perform
evidence lookups or assign criteria.
