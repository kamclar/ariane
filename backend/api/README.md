# API transport

This package maps HTTP requests and authentication context to application
services. It may set response headers, enforce request quotas, and translate
application failures into HTTP errors. It must not evaluate evidence or derive
classification criteria.

- `classification.py` handles single and batch classification transport.
- `manual.py` handles manual evidence, PS1 reference resolution and normalization.
- `system.py` exposes health, rules, resources and cache maintenance.
- `public.py` defines the stable versioned public API projection.
- `routing.py` registers equivalent browser and public routes without duplicating
  authentication declarations.
- `auth.py` and `session.py` keep API-key and browser authentication separate.
- `admin.py` and `review.py` expose restricted operational workflows.
- `middleware.py`, `errors.py` and `audit.py` provide application-wide HTTP
  behavior.
- `frontend.py` serves the HTML shell and static assets.

Equivalent browser and public endpoints are registered with
`PairedApiRoutes`. A declaration such as `@paired.get("/rules")` creates both
`/ui-api/rules` and `/api/rules`. The helper always applies the signed browser
session to the UI route, the API key to the public route, and excludes only the
UI route from OpenAPI. Route modules do not repeat these security decorators.
