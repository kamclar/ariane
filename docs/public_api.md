# ARIANE public API v1

The public API is a versioned beta interface for programmatic variant
classification. It uses the same normalization service, evidence providers and
production DAG as the ARIANE web application. It does not contain a separate
classification implementation.

Base URL:

```text
https://ariane-app.duckdns.org/api/v1
```

The current contract version is `1.0`. Compatible additions may be made within
v1. Breaking request or response changes require a new API path.

## Authentication

Classification endpoints require an individual API key in the
`X-ARIANE-API-Key` header. The capabilities endpoint is public and reports the
authentication requirement without exposing registered key IDs.

Store the key in an environment variable instead of a script or repository:

```bash
export ARIANE_API_KEY='value-provided-by-the-administrator'
```

The server stores only a SHA-256 digest of each high-entropy key. A key can be
disabled independently. Missing, invalid and disabled keys receive HTTP 401.
If the protected registry is missing or malformed, classification fails closed
with HTTP 503.

## Capabilities

Read the active application version, supported genes, transcripts, policies and
request limits before starting a run:

```bash
curl https://ariane-app.duckdns.org/api/v1/capabilities
```

The capability response is the authoritative source for current limits. The
reference deployment currently accepts a maximum of 10 variants per synchronous
v1 batch. Five variants are recommended when the results may not yet be cached.
It also reports the daily per-key classification allowance and the concurrent
classification limits enforced by the reference deployment.

## Single classification

```bash
curl -X POST https://ariane-app.duckdns.org/api/v1/classify \
  -H "Content-Type: application/json" \
  -H "X-ARIANE-API-Key: $ARIANE_API_KEY" \
  -H "X-Request-ID: local-run-0001" \
  -d '{"gene":"BRCA1","c_notation":"c.4185G>A"}'
```

The response contains three separate sections:

- `metadata`: API, application, VCEP policy and classifier identities, request
  ID, processing time and cache status
- `classification`: the automatic Module 1 result and its evidence audit
- `manual_review`: review aids that may identify useful expert follow-up

`manual_review.affects_automatic_classification` is always `false`. A review
recommendation is not an ACMG/AMP criterion, adds no points and cannot change
the automatic class. An amended working result requires separately documented
manual evidence and backend evaluation.

## Batch classification

```bash
curl -X POST https://ariane-app.duckdns.org/api/v1/classify/batch \
  -H "Content-Type: application/json" \
  -H "X-ARIANE-API-Key: $ARIANE_API_KEY" \
  -d '{"variants":[
    {"gene":"BRCA1","c_notation":"c.4185G>A"},
    {"gene":"BRCA2","c_notation":"c.7805+9T>G"}
  ]}'
```

Input order is preserved. Each item has its original zero-based `index` and a
status of `ok` or `error`. An invalid item does not prevent other variants in
the request from being classified.

Item errors contain:

- `code`: stable machine-readable category
- `message`: explanation for a person or log
- `retryable`: whether a later identical request may reasonably succeed after
  a temporary service condition changes

Clients must not retry `invalid_request`, `invalid_variant` or
`variant_not_classifiable` without correcting the input. A retryable response
still provides no classification and must never be interpreted as negative or
benign evidence.

For a missense, confirmed in-frame, synonymous or relevant intronic variant,
SpliceAI is required by the automatic Figure 1A path. If the required score is
unavailable, the item has `status: "error"` and contains no classification.
Temporary source failures use `spliceai_temporarily_unavailable` with
`retryable: true`. Missing GRCh38 coordinates use
`spliceai_coordinates_unavailable` with `retryable: false`. Invalid or
methodologically incomplete source responses use `spliceai_result_unavailable`
with `retryable: false`.

The applicable type is determined from the normalized consequence on the
configured reference transcript. A DNA-level `delins` that encodes one amino
acid substitution is therefore a missense variant for Figure 1A and requires
SpliceAI. A `delins` that produces a frameshift or a simple nonsense consequence
uses the corresponding PTC branch instead.

Variant types whose automatic path does not require SpliceAI are not sent to
the SpliceAI service. Their classification remains complete without a score,
and the audit records `required_for_classification: false`.

When a protein PS1 reference candidate exists, the SpliceAI comparison of the
assessed and reference variants is part of the required evidence. A temporary
failure for either variant makes the classification retryable and incomplete.

The same publication gate covers the other required upstream inputs. Both
gene-policy gnomAD datasets must complete. An applicable Appendix G structural
population path must reach an explicit decision. BayesDel_noAF must be
available for a missense or in-frame variant inside a functional domain when
SpliceAI is below 0.2.
Failures use `population_evidence_unavailable`,
`structural_population_evidence_unavailable`,
`bayesdel_temporarily_unavailable`, `bayesdel_result_unavailable`,
or `protein_interval_unavailable`. These responses contain no classification.

A completed negative observation does not count as a failure. Examples include
a score below a criterion threshold, an absent variant in a successfully
queried dataset, a filtered gnomAD record, and a policy-defined not-applicable
branch. Live ClinVar and ClinGen ERepo lookups remain external comparisons, so
their failure does not invalidate the automatic classification. Automatic
PVS1 RNA can use only an exact assertion in the versioned, checksum-validated
local ERepo registry. A live ERepo response is not a fallback classification
source.

An unresolved pathogenic-founder check is a declared expert-review state, not
a provider failure. The automatic result contains no BA1 or BS1 from that
frequency value and records the reason. Other independently completed evidence
may still produce the Module 1 result.

## Request and traffic limits

The reference deployment currently permits 5 HTTP requests per second from one
client IP with a temporary burst of 20 requests. This is a safety ceiling, not a
recommended sustained rate. Authenticated requests are additionally limited to
30 HTTP requests per minute for each API key with a burst of 3. Clients should
send batches sequentially and wait for each response.

Each API key may reserve at most 5,000 variant classifications per UTC day. A
single request consumes one unit and a batch consumes one unit for every submitted
variant that passes request validation, including an item that later fails evidence
retrieval. Invalid batch items do not consume classification quota, but remain
subject to the HTTP request-rate limit. The reservation is made atomically before
classification starts, so parallel requests and batches cannot exceed the
allowance. Successful responses include `X-RateLimit-Limit`,
`X-RateLimit-Remaining` and `X-RateLimit-Reset`. An exhausted allowance returns
HTTP 429 with `daily_classification_quota_exceeded` and a `Retry-After` header.

The reverse proxy permits at most four concurrent classification requests from
one IP address and at most two concurrent requests for one API key. These limits
apply to classification endpoints, not to static files or the health endpoint.
They bound active work while allowing normal browser page loading.

The proxy returns HTTP 429 when the request rate is exceeded. Respect
`Retry-After` when present. Use bounded retries only for HTTP 429, 502, 503 and
504. Do not retry HTTP 422 unchanged.

The server returns `X-Request-ID`, `X-ARIANE-API-Version` and `X-ARIANE-Version`
headers. A client may supply an `X-Request-ID` containing a short local run ID.
Record the returned request ID with each result for troubleshooting.

The reference client reads the key from `ARIANE_API_KEY`. The command-line
`--api-key` option is available, but the environment variable is preferred
because command arguments may be retained in shell history or process lists.

## Input policy

Submit reference-transcript c. HGVS or a supported genomic HGVS form. The
normalization engine derives and verifies the reference transcript and protein
consequence. Do not manufacture a protein consequence when it is unknown.

The API accepts variant information only. Do not submit patient names, dates of
birth, sample identifiers or clinical notes.

## Legacy endpoints

The unversioned `/api/classify` and `/api/classify/batch` routes require the
same API key as v1, so they cannot bypass authentication. They remain only for
compatibility and do not have the documented v1 response contract. External
integrations must use `/api/v1`.

The web application does not contain an API key. It uses internal `/ui-api`
routes with a short-lived, server-signed HttpOnly session created when the main
page loads. Mutating requests must also carry the browser's matching Origin
header. These routes are intended only for the same-origin browser and are
limited separately by client IP. Obtaining a browser session is not user
authentication, so request limits remain necessary.

There is deliberately no quota keyed only by the short-lived UI session. The
main page can issue a new anonymous session, so such a quota would be easy to
reset and would not protect the service. Interactive classification is limited
by IP address and by the concurrent classification ceiling instead.

For local programmatic testing, create a development key in the default ignored
runtime-data directory before starting ARIANE:

```powershell
python scripts/manage_api_keys.py create --id local-development
$env:ARIANE_API_KEYS_FILE = ".runtime-data/api_keys.json"
```

The command prints the plaintext key once. The registry contains only its hash.
Local UI use does not require an API key, but it still requires a configured
`ARIANE_UI_SESSION_SECRET`.

## Python example

The repository contains `scripts/ariane_api_client.py`. It reads a tab-separated
file with `gene` and `c_notation` columns, checks server capabilities, submits
small sequential batches and writes one JSON object per input row. The added
`input_index` is the zero-based row position across the complete input file.
The client repeats only batch items carrying `retryable: true`, with a bounded
number of attempts. Other errors are written directly to the output.

```bash
python scripts/ariane_api_client.py variants.tsv --output results.jsonl
```
