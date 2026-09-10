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

Variant types whose automatic path does not require SpliceAI are not sent to
the SpliceAI service. Their classification remains complete without a score,
and the audit records `required_for_classification: false`.

When a protein PS1 reference candidate exists, the SpliceAI comparison of the
assessed and reference variants is part of the required evidence. A temporary
failure for either variant makes the classification retryable and incomplete.

## Request and traffic limits

The reference deployment currently permits 5 HTTP requests per second from one
client IP with a temporary burst of 20 requests. This is a safety ceiling, not a
recommended sustained rate. Clients should send batches sequentially and wait
for each response.

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

The unversioned `/api/classify` and `/api/classify/batch` routes remain available
for the current web interface and compatibility. External integrations should
use `/api/v1`. Only the versioned route has the documented public response
contract.

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
