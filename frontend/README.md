# Frontend structure

ARIANE uses Alpine.js directly in the browser. There is no compilation or
bundling step. The server assembles the HTML templates and adds a content-based
version to every local asset URL.

## JavaScript

Each feature registers one state factory or one or more method objects on
`window.ArianeFrontend`. `app.js` passes named contributors to
`composition.js`. A state or method name may have exactly one owner. Duplicate
keys stop component creation with an error naming both contributors.

- `core.js`: application metadata and shared resources
- `classification.js`: single-variant input and classification request
- `batch.js`: batch parsing, execution and export
- `rules.js` and `graphs.js`: rule explorer and decision graphs
- `formatters.js`: presentation-only transformations
- `manual-review.js`: review state, navigation and safe prefill
- `manual-review-api.js`: evidence payloads, status and amended result request
- `manual-review-ps1.js`: protein and splice PS1 reference assistance
- `manual-review-persistence.js`: draft, approval and audit export

Feature modules must not redefine another module's state or method. They may
call a method owned by another feature after the complete component has been
composed.

## Templates

`index.html` is the page shell. Top-level templates live in `templates/`.
Manual review is further divided by evidence family and result lifecycle. The
explicit include list in `backend/presentation/frontend.py` determines assembly
order and rejects missing, duplicate or unknown include markers.

## Styles

CSS files own a feature or a reusable concern:

- `base.css`, `forms.css`, `layout.css`: shared tokens, controls and page layout
- `results.css`, `evidence.css`: classification and evidence presentation
- `manual-review.css`: manual evidence workflow
- `navigation.css`: primary mode navigation
- `batch.css`: batch mode
- `decision-paths.css`: embedded decision paths
- `rules.css`: rule explorer
- `sources.css`: reference source cards

New feature-specific selectors belong in the corresponding feature file. Shared
tokens belong in `base.css`; shared form controls belong in `forms.css`.

## Tests

`tests/test_frontend_javascript.py` executes the JavaScript modules in Node and
checks strict composition, classification loading, manual-review prefill and the
amended-result request. Template, layout, security-header and batch-parser
contracts are covered by the other `test_frontend_*` and `test_ui_api_auth.py`
tests. External services are mocked in frontend behavior tests.
