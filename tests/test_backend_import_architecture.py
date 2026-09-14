"""Architecture checks for explicit backend dependency wiring."""

from __future__ import annotations

import ast
import importlib
from pathlib import Path
import subprocess
import sys

from backend.domain.classification import ClassificationInputs


BACKEND_ROOT = Path(__file__).resolve().parents[1] / "backend"
PROJECT_ROOT = BACKEND_ROOT.parent


class _FunctionLocalBackendImportVisitor(ast.NodeVisitor):
    def __init__(self, relative_path: Path) -> None:
        self.relative_path = relative_path
        self.function_depth = 0
        self.violations: list[str] = []

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self.function_depth += 1
        self.generic_visit(node)
        self.function_depth -= 1

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        self.function_depth += 1
        self.generic_visit(node)
        self.function_depth -= 1

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        if self.function_depth and (node.module or "").startswith("backend"):
            self.violations.append(
                f"{self.relative_path}:{node.lineno}: from {node.module} import ..."
            )

    def visit_Import(self, node: ast.Import) -> None:
        if self.function_depth:
            for alias in node.names:
                if alias.name == "backend" or alias.name.startswith("backend."):
                    self.violations.append(
                        f"{self.relative_path}:{node.lineno}: import {alias.name}"
                    )


def test_backend_dependencies_are_not_imported_inside_functions() -> None:
    violations: list[str] = []
    for path in sorted(BACKEND_ROOT.rglob("*.py")):
        relative_path = path.relative_to(BACKEND_ROOT.parent)
        visitor = _FunctionLocalBackendImportVisitor(relative_path)
        visitor.visit(ast.parse(path.read_text(encoding="utf-8"), filename=str(path)))
        violations.extend(visitor.violations)

    assert not violations, (
        "Project-local backend imports must be declared at module level or bound in "
        "the production composition root. Hidden function-level dependencies found:\n"
        + "\n".join(violations)
    )


def test_classification_inputs_belong_to_domain_layer() -> None:
    assert ClassificationInputs.__module__ == "backend.domain.classification"


def _backend_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("backend"):
            imports.add(node.module or "")
        elif isinstance(node, ast.Import):
            imports.update(
                alias.name
                for alias in node.names
                if alias.name == "backend" or alias.name.startswith("backend.")
            )
    return imports


def _backend_package_edges() -> dict[str, set[str]]:
    packages = {
        path.name
        for path in BACKEND_ROOT.iterdir()
        if path.is_dir() and (path / "__init__.py").is_file()
    }
    edges = {package: set() for package in packages}
    for path in BACKEND_ROOT.rglob("*.py"):
        source = path.relative_to(BACKEND_ROOT).parts[0]
        if source not in packages:
            continue
        for imported in _backend_imports(path):
            parts = imported.split(".")
            if len(parts) > 1 and parts[1] in packages and parts[1] != source:
                edges[source].add(parts[1])
    return edges


def test_backend_package_dependencies_are_acyclic() -> None:
    edges = _backend_package_edges()
    reachable = {
        (source, dependency)
        for source, dependencies in edges.items()
        for dependency in dependencies
    }
    for _ in range(len(edges)):
        reachable |= {
            (source, target)
            for source, intermediate in reachable
            for candidate, target in reachable
            if intermediate == candidate
        }
    cycles = sorted(
        (source, target)
        for source, target in reachable
        if source < target and (target, source) in reachable
    )
    assert not cycles, "Backend package dependency cycles: " + ", ".join(
        f"{source} <-> {target}" for source, target in cycles
    )


def test_runtime_cache_owners_register_their_clearers() -> None:
    for module_name in (
        "backend.classification_runtime.cache",
        "backend.lookups.bayesdel",
        "backend.lookups.clingen",
        "backend.lookups.clinvar",
        "backend.lookups.spliceai",
    ):
        importlib.import_module(module_name)

    registry = importlib.import_module("backend.infrastructure.cache_registry")
    assert set(registry.registered_runtime_caches()) == {
        "bayesdel",
        "classification_results",
        "clingen_erepo",
        "clinvar",
        "spliceai",
    }


def test_domain_and_policy_do_not_depend_on_application_layers() -> None:
    forbidden = (
        "backend.api",
        "backend.classification_dag",
        "backend.criteria",
        "backend.lookups",
        "backend.presentation",
        "backend.reference_data",
        "backend.review",
        "backend.services",
    )
    violations: list[str] = []
    for package in ("domain", "policy"):
        for path in sorted((BACKEND_ROOT / package).glob("*.py")):
            for imported in sorted(_backend_imports(path)):
                if imported.startswith(forbidden):
                    violations.append(f"{path.name}: {imported}")
    assert not violations, "Shared layers import an application layer:\n" + "\n".join(violations)


def test_criteria_do_not_import_lookups_or_presentation() -> None:
    forbidden = (
        "backend.api",
        "backend.lookups",
        "backend.presentation",
        "backend.review",
    )
    violations: list[str] = []
    for path in sorted((BACKEND_ROOT / "criteria").glob("*.py")):
        for imported in sorted(_backend_imports(path)):
            if imported.startswith(forbidden):
                violations.append(f"{path.name}: {imported}")
    assert not violations, "Criterion modules contain an inverted dependency:\n" + "\n".join(violations)


def test_lookups_do_not_import_criteria_or_review_layers() -> None:
    forbidden = (
        "backend.api",
        "backend.criteria",
        "backend.presentation",
        "backend.review",
    )
    violations: list[str] = []
    for path in sorted((BACKEND_ROOT / "lookups").glob("*.py")):
        for imported in sorted(_backend_imports(path)):
            if imported.startswith(forbidden):
                violations.append(f"{path.name}: {imported}")
    assert not violations, "Lookup modules contain an inverted dependency:\n" + "\n".join(violations)


def test_reference_data_does_not_import_criteria() -> None:
    violations: list[str] = []
    for path in sorted((BACKEND_ROOT / "reference_data").glob("*.py")):
        for imported in sorted(_backend_imports(path)):
            if imported.startswith("backend.criteria"):
                violations.append(f"{path.name}: {imported}")
    assert not violations, "Reference-data modules import criteria:\n" + "\n".join(violations)


def test_infrastructure_does_not_import_cache_owners() -> None:
    forbidden = (
        "backend.classification_runtime",
        "backend.criteria",
        "backend.lookups",
        "backend.reference_data",
        "backend.review",
        "backend.services",
    )
    violations: list[str] = []
    for path in sorted((BACKEND_ROOT / "infrastructure").glob("*.py")):
        for imported in sorted(_backend_imports(path)):
            if imported.startswith(forbidden):
                violations.append(f"{path.name}: {imported}")
    assert not violations, "Infrastructure imports a cache owner:\n" + "\n".join(violations)


def test_variant_processing_does_not_import_lookup_or_rule_layers() -> None:
    forbidden = (
        "backend.api",
        "backend.classification_dag",
        "backend.criteria",
        "backend.lookups",
        "backend.presentation",
        "backend.review",
        "backend.services",
    )
    violations: list[str] = []
    for path in sorted((BACKEND_ROOT / "variant_processing").glob("*.py")):
        for imported in sorted(_backend_imports(path)):
            if imported.startswith(forbidden):
                violations.append(f"{path.name}: {imported}")
    assert not violations, "Variant processing contains an inverted dependency:\n" + "\n".join(violations)


def test_application_and_data_layers_do_not_import_api_transport() -> None:
    violations: list[str] = []
    for package in (
        "classification_dag",
        "classification_runtime",
        "population_frequency",
        "reference_data",
        "review",
        "services",
    ):
        for path in sorted((BACKEND_ROOT / package).rglob("*.py")):
            for imported in sorted(_backend_imports(path)):
                if imported.startswith("backend.api"):
                    violations.append(f"{path.relative_to(BACKEND_ROOT)}: {imported}")
    assert not violations, "A lower layer imports API transport:\n" + "\n".join(violations)


def test_reference_and_lookup_modules_do_not_read_or_create_files_on_import() -> None:
    modules = (
        "backend.reference_data.table4",
        "backend.reference_data.table9",
        "backend.reference_data.erepo_pvs1_rna",
        "backend.policy.spliceai_profile",
        "backend.lookups.founder_variants",
        "backend.lookups.coordinates",
        "backend.lookups.bayesdel",
        "backend.lookups.spliceai",
    )
    script = f"""
import builtins
from pathlib import Path

def blocked(*args, **kwargs):
    raise AssertionError("file I/O occurred during module import")

builtins.open = blocked
Path.read_text = blocked
Path.read_bytes = blocked
Path.mkdir = blocked

for module in {modules!r}:
    __import__(module)
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_missing_required_data_keeps_health_available_and_blocks_readiness() -> None:
    script = """
import os
from pathlib import Path

os.environ["ARIANE_UI_SESSION_SECRET"] = "test-session-secret-at-least-32-bytes-long"
import backend.reference_data.paths as paths
paths.TABLE9_PATH = Path("does-not-exist-table9.json")

from backend import main
from fastapi.testclient import TestClient

response = TestClient(main.app).get("/api/health")
assert response.status_code == 503, response.text
payload = response.json()
assert payload["status"] == "not_ready", payload
assert payload["ready"] is False, payload
assert any(
    item["component"] == "required classification datasets"
    and "Table 9" in item["detail"]
    for item in payload["startup"]["failures"]
), payload
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def _function_call_paths(path: Path, function_name: str) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    function = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function_name
    )

    def dotted_name(node: ast.AST) -> str:
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Attribute):
            parent = dotted_name(node.value)
            return f"{parent}.{node.attr}" if parent else node.attr
        return ""

    return {
        dotted_name(node.func)
        for node in ast.walk(function)
        if isinstance(node, ast.Call)
    }


def test_classification_controllers_delegate_application_workflow() -> None:
    controller_path = BACKEND_ROOT / "api" / "classification.py"
    single_calls = _function_call_paths(controller_path, "classify_variant")
    batch_calls = _function_call_paths(controller_path, "classify_batch")

    assert "self.service.classify_single" in single_calls
    assert "self.service.prepare_batch" in batch_calls
    assert "self.service.classify_batch" in batch_calls

    forbidden = {
        "CLASSIFICATION_CACHE.get",
        "CLASSIFICATION_CACHE.put",
        "CLASSIFICATION_USAGE.record",
        "VariantRequest.model_validate",
        "asyncio.gather",
        "classification_fingerprint",
        "get_gene_policy",
    }
    assert not (single_calls & forbidden)
    assert not (batch_calls & forbidden)


def test_native_dag_has_no_runtime_engine_selector() -> None:
    runtime = (BACKEND_ROOT / "classification_dag" / "runtime.py").read_text(
        encoding="utf-8"
    )
    assert "ClassifierEngineMode" not in runtime
    assert "get_configured_engine_mode" not in runtime
    assert "ARIANE_CLASSIFIER_ENGINE" not in runtime


def test_classifier_fingerprint_covers_all_classification_layers() -> None:
    from backend.classification_runtime import fingerprint

    covered = {
        path.relative_to(fingerprint.PROJECT_ROOT).as_posix()
        for path in fingerprint._CODE_DIRECTORIES
    }
    assert {
        "backend/classification_dag",
        "backend/classification_runtime",
        "backend/contracts",
        "backend/criteria",
        "backend/domain",
        "backend/lookups",
        "backend/policy",
        "backend/population_frequency",
        "backend/reference_data",
        "backend/review",
        "backend/services",
        "backend/variant_processing",
    } <= covered
    assert "backend/modules" not in covered


def test_manual_evidence_responsibilities_remain_split() -> None:
    review_root = BACKEND_ROOT / "review"
    expected = {
        "definitions.py",
        "strength.py",
        "validation.py",
        "service.py",
        "manual_evidence.py",
    }
    assert expected.issubset({path.name for path in review_root.glob("*.py")})

    for name in expected - {"manual_evidence.py"}:
        line_count = len((review_root / name).read_text(encoding="utf-8").splitlines())
        assert line_count < 600, f"backend/review/{name} has grown to {line_count} lines"

    facade_tree = ast.parse(
        (review_root / "manual_evidence.py").read_text(encoding="utf-8")
    )
    assert not any(
        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        for node in facade_tree.body
    )


def test_backend_root_contains_only_composition_entry_points() -> None:
    root_modules = {path.name for path in BACKEND_ROOT.glob("*.py")}
    assert root_modules == {"__init__.py", "bootstrap.py", "main.py", "version.py"}


def test_transport_and_contract_boundaries_are_explicit() -> None:
    api_modules = {path.name for path in (BACKEND_ROOT / "api").glob("*.py")}
    assert {
        "admin.py",
        "auth.py",
        "classification.py",
        "manual.py",
        "public.py",
        "review.py",
        "session.py",
        "system.py",
    }.issubset(api_modules)
    assert {
        "batch.py",
        "client.py",
        "ps1.py",
        "result.py",
        "review.py",
        "variant.py",
    }.issubset({path.name for path in (BACKEND_ROOT / "contracts").glob("*.py")})

    retired = {
        "admin.py",
        "api_auth.py",
        "config.py",
        "data_health.py",
        "data_validation.py",
        "models.py",
        "public_api.py",
        "review_api.py",
        "review_records.py",
        "runtime_cache.py",
        "runtime_data.py",
        "spliceai_profile.py",
        "startup.py",
        "ui_session.py",
    }
    assert not (retired & {path.name for path in BACKEND_ROOT.glob("*.py")})


def test_main_is_a_composition_root_not_a_route_collection() -> None:
    main_path = BACKEND_ROOT / "main.py"
    lines = main_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) < 180
    source = "\n".join(lines)
    assert "create_manual_router" in source
    assert "create_system_router" in source
    assert "install_http_middleware" in source
    assert "install_exception_handlers" in source
    assert "install_frontend" in source
    assert "@app." not in source


def test_paired_routes_use_one_auth_policy_helper() -> None:
    routing = (BACKEND_ROOT / "api" / "routing.py").read_text(encoding="utf-8")
    assert "require_public_api_key" in routing
    assert "require_ui_session" in routing
    assert "include_in_schema=False" in routing

    for module_name in ("classification.py", "manual.py", "system.py"):
        source = (BACKEND_ROOT / "api" / module_name).read_text(encoding="utf-8")
        assert '"/ui-api/' not in source
    assert '@paired.post("/manual-evidence/evaluate")' in (
        BACKEND_ROOT / "api" / "manual.py"
    ).read_text(encoding="utf-8")
    assert '@paired.get("/rules")' in (
        BACKEND_ROOT / "api" / "system.py"
    ).read_text(encoding="utf-8")


def test_runtime_health_has_no_process_global_registry() -> None:
    health_path = BACKEND_ROOT / "infrastructure" / "health.py"
    tree = ast.parse(health_path.read_text(encoding="utf-8"))
    assert not any(
        isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "_ISSUES" for target in node.targets)
        for node in tree.body
    )
    assert "DataHealthRegistry()" not in health_path.read_text(encoding="utf-8")
