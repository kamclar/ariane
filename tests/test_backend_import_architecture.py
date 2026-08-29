"""Architecture checks for explicit backend dependency wiring."""

from __future__ import annotations

import ast
from pathlib import Path

from backend.classification_dag.domain import ClassificationInputs


BACKEND_ROOT = Path(__file__).resolve().parents[1] / "backend"


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
    assert ClassificationInputs.__module__ == "backend.classification_dag.domain"

