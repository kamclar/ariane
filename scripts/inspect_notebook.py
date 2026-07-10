import argparse
import json
import re
from pathlib import Path


def first_line(cell):
    source = "".join(cell.get("source", []))
    return next((line.strip() for line in source.splitlines() if line.strip()), "")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("notebook", type=Path)
    parser.add_argument("--cell", type=int, action="append", default=[])
    parser.add_argument("--functions", action="store_true")
    args = parser.parse_args()

    notebook = json.loads(args.notebook.read_text(encoding="utf-8"))
    cells = notebook["cells"]
    print(f"cells={len(cells)} nbformat={notebook.get('nbformat')}.{notebook.get('nbformat_minor')}")

    if args.functions:
        for index, cell in enumerate(cells):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            functions = re.findall(r"^def\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(", source, re.MULTILINE)
            if functions:
                print(f"{index:03d}: {', '.join(functions)}")
        return

    if args.cell:
        for index in args.cell:
            cell = cells[index]
            print(f"\n===== CELL {index:03d} {cell.get('cell_type')} =====")
            print("".join(cell.get("source", [])))
        return

    for index, cell in enumerate(cells):
        print(f"{index:03d} {cell.get('cell_type'):8s} {first_line(cell)[:150]}")


if __name__ == "__main__":
    main()
