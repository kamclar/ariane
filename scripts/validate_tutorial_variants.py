"""Run the tutorial variant panel through the production classification path."""
import asyncio
import json
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.main import _classify_one


# Expected ENIGMA classifications and evidence are transcribed from the local
# Run 13-15 webinar transcripts. None means that the variant came from the
# comparison spreadsheet, not from an explanatory tutorial transcript.
CASES = [
    ("run14", "BRCA1", "c.509G>A", "p.(Arg170Gln)", 1, {"BS1_Supporting", "BP1", "BS3"}),
    ("run14", "BRCA1", "c.1534C>T", "p.(Leu512Phe)", 1, {"BS1_Supporting", "BP1", "BS3", "BP5"}),
    ("run14", "BRCA1", "c.3668_3671dup", "p.(Cys1225fs)", 5, {"PVS1", "PM2_Supporting", "PM5_PTC"}),
    ("run14", "BRCA2", "c.9097del", "p.(Thr3033fs)", 5, {"PVS1", "PM2_Supporting", "PM5_PTC"}),
    ("run15", "BRCA1", "c.5551_5552insT", "p.(Asp1851ValfsTer29)", 3, {"PVS1"}),
    ("run15", "BRCA2", "c.(793+1_794-1)_(1909+1_1910-1)del", "p.(?)", 3, {"PM2_Supporting", "BS3"}),
    ("run15", "BRCA2", "c.6147_6149del", "p.(Val2050del)", 2, {"BP1"}),
    ("run13", "BRCA1", "c.3891_3893del", "p.(Ser1298del)", 1, {"BS3", "BP1", "BP5"}),
    ("run13", "BRCA1", "c.4185G>A", "p.(Gln1395=)", 5, {"PVS1_RNA", "PM2_Supporting", "PP4"}),
    ("run13", "BRCA1", "c.628C>T", "p.(Gln210Ter)", 3, {"PM2_Supporting"}),
    ("run13", "BRCA2", "c.8953+2T>C", "p.(?)", 3, {"PM2_Supporting"}),
    ("clinvar", "BRCA1", "c.3247A>C", "p.(Met1083Leu)", None, None),
    ("extra", "BRCA1", "c.5217T>A", "p.(Asp1739Glu)", None, None),
]


async def main() -> None:
    unavailable = {"status": "not_found"}
    with patch("backend.lookups.clinvar.clinvar_lookup", return_value=unavailable), patch(
        "backend.lookups.clingen.clingen_erepo_lookup", return_value=unavailable
    ):
        for run, gene, c_notation, p_notation, expected_class, expected_criteria in CASES:
            try:
                result = await _classify_one(gene, c_notation, p_notation, "Unknown")
                record = {
                    "run": run,
                    "variant": f"{gene} {c_notation} {p_notation}",
                    "class": result.predicted_class,
                    "label": result.predicted_label,
                    "points": result.total_points,
                    "criteria": [
                        f"{item.name}:{item.strength}:{item.points}"
                        for item in result.criteria
                    ],
                    "material_warnings": [
                        warning for warning in result.warnings
                        if not warning.startswith("FIRST PASS")
                        and not warning.startswith("Data source degraded")
                    ],
                }
                if expected_class is not None:
                    actual_criteria = {item.name for item in result.criteria}
                    record["transcript_expected_class"] = expected_class
                    record["class_matches_transcript"] = result.predicted_class == expected_class
                    record["missing_transcript_criteria"] = sorted(expected_criteria - actual_criteria)
                    record["extra_vs_transcript"] = sorted(actual_criteria - expected_criteria)
            except Exception as exc:
                record = {
                    "run": run,
                    "variant": f"{gene} {c_notation} {p_notation}",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            print(json.dumps(record, ensure_ascii=False))


if __name__ == "__main__":
    asyncio.run(main())
