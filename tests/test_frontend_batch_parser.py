import json
import shutil
import subprocess
from pathlib import Path

import pytest


BATCH_JS = (
    Path(__file__).resolve().parents[1] / "frontend" / "static" / "js" / "batch.js"
)


def parse_batch(text: str) -> dict:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node.js is required to execute the frontend batch parser test")
    script = r"""
const fs = require("fs");
global.window = { ArianeFrontend: {} };
eval(fs.readFileSync(process.argv[1], "utf8"));
const input = JSON.parse(fs.readFileSync(0, "utf8"));
const result = window.ArianeFrontend.parseBatchInput(input.text, input.genes);
process.stdout.write(JSON.stringify(result));
"""
    completed = subprocess.run(
        [node, "-e", script, str(BATCH_JS)],
        input=json.dumps(
            {
                "text": text,
                "genes": [{"symbol": "BRCA1"}, {"symbol": "BRCA2"}],
            }
        ),
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(completed.stdout)


def test_space_separated_transcript_qualified_batch_is_parsed():
    text = """BRCA1 NM_007294.4:c.509G>A p.(Arg170Gln)
BRCA1 NM_007294.4:c.1534C>T p.(Leu512Phe)
BRCA2 NM_000059.4:c.9097del p.(Thr3033fs)
BRCA1 NM_007294.4:c.5551_5552insT p.(Asp1851ValfsTer29)
BRCA2 NM_000059.4:c.6147_6149del p.(Val2050del)
BRCA1 NM_007294.4:c.3891_3893del p.(Ser1298del)
BRCA1 NM_007294.4:c.4185G>A p.(Gln1395=)
BRCA1 NM_007294.4:c.628C>T p.(Gln210Ter)
BRCA2 NM_000059.4:c.8953+2T>C p.(?)
BRCA1 NM_007294.4:c.3247A>C p.(Met1083Leu)"""

    result = parse_batch(text)

    assert result["errors"] == []
    assert len(result["parsed"]) == 10
    assert result["parsed"][0] == {
        "gene": "BRCA1",
        "c_notation": "NM_007294.4:c.509G>A",
        "p_notation": "p.(Arg170Gln)",
        "assembly": "",
        "dup_type": "Unknown",
    }
    assert result["parsed"][8]["p_notation"] == "p.(?)"


def test_comma_and_space_formats_can_be_mixed_with_genomic_assembly():
    result = parse_batch(
        "BRCA1, c.509G>A, p.(Arg170Gln)\n"
        "BRCA1 c.1534C>T p.Leu512Phe\n"
        "BRCA1 chr17:43099813:C>T hg38"
    )

    assert result["errors"] == []
    assert [item["p_notation"] for item in result["parsed"]] == [
        "p.(Arg170Gln)",
        "p.(Leu512Phe)",
        "",
    ]
    assert result["parsed"][2]["assembly"] == "GRCh38"


def test_markdown_wrappers_and_escaped_underscores_are_tolerated():
    result = parse_batch(
        "*BRCA1 NM\\_007294.4:c.5551\\_5552insT "
        "p.(Asp1851ValfsTer29)*\\"
    )

    assert result["errors"] == []
    assert result["parsed"][0]["c_notation"] == "NM_007294.4:c.5551_5552insT"
