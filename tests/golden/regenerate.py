"""Regenerate the golden output for the synthetic transcriptome.

Run from the repository root::

    python -m tests.golden.regenerate

The golden file pins RiboMetric's output on ``tests/synthetic.py`` so that CI
fails when any change alters it. An intentional change regenerates the file
and says in its commit message which numbers moved and why
(docs/PRODUCTIONISATION.md, Phase 1).
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

from tests import synthetic

GOLDEN = Path(__file__).parent / "synthetic_tripsviz_a_site.json"
RUN = "import sys; from RiboMetric.cli import main; sys.exit(main())"
ARGS = ["--offset-calculation-method", "tripsviz", "--offset-target", "a_site"]
KEYS = (
    "mode",
    "metrics",
    "scores",
    "computed_offsets",
    "read_frame_distribution",
    "read_length_distribution",
    "recommended_read_lengths",
    "mRNA_distribution",
)


def run_synthetic(workdir: Path) -> dict:
    bam, annotation = synthetic.build_fixture(workdir / "fixture")
    out = workdir / "out"
    cmd = [sys.executable, "-c", RUN, "run", "-b", str(bam), "-a", str(annotation), *ARGS]
    cmd += ["--threads", "2", "--json", "-o", str(out), "-n", "synth"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return json.loads((out / "synth_RiboMetric.json").read_text())["results"]


def golden_view(results: dict) -> dict:
    return {key: results.get(key) for key in KEYS}


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        view = golden_view(run_synthetic(Path(tmp)))
    GOLDEN.write_text(json.dumps(view, indent=1, sort_keys=True) + "\n")
    print(f"wrote {GOLDEN}")


if __name__ == "__main__":
    main()
