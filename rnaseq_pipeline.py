#!/usr/bin/env python3
"""Run the pipeline straight from a source checkout (no installation needed).

    python3 rnaseq_pipeline.py            # interactive menu
    python3 rnaseq_pipeline.py --check    # health check

After `pip install .` / `pipx install .` the same program is available anywhere as `rnaseq-pipeline`.
"""
import sys
from pathlib import Path

if sys.version_info < (3, 10):
    sys.exit(f"[ERROR] Python 3.10+ is required (found {sys.version.split()[0]}).")

try:
    import yaml  # noqa: F401
except ImportError:
    sys.exit("[ERROR] PyYAML is missing. Install it with:  python3 -m pip install --user pyyaml\n"
             "        (or install the pipeline with pipx/conda, which brings it automatically)")

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from rnaseq_pipeline.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
