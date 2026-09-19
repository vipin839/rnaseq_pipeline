#!/usr/bin/env python3
"""Bulk RNA-seq Pipeline v1 — main entry point.

    python3 rnaseq_pipeline.py            # interactive menu
    python3 rnaseq_pipeline.py --dry-run  # show what would run
    python3 rnaseq_pipeline.py --help
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

if sys.version_info < (3, 10):
    sys.exit(f"[ERROR] Python 3.10+ is required (found {sys.version.split()[0]}).")

try:
    import yaml  # noqa: F401
except ImportError:
    # fall back to the pipeline's own conda env python, which has PyYAML
    for cand in (Path.home() / "miniforge3/envs/rnaseq-tools/bin/python",):
        if cand.exists() and os.environ.get("RNASEQ_REEXEC") != "1":
            os.environ["RNASEQ_REEXEC"] = "1"
            os.execv(str(cand), [str(cand), __file__, *sys.argv[1:]])
    sys.exit("[ERROR] PyYAML is missing. Install it with:  python3 -m pip install --user pyyaml")

sys.path.insert(0, str(ROOT / "python"))

from rnaseq.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
