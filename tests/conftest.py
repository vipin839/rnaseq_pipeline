import gzip
import os
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

# Scientific tools: $RNASEQ_TOOLS_BIN, else the default conda env, else whatever is on PATH (e.g. CI env).
TOOLS_BIN = Path(os.environ.get("RNASEQ_TOOLS_BIN", Path.home() / "miniforge3" / "envs" / "rnaseq-tools" / "bin"))
if TOOLS_BIN.is_dir():
    os.environ["PATH"] = f"{TOOLS_BIN}{os.pathsep}{os.environ['PATH']}"


def _find_rscript():
    cand = os.environ.get("RNASEQ_RSCRIPT") or str(Path.home() / "miniforge3" / "envs" / "rnaseq-r" / "bin" / "Rscript")
    if Path(cand).is_file():
        return Path(cand)
    found = shutil.which("Rscript")
    return Path(found) if found else None


RSCRIPT = _find_rscript()


def have(*tools):
    return all(shutil.which(t) for t in tools)


def write_fastq(path, records, gz=True):
    """records: list of (id, seq, qual)."""
    opener = gzip.open if gz else open
    with opener(path, "wt") as f:
        for rid, seq, qual in records:
            f.write(f"@{rid}\n{seq}\n+\n{qual}\n")
    return path


@pytest.fixture
def good_pair(tmp_path):
    recs1 = [(f"r{i}/1", "ACGT" * 10, "I" * 40) for i in range(100)]
    recs2 = [(f"r{i}/2", "TTGA" * 10, "I" * 40) for i in range(100)]
    return (write_fastq(tmp_path / "S1_R1.fastq.gz", recs1), write_fastq(tmp_path / "S1_R2.fastq.gz", recs2))
