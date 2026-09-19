import gzip
import os
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

TOOLS_BIN = Path.home() / "miniforge3" / "envs" / "rnaseq-tools" / "bin"
if TOOLS_BIN.is_dir():
    os.environ["PATH"] = f"{TOOLS_BIN}{os.pathsep}{os.environ['PATH']}"


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
