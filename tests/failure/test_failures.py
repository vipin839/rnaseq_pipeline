"""Intentional-failure tests: the pipeline must stop, explain, log, and never accept broken outputs."""
import json
import subprocess
import threading
import time

import pytest

from conftest import RSCRIPT, have
from rnaseq_pipeline import PipelineError, bam_manager, logger, runner, storage
from rnaseq_pipeline import config as C


@pytest.fixture(autouse=True)
def logs(tmp_path):
    logger.attach_project(tmp_path / "logs")
    yield tmp_path / "logs"
    logger.detach_project()


def test_missing_tool(tmp_path):
    with pytest.raises(PipelineError, match="program not found"):
        runner.run(["definitely-not-a-tool-xyz", "--version"], stage="t")


def test_pipefail_semantics(tmp_path, logs):
    # upstream fails, downstream succeeds -> must be reported as failure (like set -o pipefail)
    with pytest.raises(PipelineError) as e:
        runner.run_pipeline([["false"], ["cat"]], stage="t", stdout_file=tmp_path / "o.txt",
                            log_file=tmp_path / "l.log")
    assert "exit codes [1, 0]" in str(e.value)
    rec = [json.loads(l) for l in (logs / "commands.jsonl").read_text().splitlines()]
    assert rec[-1]["status"] == "failed" and rec[-1]["exit_codes"] == [1, 0]
    assert "false | cat" in (logs / "command_history.log").read_text()


def test_empty_argument_rejected():
    with pytest.raises(PipelineError, match="empty argument"):
        runner.run(["echo", ""], stage="t")


def test_no_shell_expansion(tmp_path):
    out = tmp_path / "o.txt"
    runner.run(["echo", "$(id) ; rm -rf / *"], stage="t", stdout_file=out)
    assert out.read_text().strip() == "$(id) ; rm -rf / *"


def test_timeout_kills_command(tmp_path):
    t0 = time.time()
    with pytest.raises(PipelineError, match="timed out"):
        runner.run(["sleep", "30"], stage="t", timeout=1)
    assert time.time() - t0 < 15


def test_interrupted_command_is_terminated(tmp_path):
    """Ctrl-C during a command: child process group is killed and the command is logged as interrupted."""
    import os
    import signal

    def interrupt():
        time.sleep(1)
        os.kill(os.getpid(), signal.SIGINT)  # what the terminal sends on Ctrl-C

    threading.Thread(target=interrupt, daemon=True).start()
    t0 = time.time()
    with pytest.raises(KeyboardInterrupt):
        runner.run(["sleep", "30"], stage="t")
    assert time.time() - t0 < 15
    assert subprocess.run(["pgrep", "-f", "^sleep 30$"], capture_output=True).returncode != 0


@pytest.mark.skipif(not have("samtools"), reason="samtools not installed")
def test_corrupted_bam_detected(tmp_path):
    sam = tmp_path / "x.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n"
                   + "".join(f"r{i}\t0\tc1\t{10 + i}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" for i in range(2000)))
    bam = tmp_path / "x.bam"
    subprocess.run(["samtools", "view", "-b", "-o", bam, sam], check=True)
    subprocess.run(["samtools", "index", bam], check=True)
    ok, probs, m = bam_manager.validate_bam(bam, 2000, paired=False)
    assert ok, probs
    ok, probs, _ = bam_manager.validate_bam(bam, 2500, paired=False)  # incomplete relative to input
    assert not ok and "incomplete" in probs[0]
    data = bam.read_bytes()
    bam.write_bytes(data[: len(data) - 40])  # remove EOF block -> truncated
    ok, probs, _ = bam_manager.validate_bam(bam, 2000, paired=False)
    assert not ok and "quickcheck" in probs[0]


def test_insufficient_disk(tmp_path, monkeypatch):
    assert storage.preflight(tmp_path, 10 ** 9, C.load()) is False


def test_invalid_statistical_thresholds():
    cfg = C.load()
    cfg["alpha"] = 5
    cfg["log2fc_threshold"] = -2
    with pytest.raises(PipelineError, match="alpha"):
        C.require_valid(cfg, 4)


@pytest.mark.skipif(RSCRIPT is None, reason="R env not installed")
@pytest.mark.parametrize("case", ["missing_sample", "confounded", "no_replicates", "bad_level"])
def test_r_rejects_bad_design(tmp_path, case):
    from rnaseq_pipeline import r_bridge
    counts = tmp_path / "c.tsv"
    counts.write_text("Gene_ID\tA\tB\tC\tD\ng1\t10\t20\t30\t40\ng2\t5\t6\t7\t8\n")
    meta = tmp_path / "m.tsv"
    rows = {"missing_sample": "sample\tcondition\nA\tC\nB\tC\nC\tT\n",
            "confounded": "sample\tbatch\tcondition\nA\tb1\tC\nB\tb1\tC\nC\tb2\tT\nD\tb2\tT\n",
            "no_replicates": "sample\tcondition\nA\tC\nB\tC\nC\tC\nD\tT\n",
            "bad_level": "sample\tcondition\nA\tC\nB\tC\nC\t1T\nD\t1T\n"}[case]
    meta.write_text(rows)
    formula = "~ batch + condition" if case == "confounded" else "~ condition"
    params = {"counts": str(counts), "metadata": str(meta), "formula": formula,
              "variables": ["batch", "condition"] if case == "confounded" else ["condition"],
              "variable_of_interest": "condition", "reference_level": "C", "contrasts": [["T", "C"]]}
    pp = tmp_path / "p.json"
    pp.write_text(json.dumps(params))
    rs = RSCRIPT
    with pytest.raises(PipelineError, match="R/DESeq2 step failed"):
        r_bridge.run(rs, pp, tmp_path / "r.log", validate_only=True)
