"""Intentional-failure tests: the pipeline must stop, explain, log, and never accept broken outputs."""
import json
import os
import subprocess
import threading
from pathlib import Path
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
        runner.run(["sleep", "30.4917"], stage="t")    # a duration no other process on the machine uses
    assert time.time() - t0 < 15
    assert subprocess.run(["pgrep", "-f", r"^sleep 30\.4917$"], capture_output=True).returncode != 0


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


def _bam(tmp_path, name, n):
    sam = tmp_path / f"{name}.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:5000\n"
                   + "".join(f"r{i}\t0\tc1\t{10 + i}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" for i in range(n)))
    bam = tmp_path / f"{name}.bam"
    subprocess.run(["samtools", "view", "-b", "-o", bam, sam], check=True)
    return bam


@pytest.mark.skipif(not have("samtools"), reason="samtools not installed")
def test_bam_index_judged_by_content_not_clock(tmp_path):
    """P1/H8: the system clock can step backwards (WSL2 time sync, NTP) or files can be copied without their
    times, so an index that merely LOOKS older must not fail a correct BAM; a stale index that LOOKS newer must."""
    import os
    bam = _bam(tmp_path, "x", 2000)
    subprocess.run(["samtools", "index", bam], check=True)
    bai = Path(str(bam) + ".bai")
    t = bam.stat().st_mtime
    os.utime(bai, (t - 10, t - 10))                    # observed: the clock stepped back ~1 s between the two writes
    ok, probs, _ = bam_manager.validate_bam(bam, 2000, paired=False)
    assert ok, probs
    # stale: index of a 2000-record BAM next to a 1500-record BAM, with a newer timestamp
    other = _bam(tmp_path, "y", 1500)
    os.replace(other, bam)
    os.utime(bai, (t + 10, t + 10))
    ok, probs, _ = bam_manager.validate_bam(bam, 1500, paired=False)
    assert not ok and any("index" in x and "does not match" in x for x in probs), probs


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


# ---------------------------------------------------------------- F2 regression: SIGTERM / SIGHUP
_SIGNAL_CHILD = r"""
import sys, time
sys.path.insert(0, {src!r})
from rnaseq_pipeline import install_signal_handlers, runner, logger, Terminated
logger.attach_project({logs!r})
install_signal_handlers()
try:
    runner.run_pipeline([["sleep", "300.4917"], ["cat"]], stage="t", stdout_file={out!r})
except Terminated as e:
    print("TERMINATED", e.signame, flush=True)
    sys.exit(143)
"""


@pytest.mark.parametrize("signame", ["SIGTERM", "SIGHUP"])
def test_termination_signal_stops_child_processes(tmp_path, signame):
    """Before the fix, SIGTERM killed the pipeline but left its tools running (they live in their own session)."""
    import signal
    import sys
    from conftest import ROOT
    code = _SIGNAL_CHILD.format(src=str(ROOT / "src"), logs=str(tmp_path / "logs"), out=str(tmp_path / "o.txt"))
    proc = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, text=True)
    deadline = time.time() + 10
    while time.time() < deadline and subprocess.run(["pgrep", "-f", r"^sleep 300\.4917$"], capture_output=True).returncode:
        time.sleep(0.1)
    assert subprocess.run(["pgrep", "-f", r"^sleep 300\.4917$"], capture_output=True).returncode == 0, "child never started"
    proc.send_signal(getattr(signal, signame))
    out, _ = proc.communicate(timeout=20)
    assert proc.returncode == 143 and f"TERMINATED {signame}" in out
    time.sleep(0.5)
    assert subprocess.run(["pgrep", "-f", r"^sleep 300\.4917$"], capture_output=True).returncode != 0, "orphaned child"
    rec = (tmp_path / "logs" / "commands.jsonl").read_text()
    assert '"interrupted"' in rec


def test_fastq_validation_pool_stops_on_interrupt(tmp_path):
    """Worker processes validating large files must be terminated at once, not after finishing the file."""
    import gzip
    import multiprocessing
    import signal
    import threading
    from rnaseq_pipeline import fastq_validator as FV
    block = "".join(f"@r{n}\n" + "ACGT" * 25 + "\n+\n" + "I" * 100 + "\n" for n in range(100_000))
    jobs = []
    for i in range(2):  # 2 x 2,000,000 reads: uninterrupted validation takes >10 s (measured), so 5 s proves
        f = tmp_path / f"s{i}.fq.gz"  # the workers were terminated rather than allowed to finish
        with gzip.open(f, "wt", compresslevel=1) as fh:
            for _ in range(20):
                fh.write(block)
        jobs.append((f"s{i}", str(f), None, "ACGTNacgtn.", 33))
    threading.Timer(0.5, lambda: os.kill(os.getpid(), signal.SIGINT)).start()
    t0 = time.time()
    with pytest.raises(KeyboardInterrupt):
        FV.validate_many(jobs, workers=2)
    assert time.time() - t0 < 5, "validation workers were not terminated promptly"
    assert not multiprocessing.active_children()


# ---------------------------------------------------------------- P1/H6: intermittent hisat2-build crash
@pytest.mark.skipif(not have("hisat2-build", "hisat2-inspect", "hisat2_extract_splice_sites.py"),
                    reason="HISAT2 not installed")
@pytest.mark.parametrize("mode,calls,expect", [
    ("crash-when-multithreaded", 2, None),        # the observed hisat2-build 2.2.3 race: retried single-threaded
    ("crash-always", 2, "SIGSEGV"),               # a crash that recurs is reported after one retry
    ("out-of-memory", 1, "SIGKILL"),              # never retried: it would fail again
])
def test_index_build_crash_handling(tmp_path, monkeypatch, mode, calls, expect):
    import random
    import shutil as sh
    from rnaseq_pipeline import reference_manager as R
    real = sh.which("hisat2-build")
    fake = tmp_path / "bin" / "hisat2-build"
    fake.parent.mkdir()
    log = tmp_path / "calls.txt"
    action = {"crash-when-multithreaded": f'if [ "$2" != "1" ]; then kill -SEGV $$; fi\nexec "{real}" "$@"',
              "crash-always": "kill -SEGV $$", "out-of-memory": "kill -KILL $$"}[mode]
    fake.write_text(f'#!/bin/sh\necho "$@" >> "{log}"\n{action}\n')
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setitem(runner._tool_env, "path_prefix", [])
    rng = random.Random(3)
    genome = "".join(rng.choice("ACGT") for _ in range(3000))
    paths = R.RefPaths(tmp_path / "ref")
    paths.genome.parent.mkdir(parents=True)
    paths.gtf.parent.mkdir(parents=True)
    paths.genome.write_text(">c1\n" + "\n".join(genome[i:i + 60] for i in range(0, 3000, 60)) + "\n")
    paths.gtf.write_text('c1\tt\texon\t201\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                         'c1\tt\texon\t1001\t1600\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    args = (paths, {"names": ["c1"]}, 4, 8, True, tmp_path / "ref.log")
    if expect is None:
        R.build_index(*args)
        assert R.index_valid(paths.index_prefix, ["c1"])[0]
    else:
        with pytest.raises(PipelineError) as e:
            R.build_index(*args)
        assert expect in str(e.value)
        assert not paths.index_dir.exists()          # nothing partial is ever put in place
    lines = log.read_text().splitlines()
    assert len(lines) == calls, lines
    if calls == 2:
        assert lines[0].startswith("-p 4") and lines[1].startswith("-p 1")


def test_interrupted_informational_command_leaves_no_orphans():
    """P1/H10: quick commands (version probes, flagstat, quickcheck) run through runner.tool_output/tool_run. Many
    tools are wrapper scripts (hisat2, fastqc) whose real program is a grandchild; Ctrl-C/SIGTERM must stop it too."""
    import signal

    def interrupt():
        time.sleep(1)
        os.kill(os.getpid(), signal.SIGINT)
    for fn in (runner.tool_output, runner.tool_run):
        threading.Thread(target=interrupt, daemon=True).start()
        t0 = time.time()
        with pytest.raises(KeyboardInterrupt):
            fn(["sh", "-c", "sleep 30.4919; true"], timeout=60)      # the wrapper's child is the real work
        assert time.time() - t0 < 15
        time.sleep(0.3)
        assert subprocess.run(["pgrep", "-f", r"^sleep 30\.4919$"], capture_output=True).returncode != 0, fn.__name__
