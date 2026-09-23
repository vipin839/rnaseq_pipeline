"""F3/F4/F5 regression: downloads are never accepted just because a file exists."""
import hashlib
import http.server
import os
import threading
from pathlib import Path

import pytest

from rnaseq_pipeline import PipelineError, logger, net

DATA = os.urandom(2_000_000)
MD5 = hashlib.md5(DATA).hexdigest()


class Handler(http.server.BaseHTTPRequestHandler):
    ranges = True  # toggled per test

    def log_message(self, *a):
        pass

    stalls_left = 0          # /stall.bin: send part of the file, then stop sending (a dead connection)
    stall_seconds = 120      # far longer than the abort + resume path (~13 s), so timing cannot be ambiguous

    def do_GET(self):
        if self.path in ("/stall.bin", "/flaky.bin", "/dead.bin"):
            return self._partial()
        body = {"/data.bin": DATA, "/empty.bin": b""}.get(self.path)
        if body is None:
            self.send_error(404)
            return
        rng = self.headers.get("Range")
        if rng and Handler.ranges:
            start = int(rng.split("=")[1].split("-")[0])
            chunk = body[start:]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(body) - 1}/{len(body)}")
        else:
            chunk = body
            self.send_response(200)
        self.send_header("Content-Length", str(len(chunk)))
        self.end_headers()
        self.wfile.write(chunk)


    def _partial(self):
        import time
        rng = self.headers.get("Range")
        start = int(rng.split("=")[1].split("-")[0]) if rng else 0
        if self.path == "/dead.bin":            # the connection breaks before any data arrives, every time
            self.connection.close()
            return
        chunk = DATA[start:]
        self.send_response(206 if start else 200)
        if start:
            self.send_header("Content-Range", f"bytes {start}-{len(DATA) - 1}/{len(DATA)}")
        self.send_header("Content-Length", str(len(chunk)))
        self.end_headers()
        if self.path == "/stall.bin" and Handler.stalls_left > 0:
            Handler.stalls_left -= 1
            self.wfile.write(chunk[:300_000])
            self.wfile.flush()
            time.sleep(Handler.stall_seconds)    # nothing more arrives; the socket stays open
            return
        if self.path == "/flaky.bin":           # each connection delivers 300 kB, then drops (progress each time)
            self.wfile.write(chunk[:300_000])
            return
        self.wfile.write(chunk)


@pytest.fixture(scope="module")
def server():
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()


@pytest.fixture(autouse=True)
def logs(tmp_path):
    logger.attach_project(tmp_path / "logs")
    Handler.ranges = True
    yield
    logger.detach_project()


def dl(url, dest, **kw):
    return net.download(url, dest, stage="t", log_file=dest.parent / "dl.log", retries=2, **kw)


def test_verified_download(server, tmp_path):
    out = dl(f"{server}/data.bin", tmp_path / "x.bin", expected_md5=MD5, expected_bytes=len(DATA))
    assert out.read_bytes() == DATA and not (tmp_path / "x.bin.part").exists()


def test_resume_from_partial(server, tmp_path):
    (tmp_path / "x.bin.part").write_bytes(DATA[:700_000])
    out = dl(f"{server}/data.bin", tmp_path / "x.bin", expected_md5=MD5)
    assert out.read_bytes() == DATA


def test_restart_when_server_cannot_resume(server, tmp_path):
    """F3: before the fix, curl exit 33 aborted the download for good."""
    Handler.ranges = False
    (tmp_path / "x.bin.part").write_bytes(DATA[:700_000])
    out = dl(f"{server}/data.bin", tmp_path / "x.bin", expected_md5=MD5)
    assert out.read_bytes() == DATA


def test_empty_download_rejected(server, tmp_path):
    """F4: before the fix, a 0-byte file without a checksum was accepted."""
    with pytest.raises(PipelineError, match="empty"):
        dl(f"{server}/empty.bin", tmp_path / "x.bin")
    assert not (tmp_path / "x.bin").exists() and (tmp_path / "x.bin.part.rejected").exists()


def test_size_mismatch_rejected(server, tmp_path):
    with pytest.raises(PipelineError, match="size is"):
        dl(f"{server}/data.bin", tmp_path / "x.bin", expected_bytes=len(DATA) + 1)
    assert not (tmp_path / "x.bin").exists()


def test_checksum_mismatch_rejected_and_kept(server, tmp_path):
    with pytest.raises(PipelineError, match="MD5"):
        dl(f"{server}/data.bin", tmp_path / "x.bin", expected_md5="0" * 32)
    assert not (tmp_path / "x.bin").exists() and (tmp_path / "x.bin.part.rejected").exists()


def test_http_404_explained_and_not_retried(server, tmp_path):
    """F5: plain-language reason; a permanent error is not retried."""
    with pytest.raises(PipelineError) as e:
        dl(f"{server}/missing.bin", tmp_path / "x.bin")
    assert "HTTP error" in str(e.value) and "exit code" not in str(e.value)
    assert (tmp_path / "dl.log").read_text().count("### exit codes") == 1


def test_unreachable_server_explained(tmp_path):
    with pytest.raises(PipelineError, match="could not connect"):
        dl("http://127.0.0.1:9/x.bin", tmp_path / "x.bin")


def test_existing_invalid_file_is_never_overwritten(server, tmp_path):
    (tmp_path / "x.bin").write_bytes(b"garbage")
    with pytest.raises(PipelineError, match="not valid"):
        dl(f"{server}/data.bin", tmp_path / "x.bin", expected_md5=MD5)
    assert (tmp_path / "x.bin").read_bytes() == b"garbage"


# ---------------------------------------------------------------- F18: pilot streaming (first N reads)
import gzip as _gzip  # noqa: E402

FQ = "".join(f"@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(1000)).encode()


class PilotHandler(http.server.BaseHTTPRequestHandler):
    fail_first = 0      # number of requests to answer with HTTP 500
    requests = 0

    def log_message(self, *a):
        pass

    def do_GET(self):
        PilotHandler.requests += 1
        if PilotHandler.requests <= PilotHandler.fail_first:
            self.send_error(500)
            return
        body = {"/reads.fq.gz": _gzip.compress(FQ), "/cut.fq.gz": _gzip.compress(FQ)[:300]}.get(self.path)
        if body is None:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


@pytest.fixture
def pilot_server(monkeypatch):
    from rnaseq_pipeline import ena_manager
    monkeypatch.setattr(ena_manager, "RETRY_DELAY", 0)
    PilotHandler.requests, PilotHandler.fail_first = 0, 0
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), PilotHandler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()


def _pilot(url, dest, n):
    from rnaseq_pipeline import ena_manager
    return ena_manager._download_head(url, dest, n, "SRRTEST", dest.parent / "dl.log")


def _records(path):
    return _gzip.open(path, "rt").read().count("\n") // 4


def test_pilot_first_n_reads(pilot_server, tmp_path):
    out = _pilot(f"{pilot_server}/reads.fq.gz", tmp_path / "x.fq.gz", 100)
    assert _records(out) == 100


def test_pilot_retries_transient_server_errors(pilot_server, tmp_path):
    """Before the fix, one failed request = 'incomplete (0 lines)' and the sample was marked FAILED."""
    PilotHandler.fail_first = 2
    out = _pilot(f"{pilot_server}/reads.fq.gz", tmp_path / "x.fq.gz", 100)
    assert _records(out) == 100 and PilotHandler.requests == 3


def test_pilot_short_file_is_accepted_whole(pilot_server, tmp_path):
    out = _pilot(f"{pilot_server}/reads.fq.gz", tmp_path / "x.fq.gz", 5000)
    assert _records(out) == 1000


def test_pilot_truncated_stream_fails_with_reason(pilot_server, tmp_path):
    # the cut stream holds ~108 complete reads; asking for 500 makes the truncation land inside the request
    with pytest.raises(PipelineError) as e:
        _pilot(f"{pilot_server}/cut.fq.gz", tmp_path / "x.fq.gz", 500)
    assert "3 attempts" in str(e.value) and "ended unexpectedly" in (e.value.cause or "")
    assert not (tmp_path / "x.fq.gz").exists() and not (tmp_path / "x.fq.gz.part").exists()


def test_pilot_permanent_error_explained(pilot_server, tmp_path):
    with pytest.raises(PipelineError) as e:
        _pilot(f"{pilot_server}/missing.fq.gz", tmp_path / "x.fq.gz", 100)
    assert "HTTP error" in (e.value.cause or "")


def test_ena_failure_falls_back_to_sra_without_mixing_mates(tmp_path, monkeypatch):
    """F18b: a run ENA cannot deliver is fetched from NCBI SRA; a mate already fetched from ENA is set aside."""
    from rnaseq_pipeline import ena_manager, runner, sra_manager, workflow
    from rnaseq_pipeline import environment_manager as EM, system_check as SC
    from rnaseq_pipeline.project import Project
    p = Project.create("fb", tmp_path)
    for acc in ("SRR0000001", "SRR0000002"):
        p.add_sample(acc, source="ena", accession=acc, layout="PAIRED", r1=None, r2=None,
                     download={"files": [["https://x/1.gz", "m1", 10], ["https://x/2.gz", "m2", 10]]})
    p.state["read_type"] = "paired"
    p.save()
    fq = p.path("data", "fastq")
    (fq / "SRR0000002_1.fastq.gz").write_bytes(b"from ENA")      # half a pair from a failed ENA attempt

    def ena(run, out_dir, **kw):
        if run["run_accession"] == "SRR0000002":
            raise PipelineError("download of SRR0000002_2.fastq.gz failed: the network connection was interrupted")
        return [(Path(out_dir) / f"{run['run_accession']}_{i}.fastq.gz") for i in (1, 2)]

    def sra(acc, raw, fastq_dir, temp, threads, paired, **kw):
        outs = [Path(fastq_dir) / f"{acc}_{i}.fastq.gz" for i in (1, 2)]
        for o in outs:
            o.write_bytes(b"from SRA")
        return outs

    for i in (1, 2):
        (fq / f"SRR0000001_{i}.fastq.gz").write_bytes(b"ok")
    monkeypatch.setattr(ena_manager, "download_run", lambda run, out_dir, **kw: tuple(ena(run, out_dir, **kw)))
    monkeypatch.setattr(sra_manager, "download_run", lambda *a, **kw: tuple(sra(*a, **kw)))
    monkeypatch.setattr(runner, "which", lambda t: "/bin/true")
    cfg = p.config()
    ctx = workflow.Context(p, cfg, EM.Environments(cfg), SC.collect(p.root))
    workflow.DataStage().execute(ctx)
    rec = p.samples["SRR0000002"]
    assert rec["source"] == "sra" and rec["status"] == "OK"
    assert (fq / "SRR0000002_1.fastq.gz").read_bytes() == b"from SRA"
    assert (fq / "SRR0000002_2.fastq.gz").read_bytes() == b"from SRA"
    assert (fq / "SRR0000002_1.fastq.gz.other_archive").read_bytes() == b"from ENA"


# ---------------------------------------------------------------- P2 finding H11: stalled / flaky connections
def test_stalled_transfer_is_aborted_and_resumed(server, tmp_path):
    """Observed in the real-data run: the connection stopped delivering data and curl waited indefinitely."""
    import time
    Handler.stalls_left = 1
    t0 = time.time()
    out = net.download(f"{server}/stall.bin", tmp_path / "s.bin", stage="t", log_file=tmp_path / "dl.log",
                       retries=2, expected_md5=MD5, stall_timeout=2)
    assert out.read_bytes() == DATA
    # aborting (2 s) + curl's retry delay (5 s) + the rest takes ~13 s; waiting out the stall would take >= 120 s
    assert time.time() - t0 < 60, "the stalled transfer was not aborted"


def test_drops_with_progress_do_not_exhaust_retries(server, tmp_path):
    """A flaky connection that keeps delivering data must finish, however many times it drops."""
    out = net.download(f"{server}/flaky.bin", tmp_path / "f.bin", stage="t", log_file=tmp_path / "dl.log",
                       retries=2, expected_md5=MD5, stall_timeout=5)
    assert out.read_bytes() == DATA          # 2 MB in 300 kB pieces: ~7 connections with retries=2


def test_no_progress_still_gives_up(server, tmp_path):
    with pytest.raises(PipelineError) as e:
        net.download(f"{server}/dead.bin", tmp_path / "d.bin", stage="t", log_file=tmp_path / "dl.log",
                     retries=2, stall_timeout=2)
    assert "no data arrived" in str(e.value) or "failed" in str(e.value)
