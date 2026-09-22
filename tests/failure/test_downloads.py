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

    def do_GET(self):
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
