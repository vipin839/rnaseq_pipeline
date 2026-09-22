"""F3/F4/F5 regression: downloads are never accepted just because a file exists."""
import hashlib
import http.server
import os
import threading

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
