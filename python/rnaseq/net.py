"""HTTP helpers: small API requests (urllib) and resumable, verified file downloads (curl)."""
import hashlib
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from . import PipelineError, __version__, runner, ui

UA = f"rnaseq-pipeline/{__version__} (research; python-urllib)"


def get_text(url, params=None, retries=3, timeout=60):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (400, 404):
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = e
        time.sleep(min(2 ** attempt, 15))
    raise PipelineError(f"network request failed: {url}", cause=str(last),
                        remedy="check the internet connection / accession and try again")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest, *, stage, sample=None, expected_md5=None, log_file=None, retries=3):
    """Resumable download to dest.part, verify MD5 if known, then atomically rename.
    An existing, verified dest is never overwritten."""
    dest = Path(dest)
    if dest.exists():
        if expected_md5 and md5(dest) != expected_md5:
            raise PipelineError(f"existing file has wrong checksum: {dest}", stage=stage, sample=sample,
                                remedy="move the file aside and download again")
        ui.skipped(f"already downloaded: {dest.name}")
        return dest
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["curl", "-fL", "--retry", str(retries), "--retry-delay", "5", "--connect-timeout", "30",
           "-C", "-", "-o", str(part), url]
    for attempt in range(1, retries + 1):
        try:
            runner.run(cmd, stage=stage, sample=sample, log_file=log_file,
                       description=f"Downloading {dest.name}" + (f" (attempt {attempt})" if attempt > 1 else ""))
            break
        except PipelineError:
            if attempt == retries:
                raise
            ui.warn(f"download interrupted; resuming ({attempt}/{retries})")
    if runner.DRY_RUN:
        return dest
    if expected_md5:
        got = md5(part)
        if got != expected_md5:
            bad = part.with_name(part.name + ".badmd5")
            os.replace(part, bad)
            raise PipelineError(f"checksum mismatch for {dest.name}", stage=stage, sample=sample,
                                cause=f"expected MD5 {expected_md5}, got {got} (file kept as {bad.name})",
                                remedy="re-download (choose 'retry' when prompted)")
        ui.ok(f"MD5 verified: {dest.name}")
    os.replace(part, dest)
    return dest
