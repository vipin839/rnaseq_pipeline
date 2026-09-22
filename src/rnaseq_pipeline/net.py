"""HTTP helpers: small API requests (urllib) and resumable, verified file downloads (curl)."""
import hashlib
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from . import PipelineError, __version__, runner, ui
from .secrets import scrub

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
    raise PipelineError(f"network request failed: {scrub(url)}", cause=scrub(str(last)),
                        remedy="check the internet connection / accession and try again")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# curl exit codes -> plain-language explanation, and whether trying again can help.
CURL_ERRORS = {
    6: ("could not resolve the server name (no internet connection or DNS problem)", True),
    7: ("could not connect to the server (down, blocked by a firewall, or no internet)", True),
    18: ("the transfer ended before the whole file arrived", True),
    22: ("the server returned an HTTP error (for example 404 not found or 403 forbidden)", False),
    23: ("could not write the file (disk full or no write permission)", False),
    28: ("the connection timed out", True),
    33: ("the server does not support resuming a partial download", True),
    35: ("the secure (TLS/SSL) connection could not be established", True),
    52: ("the server closed the connection without sending data", True),
    56: ("the network connection was interrupted", True),
}


def download(url, dest, *, stage, sample=None, expected_md5=None, expected_bytes=None, log_file=None, retries=3):
    """Download url to dest safely.

    * data go to dest.part; an interrupted transfer is resumed (or restarted if the server cannot resume)
    * the result must be non-empty, have the expected size (if known) and MD5 (if known)
    * only then is it renamed to dest; an existing dest is re-verified and never overwritten
    """
    dest = Path(dest)
    if dest.exists():
        problem = _verify(dest, expected_md5, expected_bytes)
        if problem:
            raise PipelineError(f"existing file {dest.name} is not valid: {problem}", stage=stage, sample=sample,
                                remedy=f"move {dest} aside (it is never overwritten automatically) and resume")
        ui.skipped(f"already downloaded and verified: {dest.name}")
        return dest
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    resume = True
    for attempt in range(1, retries + 1):
        cmd = ["curl", "-fL", "--retry", str(retries), "--retry-delay", "5", "--connect-timeout", "30",
               *(["-C", "-"] if resume else []), "-o", str(part), url]
        res = runner.run(cmd, stage=stage, sample=sample, log_file=log_file, check=False,
                         description=f"Downloading {dest.name}" + (f" (attempt {attempt})" if attempt > 1 else ""))
        if res.returncode == 0:
            break
        why, transient = CURL_ERRORS.get(res.returncode, (f"curl exit code {res.returncode}", True))
        if res.returncode == 33 or (part.exists() and res.returncode == 22 and part.stat().st_size and resume):
            resume = False  # server cannot resume (or the partial file is stale): start this file again cleanly
            part.unlink(missing_ok=True)
            ui.warn(f"{dest.name}: {why}; restarting the download from the beginning")
            continue
        if not transient or attempt == retries:
            raise PipelineError(f"download of {dest.name} failed: {why}", stage=stage, sample=sample,
                                cause=f"URL: {scrub(url)}" + (f"\n  Log: {log_file}" if log_file else ""),
                                remedy="check the internet connection and that the data are public, then resume "
                                       "(completed parts are kept)" if transient else
                                       "the file cannot be fetched from this address; check the accession/URL")
        ui.warn(f"{dest.name}: {why}; retrying ({attempt}/{retries}), keeping what was already downloaded")
    if runner.DRY_RUN:
        return dest
    problem = _verify(part, expected_md5, expected_bytes)
    if problem:
        bad = part.with_name(part.name + ".rejected")
        os.replace(part, bad)
        raise PipelineError(f"downloaded {dest.name} failed verification: {problem}", stage=stage, sample=sample,
                            cause=f"the file was kept for inspection as {bad.name}",
                            remedy="resume the project to download it again")
    if expected_md5:
        ui.ok(f"MD5 verified: {dest.name}")
    os.replace(part, dest)
    return dest


def _verify(path, expected_md5=None, expected_bytes=None):
    """Return a problem description, or None if the file passes every check available."""
    size = Path(path).stat().st_size
    if size == 0:
        return "the file is empty (0 bytes)"
    if expected_bytes and size != int(expected_bytes):
        return f"size is {size:,} bytes, the archive lists {int(expected_bytes):,}"
    if expected_md5:
        got = md5(path)
        if got != expected_md5:
            return f"MD5 checksum {got} does not match the published {expected_md5}"
    return None
