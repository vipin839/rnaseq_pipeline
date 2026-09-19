"""Streaming FASTQ validation (constant memory, any file size).

Checks every record: 4-line structure, '@' header, '+' separator, sequence
alphabet, quality alphabet, len(quality) == len(sequence), truncation, gzip
integrity (decompressor exit status), and — for pairs — identical read counts
and matching mate identifiers, read in lockstep.
"""
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPORT_COLUMNS = ["sample", "file", "mate", "status", "reads", "bases", "min_length", "max_length",
                  "mean_length", "min_quality_char", "max_quality_char", "problem", "record", "reason",
                  "recommended_action"]


class FastqError(Exception):
    def __init__(self, problem, record, reason, action, file=None):
        super().__init__(problem)
        self.problem, self.record, self.reason, self.action, self.file = problem, record, reason, action, file


def _open(path):
    """Return (binary stream, process or None). gzip is decoded by an external process so that
    corruption/truncation surfaces as a non-zero exit code."""
    path = Path(path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        dec = shutil.which("pigz") or shutil.which("gzip")
        if dec is None:
            import gzip
            return gzip.open(path, "rb"), None
        p = subprocess.Popen([dec, "-dc", str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             bufsize=1 << 20)
        return p.stdout, p
    if path.suffix in (".bz2", ".xz", ".zst"):
        raise FastqError("unsupported compression", 0, f"{path.suffix} is not supported",
                         "recompress with gzip (bgzip/pigz)")
    return open(path, "rb", buffering=1 << 20), None


def _mate_id(header):
    tok = header[1:].split(None, 1)[0] if len(header) > 1 else b""
    if tok.endswith((b"/1", b"/2")):
        tok = tok[:-2]
    return tok


class _Reader:
    def __init__(self, path, allowed, qmin, qmax):
        self.path = Path(path)
        self.stream, self.proc = _open(path)
        self.allowed = allowed
        self.qual_ok = bytes(range(qmin, qmax + 1))
        self.n = 0
        self.bases = 0
        self.minlen = None
        self.maxlen = 0
        self.qlo = 255
        self.qhi = 0

    def err(self, problem, reason, action):
        return FastqError(problem, self.n + 1, reason, action, file=str(self.path))

    def next(self):
        """Return mate id of the next record, or None at clean EOF."""
        rl = self.stream.readline
        h = rl()
        if not h:
            return None
        s, p, q = rl(), rl(), rl()
        if not q:
            raise self.err("truncated file", "file ends in the middle of a FASTQ record",
                           "re-download or re-copy the file; it is incomplete")
        if not q.endswith(b"\n"):
            # last line of file without newline is acceptable only if complete
            pass
        h, s, p, q = h.rstrip(b"\r\n"), s.rstrip(b"\r\n"), p.rstrip(b"\r\n"), q.rstrip(b"\r\n")
        if not h.startswith(b"@") or len(h) < 2:
            raise self.err("invalid header", f"line {self.n * 4 + 1} does not start with '@': {h[:60]!r}",
                           "file is not FASTQ, or records are misaligned (corrupted)")
        if not p.startswith(b"+"):
            raise self.err("invalid separator", f"line {self.n * 4 + 3} does not start with '+': {p[:60]!r}",
                           "records are misaligned: a line is missing or extra (corrupted file)")
        if len(p) > 1 and p[1:] != h[1:]:
            raise self.err("separator mismatch", "'+' line repeats a different identifier than the header",
                           "file is corrupted")
        if len(s) != len(q):
            raise self.err("length mismatch", f"sequence length {len(s)} != quality length {len(q)}",
                           "file is corrupted or truncated")
        if s.translate(None, self.allowed):
            bad = sorted(set(s.translate(None, self.allowed)))[:5]
            raise self.err("invalid sequence character", f"unexpected characters {bytes(bad)!r} in sequence",
                           "file is corrupted or not nucleotide FASTQ")
        if q.translate(None, self.qual_ok):
            bad = sorted(set(q.translate(None, self.qual_ok)))[:5]
            raise self.err("invalid quality character", f"quality characters {bytes(bad)!r} outside the "
                           "allowed Phred range", "file is corrupted or uses an unsupported quality encoding")
        L = len(s)
        self.n += 1
        self.bases += L
        if self.minlen is None or L < self.minlen:
            self.minlen = L
        if L > self.maxlen:
            self.maxlen = L
        if q:
            lo, hi = min(q), max(q)
            if lo < self.qlo:
                self.qlo = lo
            if hi > self.qhi:
                self.qhi = hi
        return _mate_id(h)

    def close(self):
        """Close and check decompressor status."""
        rc = 0
        err = b""
        if self.proc:
            try:
                self.stream.read()  # drain
            except Exception:
                pass
            err = self.proc.stderr.read() if self.proc.stderr else b""
            rc = self.proc.wait()
        self.stream.close()
        if rc != 0:
            raise FastqError("gzip integrity failure", self.n + 1,
                             f"decompression failed ({err.decode(errors='replace').strip()[:200]})",
                             "the compressed file is corrupted or truncated; re-download it", file=str(self.path))

    def abort(self):
        try:
            if self.proc and self.proc.poll() is None:
                self.proc.kill()
                self.proc.wait()
            self.stream.close()
        except Exception:
            pass

    def stats(self):
        return {"reads": self.n, "bases": self.bases, "min_length": self.minlen or 0,
                "max_length": self.maxlen,
                "mean_length": round(self.bases / self.n, 1) if self.n else 0,
                "min_quality_char": chr(self.qlo) if self.n else "",
                "max_quality_char": chr(self.qhi) if self.n else ""}


def validate_sample(sample, r1, r2=None, allowed="ACGTNacgtn.", offset=33, min_reads=1):
    """Validate a single-end file or a pair in lockstep. Returns list of report rows (one per file)."""
    allowed_b = allowed.encode()
    qmin, qmax = offset, 126
    files = [("R1", r1)] + ([("R2", r2)] if r2 else [])
    rows = {m: {"sample": sample, "file": str(f), "mate": m, "status": "PASS", "problem": "", "record": "",
                "reason": "", "recommended_action": ""} for m, f in files}
    readers = {}
    try:
        for m, f in files:
            fp = Path(f)
            if not fp.exists():
                raise FastqError("file not found", 0, f"{fp} does not exist", "check the path / re-import", str(fp))
            if not os.access(fp, os.R_OK):
                raise FastqError("file not readable", 0, "permission denied", "fix file permissions", str(fp))
            if fp.stat().st_size == 0:
                raise FastqError("empty file", 0, "file size is 0 bytes", "re-download / re-copy", str(fp))
            readers[m] = _Reader(fp, allowed_b, qmin, qmax)
        a = readers["R1"]
        b = readers.get("R2")
        while True:
            ida = a.next()
            if b is None:
                if ida is None:
                    break
                continue
            idb = b.next()
            if ida is None and idb is None:
                break
            if ida is None or idb is None:
                short = a if ida is None else b
                raise FastqError("paired read count mismatch",
                                 short.n + 1,
                                 f"R1 and R2 have different numbers of reads ({'R1' if ida is None else 'R2'} ended "
                                 f"after {short.n} reads)",
                                 "files are truncated or not mates; re-download both files", str(short.path))
            if ida != idb:
                raise FastqError("mate ID mismatch", a.n,
                                 f"R1 id {ida[:60].decode(errors='replace')!r} != R2 id "
                                 f"{idb[:60].decode(errors='replace')!r}",
                                 "R1/R2 files are not mates or are out of order", str(b.path))
        for r in readers.values():
            r.close()
        for m, r in readers.items():
            rows[m].update(r.stats())
            if r.n < min_reads:
                raise FastqError("no reads", 0, "file contains no FASTQ records", "re-download / re-copy",
                                 str(r.path))
        # Phred+64 heuristic: no character below ';' but many above 'J'
        for m, r in readers.items():
            if r.n and r.qlo >= 64 and r.qhi > 75 and offset == 33:
                rows[m]["status"] = "WARNING"
                rows[m]["problem"] = "possible Phred+64 encoding"
                rows[m]["reason"] = f"lowest quality character is {chr(r.qlo)!r}"
                rows[m]["recommended_action"] = "confirm encoding; convert to Phred+33 before alignment"
    except FastqError as e:
        for r in readers.values():
            r.abort()
        for m, r in readers.items():
            rows[m].update(r.stats())
        target = next((m for m, f in files if str(f) == e.file), "R1")
        for m in rows:
            rows[m]["status"] = "FAILED"
        rows[target].update({"problem": e.problem, "record": e.record, "reason": e.reason,
                             "recommended_action": e.action})
        for m in rows:
            if m != target:
                rows[m].update({"problem": "mate failed", "reason": f"{target} failed validation",
                                "recommended_action": e.action})
    return list(rows.values())


def _job(args):
    return validate_sample(*args)


def validate_many(jobs, workers=1, progress=None):
    """jobs: list of (sample, r1, r2, allowed, offset). Returns rows in job order."""
    results = {}
    if workers <= 1 or len(jobs) == 1:
        for j in jobs:
            results[j[0]] = _job(j)
            if progress:
                progress(j[0], results[j[0]])
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_job, j): j[0] for j in jobs}
            for f in as_completed(futs):
                results[futs[f]] = f.result()
                if progress:
                    progress(futs[f], results[futs[f]])
    rows = []
    for j in jobs:
        rows.extend(results[j[0]])
    return rows


def write_report(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(REPORT_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in REPORT_COLUMNS) + "\n")
    return path
