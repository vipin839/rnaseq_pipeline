"""Stage checkpoints.

A checkpoint records, for a completed stage: the outputs (size + mtime, and a
SHA-256 for small files), the parameters used, and the fingerprints of the
upstream checkpoints it depended on. A checkpoint is trusted only if:

  1. every recorded output still exists with the same size/mtime (and hash),
  2. every upstream checkpoint is still valid and unchanged, and
  3. the stage's own re-validation function (e.g. BAM quickcheck) passes.
"""
import hashlib
import json
import os
from pathlib import Path

from . import __version__
from .project import now

HASH_LIMIT = 50 * 1024 * 1024  # hash files up to 50 MB; larger ones use size+mtime + stage validator


def sha256(path, limit=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


QUICK_BLOCK = 1 << 20


def fingerprint(path):
    """Identity of a file for the manifest: full SHA-256 up to HASH_LIMIT, otherwise a quick fingerprint
    (SHA-256 of size + first and last MiB) that detects replacement/truncation without reading 100 GB."""
    p = Path(path)
    size = p.stat().st_size
    if size <= HASH_LIMIT:
        return {"size": size, "sha256": sha256(p)}
    h = hashlib.sha256(str(size).encode())
    with open(p, "rb") as f:
        h.update(f.read(QUICK_BLOCK))
        f.seek(max(0, size - QUICK_BLOCK))
        h.update(f.read(QUICK_BLOCK))
    return {"size": size, "quick_sha256": h.hexdigest(), "note": "size + first/last MiB"}


def describe_file(project, path):
    p = Path(path)
    st = p.stat()
    rec = {"path": project.rel(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    if st.st_size <= HASH_LIMIT:
        rec["sha256"] = sha256(p)
    return rec


class Checkpoints:
    def __init__(self, project):
        self.project = project
        self.dir = project.path("checkpoints")

    def file(self, stage):
        return self.dir / f"{stage}.json"

    def fingerprint(self, stage):
        """Content fingerprint: outputs (hash, or size+mtime for large files) + params, never timestamps,
        so re-running a stage that reproduces identical outputs does not invalidate downstream stages."""
        data = self.read(stage)
        if data is None:
            return None
        outs = [(o["path"], o["size"], o.get("sha256") or o["mtime_ns"]) for o in data["outputs"]]
        blob = json.dumps({"outputs": outs, "params": data.get("params")}, sort_keys=True, default=str)
        return hashlib.sha256(blob.encode()).hexdigest()

    def write(self, stage, outputs, params=None, depends_on=(), summary=None):
        outs = []
        for o in outputs:
            o = Path(o)
            if o.is_dir():
                for sub in sorted(o.rglob("*")):
                    if sub.is_file():
                        outs.append(describe_file(self.project, sub))
            elif o.exists():
                outs.append(describe_file(self.project, o))
            else:
                raise FileNotFoundError(f"checkpoint output missing: {o}")
        data = {
            "stage": stage, "completed_at": now(), "pipeline_version": __version__,
            "outputs": outs, "params": params or {}, "summary": summary or {},
            "depends_on": {d: self.fingerprint(d) for d in depends_on},
        }
        tmp = self.file(stage).with_suffix(".json.tmp")
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, default=str)
        os.replace(tmp, self.file(stage))
        return data

    def read(self, stage):
        f = self.file(stage)
        if not f.exists():
            return None
        try:
            with open(f, encoding="utf-8") as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            return None

    def exists(self, stage):
        return self.file(stage).exists()

    def verify_files(self, stage, deep=False):
        """Check recorded outputs. Returns list of problems."""
        data = self.read(stage)
        if data is None:
            return ["checkpoint missing or unreadable"]
        problems = []
        for o in data["outputs"]:
            p = self.project.abs(o["path"])
            if not p.exists():
                problems.append(f"missing output: {o['path']}")
                continue
            st = p.stat()
            if st.st_size != o["size"]:
                problems.append(f"size changed: {o['path']}")
            elif "sha256" in o and (deep or st.st_mtime_ns != o["mtime_ns"]):
                if sha256(p) != o["sha256"]:
                    problems.append(f"content changed: {o['path']}")
            elif "sha256" not in o and st.st_mtime_ns != o["mtime_ns"]:
                problems.append(f"modified since checkpoint: {o['path']}")
        for dep, fp in data.get("depends_on", {}).items():
            if self.fingerprint(dep) != fp:
                problems.append(f"upstream stage '{dep}' changed or was re-run")
        return problems

    def invalidate(self, stage, reason):
        """Move a checkpoint aside (never deletes outputs)."""
        f = self.file(stage)
        if f.exists():
            dest = f.with_suffix(f".invalid.{now().replace(':', '')}.json")
            os.replace(f, dest)
            self.project.log_event(f"checkpoint {stage} invalidated: {reason}")
