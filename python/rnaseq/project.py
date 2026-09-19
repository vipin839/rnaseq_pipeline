"""Project directory layout and persistent project state (project.json)."""
import fcntl
import json
import os
import uuid
from datetime import datetime
from pathlib import Path

from . import PIPELINE_ROOT, PipelineError, __version__
from . import config as C
from . import validators as V

DEFAULT_PROJECTS_DIR = PIPELINE_ROOT / "projects"

LAYOUT = [
    "config",
    "data/raw", "data/fastq", "data/trimmed", "data/metadata",
    "reference/genome", "reference/annotation", "reference/transcriptome", "reference/index",
    "reference/checksums", "reference/logs",
    "qc/fastqc_raw", "qc/multiqc_raw", "qc/fastqc_trimmed", "qc/multiqc_trimmed", "qc/assessment",
    "alignment/sam", "alignment/bam", "alignment/index", "alignment/reports", "alignment/logs",
    "stringtie/abundance", "stringtie/merged",
    "featurecounts/logs",
    "counts",
    "results/deseq2", "results/upregulated", "results/downregulated", "results/plots", "results/tables",
    "logs", "reports", "checkpoints", "temp", "pipeline_manifest",
]


def now():
    return datetime.now().isoformat(timespec="seconds")


class Project:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.state_file = self.root / "project.json"
        self.state = {}
        self._lock_fh = None

    # ---------- creation / loading ----------
    @classmethod
    def create(cls, name, parent=None, base_config=None):
        name = V.project_name(name)
        parent = Path(parent or DEFAULT_PROJECTS_DIR).expanduser().resolve()
        root = parent / f"RNAseq_{name}" if not name.startswith("RNAseq_") else parent / name
        if root.exists() and any(root.iterdir()):
            raise PipelineError(f"project directory already exists and is not empty: {root}",
                                remedy="choose another name or use 'Resume Existing Project'")
        p = cls(root)
        for d in LAYOUT:
            (root / d).mkdir(parents=True, exist_ok=True)
        p.state = {
            "project_id": str(uuid.uuid4()), "name": name, "created": now(),
            "pipeline_version": __version__, "data_source": None, "accessions": [],
            "samples": {}, "read_type": None, "reference": None, "strandedness": None,
            "qc_decision": None, "design": None, "selected_samples": None, "history": [],
        }
        C.save_yaml(base_config or C.load_yaml(C.DEFAULT_CONFIG), root / "config" / "project_config.yaml")
        p.save()
        return p

    @classmethod
    def open(cls, root):
        p = cls(root)
        if not p.state_file.exists():
            raise PipelineError(f"not a pipeline project (no project.json): {p.root}")
        with open(p.state_file, encoding="utf-8") as f:
            p.state = json.load(f)
        for d in LAYOUT:  # repair missing dirs silently (never deletes)
            (p.root / d).mkdir(parents=True, exist_ok=True)
        return p

    @staticmethod
    def list_projects(parent=None):
        parent = Path(parent or DEFAULT_PROJECTS_DIR).expanduser()
        if not parent.is_dir():
            return []
        return sorted(d for d in parent.iterdir() if (d / "project.json").exists())

    def save(self):
        tmp = self.state_file.with_suffix(".json.tmp")
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(self.state, f, indent=2, default=str)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, self.state_file)

    # ---------- locking (one pipeline process per project) ----------
    def lock(self):
        self._lock_fh = open(self.root / ".lock", "w")
        try:
            fcntl.flock(self._lock_fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise PipelineError("this project is already open in another pipeline process",
                                remedy="close the other session first")
        self._lock_fh.write(str(os.getpid()))
        self._lock_fh.flush()

    def unlock(self):
        if self._lock_fh:
            fcntl.flock(self._lock_fh, fcntl.LOCK_UN)
            self._lock_fh.close()
            self._lock_fh = None

    # ---------- paths ----------
    def path(self, *parts):
        return self.root.joinpath(*parts)

    @property
    def name(self):
        return self.state["name"]

    @property
    def config_path(self):
        return self.path("config", "project_config.yaml")

    def config(self):
        return C.load(self.config_path)

    def save_config(self, cfg):
        C.save_yaml(cfg, self.config_path)

    def rel(self, p):
        p = Path(p)
        try:
            return str(p.resolve().relative_to(self.root))
        except ValueError:
            return str(p)

    def abs(self, p):
        p = Path(p)
        return p if p.is_absolute() else self.root / p

    # ---------- samples ----------
    @property
    def samples(self):
        return self.state["samples"]

    def add_sample(self, sid, **fields):
        sid = V.sample_name(sid)
        if sid in self.samples:
            raise PipelineError(f"duplicate sample name: {sid}")
        rec = {"id": sid, "status": "OK", "fail_reason": None, "metadata": {}}
        rec.update(fields)
        self.samples[sid] = rec
        return rec

    def active_samples(self):
        """Samples that are selected for processing and not FAILED/EXCLUDED, in stable order."""
        sel = self.state.get("selected_samples")
        out = []
        for sid, rec in self.samples.items():
            if sel is not None and sid not in sel:
                continue
            if rec.get("status") in ("FAILED", "EXCLUDED"):
                continue
            out.append(sid)
        return out

    def mark_failed(self, sid, reason, stage):
        rec = self.samples[sid]
        rec["status"] = "FAILED"
        rec["fail_reason"] = reason
        rec["failed_stage"] = stage
        self.save()

    def fastqs(self, sid, trimmed=None):
        """Return (r1, r2|None) to use for sample. trimmed=None -> use trimmed if trimming completed."""
        rec = self.samples[sid]
        use_trim = rec.get("trimmed") and (trimmed is None or trimmed)
        if trimmed is False:
            use_trim = False
        src = rec["trimmed"] if use_trim else rec
        r1 = self.abs(src["r1"])
        r2 = self.abs(src["r2"]) if src.get("r2") else None
        return r1, r2

    def is_paired(self):
        return self.state.get("read_type") == "paired"

    def log_event(self, event):
        self.state.setdefault("history", []).append({"time": now(), "event": event})
        self.save()
