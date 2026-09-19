"""Project log files: pipeline.log, pipeline_errors.log, command_history.log, commands.jsonl."""
import json
import logging
import threading
from datetime import datetime
from pathlib import Path

log = logging.getLogger("rnaseq")
_lock = threading.Lock()
_state = {"logs_dir": None, "handlers": []}

FMT = "%(asctime)s %(levelname)-7s %(message)s"


def init_base_logging(verbose=False):
    log.setLevel(logging.DEBUG)
    log.propagate = False


def attach_project(logs_dir):
    """Route all log records to files inside <project>/logs."""
    logs_dir = Path(logs_dir)
    logs_dir.mkdir(parents=True, exist_ok=True)
    detach_project()
    main = logging.FileHandler(logs_dir / "pipeline.log", encoding="utf-8")
    main.setLevel(logging.DEBUG)
    main.setFormatter(logging.Formatter(FMT))
    errs = logging.FileHandler(logs_dir / "pipeline_errors.log", encoding="utf-8")
    errs.setLevel(logging.WARNING)
    errs.setFormatter(logging.Formatter(FMT))
    for h in (main, errs):
        log.addHandler(h)
    _state["handlers"] = [main, errs]
    _state["logs_dir"] = logs_dir


def detach_project():
    for h in _state["handlers"]:
        log.removeHandler(h)
        h.close()
    _state["handlers"] = []
    _state["logs_dir"] = None


def logs_dir():
    return _state["logs_dir"]


def record_command(entry):
    """Append a structured command record (dict) to command_history.log and commands.jsonl."""
    d = _state["logs_dir"]
    if d is None:
        return
    entry = dict(entry)
    entry.setdefault("timestamp", datetime.now().isoformat(timespec="seconds"))
    with _lock:
        with open(d / "commands.jsonl", "a", encoding="utf-8") as f:
            f.write(json.dumps(entry, default=str) + "\n")
        with open(d / "command_history.log", "a", encoding="utf-8") as f:
            f.write(
                f"[{entry['timestamp']}] stage={entry.get('stage')} sample={entry.get('sample')} "
                f"status={entry.get('status')} exit={entry.get('exit_codes')} "
                f"duration={entry.get('duration_s')}s\n"
                f"  $ {entry.get('command_str')}\n"
                + (f"  log: {entry.get('log_file')}\n" if entry.get("log_file") else "")
            )
