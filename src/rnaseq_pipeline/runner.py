"""Safe external command execution.

* Commands are argument lists; nothing is ever passed through a shell.
* Pipelines (a | b | c) are built from Popen objects; every exit code is
  checked (equivalent of `set -o pipefail`).
* stderr of every process goes to a per-step log file; the terminal only shows
  a concise status line (unless verbose).
* Every command is recorded in logs/command_history.log + commands.jsonl.
* Ctrl-C terminates the whole process group and raises KeyboardInterrupt.
"""
import os
import shlex
import signal
import subprocess
import time
from pathlib import Path

from . import PipelineError, logger, ui

DRY_RUN = False
_tool_env = {"path_prefix": [], "extra": {}}


def configure_tool_path(bin_dirs, extra_env=None):
    """Prepend conda env bin dirs to PATH for all child processes."""
    _tool_env["path_prefix"] = [str(b) for b in bin_dirs if b]
    _tool_env["extra"] = dict(extra_env or {})


def child_env():
    env = os.environ.copy()
    if _tool_env["path_prefix"]:
        env["PATH"] = os.pathsep.join(_tool_env["path_prefix"] + [env.get("PATH", "")])
    env.update(_tool_env["extra"])
    env["PYTHONNOUSERSITE"] = "1"  # keep ~/.local packages out of tool envs
    env.setdefault("LC_ALL", "C.UTF-8")
    return env


def which(tool):
    for d in child_env()["PATH"].split(os.pathsep):
        p = Path(d) / tool
        if p.is_file() and os.access(p, os.X_OK):
            return p
    return None


def fmt(cmd):
    return " ".join(shlex.quote(str(c)) for c in cmd)


def fmt_pipeline(cmds, stdout_file=None):
    s = " | ".join(fmt(c) for c in cmds)
    return s + (f" > {shlex.quote(str(stdout_file))}" if stdout_file else "")


def _check_args(cmd):
    if not cmd or not all(isinstance(c, (str, Path, int, float)) for c in cmd):
        raise PipelineError(f"internal error: malformed command {cmd!r}")
    for c in cmd:
        if isinstance(c, str) and c == "":
            raise PipelineError(f"internal error: empty argument in command {fmt(cmd)}",
                                cause="an unset variable produced an empty argument")


def run(cmd, *, stage, sample=None, log_file=None, stdout_file=None, cwd=None,
        capture=False, check=True, description=None, expected_codes=(0,), timeout=None):
    """Run one command. Returns CompletedProcess-like object (returncode, stdout)."""
    return run_pipeline([cmd], stage=stage, sample=sample, log_file=log_file, stdout_file=stdout_file,
                        cwd=cwd, capture=capture, check=check, description=description,
                        expected_codes=expected_codes, timeout=timeout)


class Result:
    def __init__(self, codes, stdout, duration, command_str):
        self.returncodes = codes
        self.returncode = next((c for c in codes if c != 0), 0)
        self.stdout = stdout
        self.duration = duration
        self.command_str = command_str


def run_pipeline(cmds, *, stage, sample=None, log_file=None, stdout_file=None, cwd=None,
                 capture=False, check=True, description=None, expected_codes=(0,), timeout=None):
    for c in cmds:
        _check_args(c)
    cmds = [[str(x) for x in c] for c in cmds]
    command_str = fmt_pipeline(cmds, stdout_file)
    if description:
        ui.running(description)
    ui.debug(f"$ {command_str}")

    if DRY_RUN:
        ui.status("DRY-RUN", command_str)
        return Result([0] * len(cmds), "", 0.0, command_str)

    log_fh = None
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        log_fh = open(log_file, "ab")
        log_fh.write(f"\n### {time.strftime('%Y-%m-%d %H:%M:%S')}  $ {command_str}\n".encode())
        log_fh.flush()
    out_fh = open(stdout_file, "wb") if stdout_file else None

    procs = []
    start = time.time()
    status = "ok"
    stdout_data = ""
    try:
        prev = None
        for i, c in enumerate(cmds):
            last = i == len(cmds) - 1
            if last:
                stdout = out_fh if out_fh else (subprocess.PIPE if capture else (log_fh or subprocess.DEVNULL))
            else:
                stdout = subprocess.PIPE
            try:
                p = subprocess.Popen(
                    c, stdin=prev.stdout if prev else subprocess.DEVNULL, stdout=stdout,
                    stderr=log_fh or subprocess.DEVNULL, cwd=cwd, env=child_env(),
                    start_new_session=True)
            except FileNotFoundError:
                raise PipelineError(f"program not found: {c[0]}", stage=stage, sample=sample,
                                    cause="required tool is not installed or not on PATH",
                                    remedy="use main menu 4 (Manage Dependencies) to install it") from None
            if prev is not None:
                prev.stdout.close()  # let upstream receive SIGPIPE if downstream exits
            procs.append(p)
            prev = p
        if capture and not out_fh:
            data, _ = procs[-1].communicate(timeout=timeout)
            stdout_data = data.decode(errors="replace") if data else ""
        for p in procs:
            p.wait(timeout=timeout)
    except KeyboardInterrupt:
        status = "interrupted"
        _terminate(procs)
        raise
    except subprocess.TimeoutExpired:
        status = "timeout"
        _terminate(procs)
        raise PipelineError(f"command timed out after {timeout}s: {command_str}", stage=stage, sample=sample) from None
    finally:
        codes = [p.poll() for p in procs]
        duration = round(time.time() - start, 2)
        if out_fh:
            out_fh.close()
        if log_fh:
            log_fh.write(f"### exit codes {codes}  duration {duration}s\n".encode())
            log_fh.close()
        if status == "ok" and any(c not in expected_codes for c in codes):
            status = "failed"
        logger.record_command({
            "stage": stage, "sample": sample, "command": cmds, "command_str": command_str,
            "exit_codes": codes, "duration_s": duration, "status": status,
            "log_file": str(log_file) if log_file else None, "cwd": str(cwd) if cwd else None,
        })

    result = Result(codes, stdout_data, duration, command_str)
    if check and any(c not in expected_codes for c in codes):
        tail = tail_file(log_file) if log_file else ""
        bad = [i for i, c in enumerate(codes) if c not in expected_codes][0]
        meaning = describe_exit(codes[bad])
        raise PipelineError(
            f"command failed (exit codes {codes}): {cmds[bad][0]}" + (f" — {meaning}" if meaning else ""),
            stage=stage, sample=sample,
            cause=(tail or "see log file") + (f"\n  Full log: {log_file}" if log_file else ""),
            remedy=("run 'rnaseq-pipeline --check': it names a build of the program that runs on this CPU"
                    if "SIGILL" in meaning else
                    "inspect the log above; fix the input/tool problem and resume the project"))
    return result


def _terminate(procs):
    for p in procs:
        if p.poll() is None:
            try:
                os.killpg(p.pid, signal.SIGTERM)
            except (ProcessLookupError, PermissionError):
                pass
    deadline = time.time() + 10
    for p in procs:
        try:
            p.wait(timeout=max(0.1, deadline - time.time()))
        except subprocess.TimeoutExpired:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass


def tail_file(path, n=12):
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 8000))
            lines = f.read().decode(errors="replace").splitlines()
        lines = [l for l in lines if l.strip() and not l.startswith("###")]
        return "\n".join("    " + l for l in lines[-n:])
    except OSError:
        return ""


def signal_of(code):
    """The signal that killed a process, from a Popen return code (-N) or a shell-style one (128 + N); else None."""
    if code is None:
        return None
    return -code if code < 0 else (code - 128 if 128 < code < 160 else None)


def describe_exit(code):
    """Plain-language meaning of an exit code that says the program was killed by a signal."""
    sig = signal_of(code)
    if sig == signal.SIGILL:
        from . import system_check
        missing = system_check.missing_x86_64_v3()
        return ("killed by SIGILL (illegal CPU instruction): this build of the program needs processor features "
                + (f"this CPU lacks ({', '.join(missing)})" if missing else "this CPU does not have"))
    if sig == signal.SIGKILL:
        return "killed by SIGKILL (often the system running out of memory)"
    if sig == signal.SIGSEGV:
        return "crashed (SIGSEGV, segmentation fault)"
    return ""


def tool_run(cmd, timeout=60):
    """Like tool_output, but returns (exit code, output); (None, None) if the program cannot be started."""
    try:
        p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             stdin=subprocess.DEVNULL, env=child_env(), start_new_session=True)
    except OSError:
        return None, None
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        _terminate([p])
        return None, None
    except BaseException:
        # Ctrl-C / SIGTERM: many tools are wrapper scripts (hisat2, fastqc) whose real program is a grandchild,
        # so the whole process group is stopped, as for pipeline commands
        _terminate([p])
        raise
    return p.returncode, (out + err).decode(errors="replace")


def tool_output(cmd, timeout=60):
    """Run a quick informational command (e.g. --version); returns combined output or None."""
    return tool_run(cmd, timeout)[1]
