"""Bulk RNA-seq pipeline — Python orchestration package."""
import os
from pathlib import Path

# The ONLY place the version is defined (pyproject.toml reads it from here).
__version__ = "1.1.0"

# Data shipped inside the package: R scripts, default configuration, reference catalog, conda env files.
PACKAGE_DIR = Path(__file__).resolve().parent

# Per-user state (install/update logs, recently opened projects). Never inside the package directory,
# which may be read-only after installation.
USER_STATE_DIR = Path(os.environ.get("RNASEQ_PIPELINE_HOME", Path.home() / ".rnaseq_pipeline"))


class PipelineError(Exception):
    """A handled pipeline failure with a user-facing explanation."""

    def __init__(self, message, *, cause=None, remedy=None, sample=None, stage=None):
        super().__init__(message)
        self.cause = cause
        self.remedy = remedy
        self.sample = sample
        self.stage = stage


class UserAbort(Exception):
    """The user chose to cancel, or stdin closed."""


class Terminated(KeyboardInterrupt):
    """SIGTERM/SIGHUP received (kill, closed terminal, dropped SSH). Handled like Ctrl-C, then the program exits."""

    def __init__(self, signame):
        super().__init__(signame)
        self.signame = signame


def install_signal_handlers():
    """Turn SIGTERM/SIGHUP into Terminated so child process groups are stopped and no partial output is used."""
    import signal

    def handler(signum, _frame):
        raise Terminated(signal.Signals(signum).name)

    for sig in (signal.SIGTERM, signal.SIGHUP):
        signal.signal(sig, handler)
