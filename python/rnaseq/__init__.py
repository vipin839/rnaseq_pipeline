"""Bulk RNA-seq pipeline, Version 1 — Python orchestration package."""
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
__version__ = (PIPELINE_ROOT / "VERSION").read_text().strip()


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
