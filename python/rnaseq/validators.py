"""Validation of user-supplied names, paths, accessions and numeric parameters."""
import os
import re
from pathlib import Path

NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
SAMPLE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._-]{0,63}$")
LEVEL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.]{0,63}$")
COLUMN_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,31}$")

ACCESSION_PATTERNS = {
    "run": re.compile(r"^[SED]RR\d{6,}$"),
    "experiment": re.compile(r"^[SED]RX\d{6,}$"),
    "sample": re.compile(r"^([SED]RS\d{6,}|SAM[NED][A-Z]?\d+)$"),
    "study": re.compile(r"^([SED]RP\d{6,}|PRJ[NED][A-Z]\d+)$"),
    "geo_series": re.compile(r"^GSE\d{1,8}$"),
    "geo_sample": re.compile(r"^GSM\d{1,9}$"),
}

RESERVED_NAMES = {".", "..", "con", "nul", "prn", "aux"}


def project_name(value):
    v = (value or "").strip()
    if not v:
        raise ValueError("name is empty")
    if "/" in v or "\\" in v or ".." in v:
        raise ValueError("name must not contain path separators or '..'")
    if v.lower() in RESERVED_NAMES:
        raise ValueError("reserved name")
    if not NAME_RE.match(v):
        raise ValueError("use letters, digits, '.', '_' or '-' (max 64 chars, must start with a letter/digit)")
    return v


def sample_name(value):
    v = (value or "").strip()
    if not SAMPLE_RE.match(v):
        raise ValueError(
            f"invalid sample name {v!r}: start with a letter; use letters, digits, '.', '_' or '-' (max 64)")
    return v


def factor_level(value):
    v = (value or "").strip()
    if not LEVEL_RE.match(v):
        raise ValueError(
            f"invalid value {v!r}: start with a letter; letters, digits, '_' or '.' only "
            "(R/DESeq2 requires syntactically valid factor levels)")
    return v


def column_name(value):
    v = (value or "").strip()
    if not COLUMN_RE.match(v):
        raise ValueError(f"invalid column name {v!r}: letters, digits, '_' only, start with a letter")
    return v


def accession(value):
    v = (value or "").strip().upper()
    for kind, pat in ACCESSION_PATTERNS.items():
        if pat.match(v):
            return v, kind
    raise ValueError(f"{value!r} is not a recognised SRA/ENA/GEO accession "
                     "(e.g. SRR1234567, ERR…, SRX…, SRP…, PRJNA…, GSE12345)")


def accession_list(text):
    items = [t for t in re.split(r"[\s,;]+", text or "") if t]
    if not items:
        raise ValueError("no accessions given")
    out, seen = [], set()
    for it in items:
        acc, kind = accession(it)
        if acc not in seen:
            seen.add(acc)
            out.append((acc, kind))
    return out


def existing_file(value, must_be_readable=True):
    p = Path(os.path.expanduser(str(value).strip())).resolve()
    if not str(value).strip():
        raise ValueError("path is empty")
    if not p.exists():
        raise ValueError(f"file not found: {p}")
    if not p.is_file():
        raise ValueError(f"not a regular file: {p}")
    if must_be_readable and not os.access(p, os.R_OK):
        raise ValueError(f"file not readable: {p}")
    return p


def existing_dir(value):
    if not str(value).strip():
        raise ValueError("path is empty")
    p = Path(os.path.expanduser(str(value).strip())).resolve()
    if not p.is_dir():
        raise ValueError(f"directory not found: {p}")
    return p


def within(base, candidate):
    """Resolve candidate and ensure it stays inside base (no path traversal)."""
    base = Path(base).resolve()
    c = (base / candidate).resolve()
    if c != base and base not in c.parents:
        raise ValueError(f"path escapes project directory: {candidate}")
    return c


def positive_int(value, name="value", minimum=1, maximum=None):
    try:
        v = int(str(value).strip())
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be an integer")
    if v < minimum or (maximum is not None and v > maximum):
        rng = f">= {minimum}" + (f" and <= {maximum}" if maximum is not None else "")
        raise ValueError(f"{name} must be {rng}")
    return v


def number(value, name="value", minimum=None, maximum=None, exclusive_min=False, exclusive_max=False):
    try:
        v = float(str(value).strip())
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number")
    if v != v:  # NaN
        raise ValueError(f"{name} must be a number")
    if minimum is not None and (v < minimum or (exclusive_min and v == minimum)):
        raise ValueError(f"{name} must be {'>' if exclusive_min else '>='} {minimum}")
    if maximum is not None and (v > maximum or (exclusive_max and v == maximum)):
        raise ValueError(f"{name} must be {'<' if exclusive_max else '<='} {maximum}")
    return v


def probability(value, name="alpha"):
    return number(value, name, 0, 1, exclusive_min=True, exclusive_max=True)


def percentage(value, name="percentage"):
    return number(value, name, 0, 100)


def threads(value, max_cores):
    return positive_int(value, "threads", 1, max(1, max_cores))


DESIGN_RE = re.compile(r"^~\s*[A-Za-z][A-Za-z0-9_]*(\s*\+\s*[A-Za-z][A-Za-z0-9_]*)*\s*$")


def design_formula(formula, columns=None):
    """Validate an additive design formula like '~ batch + condition'. Returns variable list."""
    f = (formula or "").strip()
    if not DESIGN_RE.match(f):
        raise ValueError(
            f"invalid design formula {formula!r}. Version 1 supports additive formulas such as "
            "'~ condition' or '~ batch + condition' (variable of interest last).")
    variables = [v.strip() for v in f.lstrip("~").split("+")]
    if len(set(variables)) != len(variables):
        raise ValueError("design formula repeats a variable")
    if columns is not None:
        missing = [v for v in variables if v not in columns]
        if missing:
            raise ValueError(f"design variables not found in sample metadata: {', '.join(missing)}")
    return variables
