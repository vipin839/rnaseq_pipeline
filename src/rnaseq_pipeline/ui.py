"""Terminal output and prompts. All user interaction goes through here."""
import logging
import shutil
import sys

from . import UserAbort
from .secrets import scrub

log = logging.getLogger("rnaseq")

_COLORS = {
    "INFO": "\033[36m", "OK": "\033[32m", "WARNING": "\033[33m", "ERROR": "\033[31m",
    "FAILED": "\033[31;1m", "RUNNING": "\033[34m", "SKIPPED": "\033[90m", "ACTION": "\033[35m",
    "DRY-RUN": "\033[95m",
}
_RESET = "\033[0m"
_LEVELS = {"WARNING": logging.WARNING, "ERROR": logging.ERROR, "FAILED": logging.ERROR}

VERBOSE = False


def _use_color():
    return sys.stdout.isatty()


def status(tag, message):
    """Print '[TAG] message' and log it."""
    tag = tag.upper()
    label = f"[{tag}]"
    if _use_color() and tag in _COLORS:
        label = f"{_COLORS[tag]}{label}{_RESET}"
    message = scrub(str(message))
    print(f"{label} {message}", flush=True)
    log.log(_LEVELS.get(tag, logging.INFO), "[%s] %s", tag, message)


def info(m): status("INFO", m)
def ok(m): status("OK", m)
def warn(m): status("WARNING", m)
def error(m): status("ERROR", m)
def failed(m): status("FAILED", m)
def running(m): status("RUNNING", m)
def skipped(m): status("SKIPPED", m)
def action(m): status("ACTION", m)


def debug(message):
    log.debug(message)
    if VERBOSE:
        print(f"[DEBUG] {message}", flush=True)


def width():
    return min(shutil.get_terminal_size((80, 20)).columns, 100)


def rule(char="-"):
    print(char * width())


def header(title, subtitle=None):
    print()
    rule("=")
    print(title.center(width()).rstrip())
    if subtitle:
        print(subtitle.center(width()).rstrip())
    rule("=")


def section(title):
    print()
    print(title)
    rule("-")


def kv(pairs, indent=0):
    if not pairs:
        return
    w = max(len(str(k)) for k, _ in pairs)
    for k, v in pairs:
        print(" " * indent + f"{str(k) + ':':<{w + 1}} {v}")


def table(rows, headers, max_col=40):
    if not rows:
        print("  (none)")
        return
    cells = [[str(c) if c is not None else "" for c in r] for r in rows]
    cells = [[c if len(c) <= max_col else c[: max_col - 3] + "..." for c in r] for r in cells]
    widths = [max(len(h), *(len(r[i]) for r in cells)) for i, h in enumerate(headers)]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    print("  ".join("-" * w for w in widths))
    for r in cells:
        print(fmt.format(*r))


def _input(prompt):
    try:
        return input(prompt)
    except EOFError:
        print()
        raise UserAbort("input closed") from None


def ask(prompt, default=None, validator=None, allow_empty=False):
    """Ask for free text. validator(value) returns the cleaned value or raises ValueError."""
    suffix = f" [{default}]" if default not in (None, "") else ""
    while True:
        raw = _input(f"{prompt}{suffix}: ").strip()
        if not raw and default is not None:
            raw = str(default)
        if not raw and not allow_empty:
            print("  A value is required.")
            continue
        if validator:
            try:
                return validator(raw)
            except ValueError as e:
                print(f"  Invalid: {e}")
                continue
        return raw


def ask_yes_no(prompt, default=False):
    hint = "[Y/n]" if default else "[y/N]"
    while True:
        raw = _input(f"{prompt} {hint}: ").strip().lower()
        if not raw:
            return default
        if raw in ("y", "yes"):
            return True
        if raw in ("n", "no"):
            return False
        print("  Please answer y or n.")


def choose(title, options, prompt="Select", default=None):
    """Numbered menu. Returns the 0-based index of the chosen option."""
    if title:
        section(title)
    for i, opt in enumerate(options, 1):
        print(f"  {i}. {opt}")
    rule("-")
    while True:
        d = f" [{default + 1}]" if default is not None else ""
        raw = _input(f"{prompt}{d}: ").strip()
        if not raw and default is not None:
            return default
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return int(raw) - 1
        print(f"  Enter a number between 1 and {len(options)}.")


def choose_many(prompt, items):
    """Select a subset: 'all', or comma/range list like '1,3-5'. Returns list of indices."""
    for i, it in enumerate(items, 1):
        print(f"  {i}. {it}")
    while True:
        raw = _input(f"{prompt} (e.g. all, 1,3-5): ").strip().lower()
        if raw in ("", "all", "a"):
            return list(range(len(items)))
        try:
            idx = set()
            for part in raw.split(","):
                part = part.strip()
                if "-" in part:
                    a, b = part.split("-", 1)
                    idx.update(range(int(a) - 1, int(b)))
                else:
                    idx.add(int(part) - 1)
            if not idx or min(idx) < 0 or max(idx) >= len(items):
                raise ValueError
            return sorted(idx)
        except ValueError:
            print("  Invalid selection.")


def explain(what, why, consequence, options=None):
    """Context shown before a scientific decision: what is decided, why it matters, what follows."""
    print()
    kv([("Decision", what), ("Why it matters", why)] + ([("Options", options)] if options else [])
       + [("Consequence", consequence)])
    print()


def pause():
    _input("Press Enter to continue...")


def explain_failure(err):
    """Print a structured failure block for a PipelineError."""
    print()
    rule("!")
    print("STATUS: FAILED")
    pairs = []
    if getattr(err, "stage", None):
        pairs.append(("Stage", err.stage))
    if getattr(err, "sample", None):
        pairs.append(("Sample", err.sample))
    pairs.append(("Problem", scrub(str(err))))
    if getattr(err, "cause", None):
        pairs.append(("Likely cause", scrub(str(err.cause))))
    if getattr(err, "remedy", None):
        pairs.append(("Recommended action", scrub(str(err.remedy))))
    kv(pairs)
    rule("!")
    log.error("FAILED stage=%s sample=%s problem=%s cause=%s remedy=%s",
              getattr(err, "stage", None), getattr(err, "sample", None), err,
              getattr(err, "cause", None), getattr(err, "remedy", None))
