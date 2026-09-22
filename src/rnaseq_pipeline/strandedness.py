"""Library strandedness: RSeQC infer_experiment.py evidence + metadata hints + explicit user confirmation.

Mapping (dUTP/TruSeq Stranded = 'reverse'):
  forward    -> featureCounts -s 1, StringTie --fr   (reads 1++,1--,2+-,2-+ / SE ++,--)
  reverse    -> featureCounts -s 2, StringTie --rf   (reads 1+-,1-+,2++,2-- / SE +-,-+)
  unstranded -> featureCounts -s 0, no StringTie strand flag
"""
import re

from . import runner

FEATURECOUNTS_FLAG = {"unstranded": "0", "forward": "1", "reverse": "2"}
STRINGTIE_FLAG = {"unstranded": None, "forward": "--fr", "reverse": "--rf"}
HINT_WORDS = re.compile(r"(stranded|dUTP|directional|TruSeq Stranded|ScriptSeq|NEBNext Ultra II Directional)", re.I)


def run_infer(bam, bed12, n_reads, log_file, sid):
    res = runner.run(["infer_experiment.py", "-r", bed12, "-i", bam, "-s", str(n_reads)],
                     stage="strandedness", sample=sid, log_file=log_file, capture=True,
                     description=f"{sid}: RSeQC infer_experiment")
    return parse_infer(res.stdout)


def parse_infer(text):
    out = {"undetermined": None, "forward": None, "reverse": None, "paired": "PairEnd" in (text or "")}
    m = re.search(r"failed to determine:\s*([\d.]+)", text or "")
    if m:
        out["undetermined"] = float(m.group(1))
    for line in (text or "").splitlines():
        m = re.search(r'explained by "([^"]+)":\s*([\d.]+)', line)
        if not m:
            continue
        key, val = m.group(1), float(m.group(2))
        if key in ("1++,1--,2+-,2-+", "++,--"):
            out["forward"] = val
        elif key in ("1+-,1-+,2++,2--", "+-,-+"):
            out["reverse"] = val
    return out


def call(ev, stranded_min, unstranded_max_diff, max_undetermined):
    """Return (call or None, confidence, explanation)."""
    f, r, u = ev.get("forward"), ev.get("reverse"), ev.get("undetermined")
    if f is None or r is None:
        return None, "none", "could not parse infer_experiment output"
    if u is not None and u > max_undetermined:
        return None, "low", f"{u:.0%} of reads could not be assigned to a gene strand"
    denom = f + r
    if denom <= 0:
        return None, "none", "no informative reads"
    ff, rr = f / denom, r / denom
    if ff >= stranded_min:
        return "forward", "high" if ff >= 0.9 else "medium", f"{ff:.1%} of informative reads forward-stranded"
    if rr >= stranded_min:
        return "reverse", "high" if rr >= 0.9 else "medium", f"{rr:.1%} of informative reads reverse-stranded"
    if abs(ff - rr) <= unstranded_max_diff:
        return "unstranded", "high" if abs(ff - rr) <= 0.05 else "medium", \
            f"forward {ff:.1%} vs reverse {rr:.1%} (balanced)"
    return None, "low", f"ambiguous: forward {ff:.1%} vs reverse {rr:.1%}"


def consensus(calls):
    """calls: {sample: (call, confidence, why)}. Returns (call or None, confidence, note)."""
    vals = {c[0] for c in calls.values()}
    if None in vals or len(vals) != 1:
        return None, "low", "samples disagree or some could not be determined"
    conf = "high" if all(c[1] == "high" for c in calls.values()) else "medium"
    return vals.pop(), conf, "all samples agree"


def metadata_hints(project):
    hints = []
    for sid, rec in project.samples.items():
        for k, v in (rec.get("metadata") or {}).items():
            if isinstance(v, str) and HINT_WORDS.search(v):
                hints.append(f"{sid}: {k} = {v[:100]}")
    return hints[:10]
