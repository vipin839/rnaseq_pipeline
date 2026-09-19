"""GEO series -> samples (GSM) with characteristics -> linked SRA experiments -> runs (via ENA)."""
import re

from . import PipelineError, net

GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def parse_soft(text):
    """Parse GEO SOFT 'brief' text for samples. Returns list of sample dicts."""
    samples, cur = [], None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=", 1)[1].strip(), "characteristics": {}, "srx": []}
            samples.append(cur)
            continue
        if cur is None or " = " not in line:
            continue
        key, val = line.split(" = ", 1)
        key = key.lstrip("!")
        if key == "Sample_title":
            cur["title"] = val
        elif key == "Sample_source_name_ch1":
            cur["source_name"] = val
        elif key == "Sample_organism_ch1":
            cur["organism"] = val
        elif key.startswith("Sample_characteristics_ch1"):
            if ":" in val:
                k, v = val.split(":", 1)
                cur["characteristics"][k.strip()] = v.strip()
            else:
                cur["characteristics"].setdefault("characteristics", val.strip())
        elif key == "Sample_relation" and val.startswith("SRA:"):
            m = re.search(r"(SRX\d+|ERX\d+|DRX\d+)", val)
            if m:
                cur["srx"].append(m.group(1))
        elif key == "Sample_library_strategy":
            cur["library_strategy"] = val
    return samples


def series(gse):
    """Return (title, samples)."""
    head = net.get_text(GEO, {"acc": gse, "targ": "self", "form": "text", "view": "brief"})
    if "^SERIES" not in head:
        raise PipelineError(f"GEO series {gse} not found or not public",
                            remedy="check the accession on https://www.ncbi.nlm.nih.gov/geo/")
    title = next((l.split(" = ", 1)[1] for l in head.splitlines() if l.startswith("!Series_title")), "")
    text = net.get_text(GEO, {"acc": gse, "targ": "gsm", "form": "text", "view": "brief"}, timeout=180)
    samples = parse_soft(text)
    if not samples:
        raise PipelineError(f"no samples found for {gse}")
    return title, samples
