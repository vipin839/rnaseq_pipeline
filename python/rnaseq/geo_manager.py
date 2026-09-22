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


def _relations(head_text):
    """Parse !Series_relation lines -> {'reanalysis_gsm': [...], 'subseries': [...]}."""
    out = {"reanalysis_gsm": [], "subseries": []}
    for line in head_text.splitlines():
        if not line.startswith("!Series_relation"):
            continue
        val = line.split(" = ", 1)[1] if " = " in line else ""
        if val.startswith("Reanalysis of:"):
            out["reanalysis_gsm"] += re.findall(r"GSM\d+", val)
        elif val.startswith("SuperSeries of:"):
            out["subseries"] += re.findall(r"GSE\d+", val)
    return out


def _series_head(gse):
    return net.get_text(GEO, {"acc": gse, "targ": "self", "form": "text", "view": "brief"})


def _series_samples(gse):
    return parse_soft(net.get_text(GEO, {"acc": gse, "targ": "gsm", "form": "text", "view": "brief"},
                                   timeout=180))


def _resolve_gsms(gsms, notes):
    """Fetch the given GSM samples efficiently: find each one's parent series and load that series once."""
    wanted, found = set(gsms), {}
    while wanted - set(found):
        probe = sorted(wanted - set(found))[0]
        text = net.get_text(GEO, {"acc": probe, "targ": "self", "form": "text", "view": "brief"})
        parents = re.findall(r"!Sample_series_id = (GSE\d+)", text)
        loaded = False
        for parent in parents:
            batch = {s["gsm"]: s for s in _series_samples(parent) if s["gsm"] in wanted}
            if batch:
                found.update(batch)
                notes.append(f"{len(batch)} sample(s) taken from the original series {parent}")
                loaded = True
                break
        if not loaded:  # sample not listed under any parent series: parse it on its own
            one = parse_soft(text if text.startswith("^SAMPLE") else "^SAMPLE = " + probe + "\n" + text)
            found[probe] = one[0] if one else {"gsm": probe, "characteristics": {}, "srx": []}
    return [found[g] for g in gsms if g in found]


def series(gse, _depth=0):
    """Return (title, samples, notes). Follows SuperSeries and third-party-reanalysis links."""
    head = _series_head(gse)
    if "^SERIES" not in head:
        raise PipelineError(f"GEO series {gse} not found or not public",
                            remedy="check the accession on https://www.ncbi.nlm.nih.gov/geo/")
    title = next((l.split(" = ", 1)[1] for l in head.splitlines() if l.startswith("!Series_title")), "")
    notes = []
    samples = _series_samples(gse)
    if samples:
        return title, samples, notes
    rel = _relations(head)
    if rel["subseries"] and _depth < 2:
        notes.append(f"{gse} is a SuperSeries; collecting samples from its sub-series "
                     f"{', '.join(rel['subseries'])}")
        seen = {}
        for sub in rel["subseries"]:
            _, sub_samples, sub_notes = series(sub, _depth + 1)
            notes += sub_notes
            for smp in sub_samples:
                seen.setdefault(smp["gsm"], smp)
        samples = list(seen.values())
    elif rel["reanalysis_gsm"]:
        notes.append(f"{gse} is a third-party REANALYSIS: it has no samples of its own and re-uses "
                     f"{len(rel['reanalysis_gsm'])} sample(s) from other GEO series. Those original samples "
                     "(and their raw reads) will be used.")
        samples = _resolve_gsms(rel["reanalysis_gsm"], notes)
    if not samples:
        types = [l.split(" = ", 1)[1] for l in head.splitlines() if l.startswith("!Series_type")]
        raise PipelineError(f"GEO series {gse} has no samples with sequencing data",
                            cause=f"series type: {', '.join(types) or 'unknown'}; it may be non-sequencing "
                                  "(e.g. microarray) or its samples are not yet public",
                            remedy="open the series on the GEO website and use the accession of the study "
                                   "that holds the raw data")
    return title, samples, notes
