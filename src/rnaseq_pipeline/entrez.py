"""NCBI Entrez (E-utilities) client: dataset/run/assembly discovery.

Used for SEARCHING NCBI (GEO series, SRA runs, taxonomy). Downloads still go through
ENA (direct FASTQ URLs + MD5) or the SRA toolkit, which are better suited to bulk files.

NCBI asks callers to identify themselves and limits requests to 3/second (10/second with a
free API key). Both are configurable under `ncbi:` in the configuration.
"""
import csv
import io
import json
import threading
import time
import urllib.parse

from . import PipelineError, __version__, net, secrets

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = f"rnaseq-pipeline/{__version__}"

_state = {"email": "", "api_key": "", "min_interval": 1 / 3, "last": 0.0}
_lock = threading.Lock()


def configure(cfg):
    n = (cfg or {}).get("ncbi") or {}
    _state["email"] = (n.get("email") or "").strip()
    _state["api_key"] = (n.get("api_key") or "").strip()
    secrets.register(_state["api_key"])
    _state["min_interval"] = 1 / 10 if _state["api_key"] else 1 / 3


def _throttle():
    with _lock:
        wait = _state["min_interval"] - (time.time() - _state["last"])
        if wait > 0:
            time.sleep(wait)
        _state["last"] = time.time()


def _params(extra):
    p = {"tool": TOOL, **extra}
    if _state["email"]:
        p["email"] = _state["email"]
    if _state["api_key"]:
        p["api_key"] = _state["api_key"]
    return p


def _get(endpoint, params, retries=3):
    _throttle()
    return net.get_text(f"{BASE}/{endpoint}", _params(params), retries=retries)


def _json(endpoint, params):
    text = _get(endpoint, {**params, "retmode": "json"})
    try:
        return json.loads(text)
    except ValueError:
        raise PipelineError(f"unexpected (non-JSON) response from NCBI {endpoint}",
                            cause=text[:200], remedy="try again in a moment; NCBI may be rate-limiting") from None


# ------------------------------------------------------------------ einfo (searchable indexes)

def databases():
    return _json("einfo.fcgi", {})["einforesult"]["dblist"]


def fields(db):
    """Searchable indexes of a database: [{name, fullname, description}] (e.g. ORGN, ACCN, FILT)."""
    info = _json("einfo.fcgi", {"db": db})["einforesult"]["dbinfo"][0]
    return [{"name": f["name"], "fullname": f["fullname"], "description": f.get("description", "")}
            for f in info.get("fieldlist", [])], int(info.get("count", 0))


# ------------------------------------------------------------------ search / summary

def search(db, term, retmax=20, retstart=0):
    r = _json("esearch.fcgi", {"db": db, "term": term, "retmax": retmax, "retstart": retstart})["esearchresult"]
    if "ERROR" in r:
        raise PipelineError(f"NCBI rejected the query: {r['ERROR']}",
                            remedy="check the field names and quoting, e.g. \"Homo sapiens\"[Organism]")
    return {"count": int(r.get("count", 0)), "ids": r.get("idlist", []),
            "translation": r.get("querytranslation", "")}


def summaries(db, ids):
    if not ids:
        return []
    d = _json("esummary.fcgi", {"db": db, "id": ",".join(ids)}).get("result", {})
    return [d[u] for u in d.get("uids", []) if u in d]


def geo_series(term, retmax=20, retstart=0):
    """Search GEO DataSets, returning series (GSE) records only."""
    full = f"({term}) AND gse[Entry Type]" if "entry type" not in term.lower() else term
    res = search("gds", full, retmax, retstart)
    out = []
    for s in summaries("gds", res["ids"]):
        out.append({"accession": s.get("accession", ""), "title": s.get("title", ""),
                    "organism": s.get("taxon", ""), "samples": int(s.get("n_samples") or 0),
                    "type": s.get("gdstype", ""), "date": s.get("pdat", ""),
                    "summary": (s.get("summary") or "")[:400],
                    "bioproject": s.get("bioproject", "")})
    return res["count"], out, res["translation"]


RUNINFO_FIELDS = ["Run", "spots", "bases", "size_MB", "LibraryLayout", "LibraryStrategy", "LibrarySource",
                  "Platform", "Model", "ScientificName", "SampleName", "Experiment", "BioProject", "BioSample"]


def runinfo(term, retmax=500):
    """Run-level metadata for an SRA query or accession (CSV 'runinfo' report)."""
    res = _json("esearch.fcgi", {"db": "sra", "term": term, "retmax": 0, "usehistory": "y"})["esearchresult"]
    if int(res.get("count", 0)) == 0:
        return 0, []
    text = _get("efetch.fcgi", {"db": "sra", "query_key": res["querykey"], "WebEnv": res["webenv"],
                                "rettype": "runinfo", "retmode": "text", "retmax": retmax})
    rows = [r for r in csv.DictReader(io.StringIO(text)) if r.get("Run")]
    return int(res["count"]), rows


def taxonomy_suggest(name, retmax=8):
    """Resolve an organism name to NCBI taxonomy records (helps pick the right strain/species)."""
    res = search("taxonomy", f"{name}[All Names]", retmax=retmax)
    out = []
    for s in summaries("taxonomy", res["ids"]):
        out.append({"taxid": s.get("taxid") or s.get("uid"), "name": s.get("scientificname", ""),
                    "rank": s.get("rank", ""), "division": s.get("division", "")})
    return out


def build_term(parts):
    """Join ('value', 'FIELD') pairs into an Entrez query: value[FIELD] AND ..."""
    terms = []
    for value, field in parts:
        value = value.strip()
        if not value:
            continue
        if " " in value and not value.startswith('"'):
            value = f'"{value}"'
        terms.append(f"{value}[{field}]" if field else value)
    return " AND ".join(terms)


def link(dbfrom, db, ids):
    """elink: e.g. GEO series UIDs -> SRA UIDs."""
    text = _get("elink.fcgi", {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "json"})
    try:
        sets = json.loads(text).get("linksets", [])
    except ValueError:
        return []
    out = []
    for s in sets:
        for ldb in s.get("linksetdbs", []):
            out += ldb.get("links", [])
    return out


def quote(s):
    return urllib.parse.quote(str(s))
