"""ENA Portal API: run metadata for SRA/ENA accessions and direct FASTQ download with MD5 check."""
import gzip
import subprocess
from pathlib import Path

from . import PipelineError, net, runner, ui

FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ["run_accession", "experiment_accession", "sample_accession", "study_accession",
          "secondary_study_accession", "library_layout", "library_strategy", "library_source",
          "library_selection", "library_name", "instrument_platform", "instrument_model", "read_count",
          "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes", "scientific_name", "tax_id", "study_title",
          "sample_title", "sample_alias", "experiment_title", "run_alias", "library_construction_protocol"]


def parse_tsv(text):
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def runs_for(accession):
    """Return list of run dicts for any run/experiment/sample/study accession."""
    text = net.get_text(FILEREPORT, {"accession": accession, "result": "read_run",
                                     "fields": ",".join(FIELDS), "format": "tsv"})
    runs = parse_tsv(text)
    for r in runs:
        r["_files"] = fastq_files(r)
    return runs


def fastq_files(run):
    """[(url, md5, bytes)] for the mate files; unpaired leftovers (SRRxxx.fastq.gz in a PE run) dropped."""
    urls = [u for u in run.get("fastq_ftp", "").split(";") if u]
    md5s = run.get("fastq_md5", "").split(";")
    sizes = run.get("fastq_bytes", "").split(";")
    files = []
    for i, u in enumerate(urls):
        files.append(("https://" + u if not u.startswith("http") else u,
                      md5s[i] if i < len(md5s) else None,
                      int(sizes[i]) if i < len(sizes) and sizes[i].isdigit() else None))
    if run.get("library_layout") == "PAIRED":
        mates = [f for f in files if f[0].endswith(("_1.fastq.gz", "_2.fastq.gz"))]
        if len(mates) == 2:
            return sorted(mates, key=lambda f: f[0])
    if run.get("library_layout") == "SINGLE" and len(files) > 1:
        return files[:1]
    return files


def download_run(run, out_dir, *, max_reads=0, log_file=None, retries=3):
    """Download a run's FASTQ files into out_dir. Returns (r1, r2|None)."""
    acc = run["run_accession"]
    files = run["_files"]
    paired = run.get("library_layout") == "PAIRED"
    if not files:
        raise PipelineError(f"ENA lists no FASTQ files for {acc}", sample=acc, stage="data",
                            remedy="try the SRA download route (config download.source_preference: sra)")
    if paired and len(files) != 2:
        raise PipelineError(f"{acc} is PAIRED but ENA lists {len(files)} FASTQ files", sample=acc, stage="data",
                            remedy="try the SRA route (fasterq-dump --split-3)")
    out = []
    for i, (url, md5, size) in enumerate(files):
        mate = f"_{i + 1}" if paired else ""
        dest = Path(out_dir) / f"{acc}{mate}.fastq.gz"
        if max_reads:
            out.append(_download_head(url, dest, max_reads, acc, log_file))
        else:
            out.append(net.download(url, dest, stage="data", sample=acc, expected_md5=md5 or None,
                                    expected_bytes=size or None,
                                    log_file=log_file, retries=retries))
    return out[0], (out[1] if paired else None)


def _download_head(url, dest, n_reads, acc, log_file):
    """PILOT MODE: stream only the first n_reads records (no full download)."""
    dest = Path(dest)
    if dest.exists():
        ui.skipped(f"already downloaded (pilot subset): {dest.name}")
        return dest
    ui.running(f"Streaming first {n_reads:,} reads of {dest.name} (PILOT subset)")
    if runner.DRY_RUN:
        ui.status("DRY-RUN", f"curl -sL {url} | head -n {4 * n_reads} > {dest}")
        return dest
    part = dest.with_name(dest.name + ".part")
    p = subprocess.Popen(["curl", "-sfL", "--retry", "3", url], stdout=subprocess.PIPE, env=runner.child_env())
    lines = 0
    try:
        with gzip.GzipFile(fileobj=p.stdout) as src, gzip.open(part, "wb", compresslevel=4) as dst:
            for line in src:
                dst.write(line)
                lines += 1
                if lines >= 4 * n_reads:
                    break
    finally:
        p.kill()
        p.wait()
    if lines % 4 != 0 or lines == 0:
        part.unlink(missing_ok=True)
        raise PipelineError(f"pilot download of {dest.name} incomplete ({lines} lines)", sample=acc, stage="data")
    part.rename(dest)
    from . import logger
    logger.record_command({"stage": "data", "sample": acc, "status": "ok", "exit_codes": [0],
                           "command_str": f"[pilot] first {n_reads} reads of {url} -> {dest}"})
    return dest
