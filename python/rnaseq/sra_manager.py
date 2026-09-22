"""NCBI SRA route: prefetch -> vdb-validate -> fasterq-dump -> compress."""
import os
import re
import shutil
from pathlib import Path

from . import PipelineError, net, runner, ui
from .ena_manager import parse_tsv

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def runinfo(accession):
    """NCBI runinfo (fallback metadata source). Returns list of dicts."""
    s = net.get_text(f"{EUTILS}/esearch.fcgi", {"db": "sra", "term": accession, "usehistory": "y",
                                                  "retmode": "json"})
    import json
    j = json.loads(s)["esearchresult"]
    if int(j.get("count", 0)) == 0:
        return []
    text = net.get_text(f"{EUTILS}/efetch.fcgi", {"db": "sra", "query_key": j["querykey"],
                                                  "WebEnv": j["webenv"], "rettype": "runinfo",
                                                  "retmode": "text"})
    rows = parse_tsv(text.replace(",", "\t"))
    return [r for r in rows if r.get("Run")]


def download_run(acc, raw_dir, fastq_dir, temp_dir, threads, paired, *, max_reads=0, log_file=None):
    raw_dir, fastq_dir, temp_dir = Path(raw_dir), Path(fastq_dir), Path(temp_dir)
    finals = ([fastq_dir / f"{acc}_1.fastq.gz", fastq_dir / f"{acc}_2.fastq.gz"] if paired
              else [fastq_dir / f"{acc}.fastq.gz"])
    if all(f.exists() for f in finals):
        ui.skipped(f"{acc}: FASTQ already present")
        return finals[0], (finals[1] if paired else None)

    work = temp_dir / f"sra_{acc}"
    if work.exists():
        shutil.rmtree(work)  # temp workspace from an interrupted run; never contains the only copy of data
    work.mkdir(parents=True)

    if max_reads:
        runner.run(["fastq-dump", "-X", str(max_reads), "--split-3", "--gzip", "-O", work, acc],
                   stage="data", sample=acc, log_file=log_file,
                   description=f"{acc}: fastq-dump first {max_reads:,} spots (PILOT subset)")
        produced = sorted(work.glob(f"{acc}*.fastq.gz"))
    else:
        sra_dir = raw_dir / acc
        runner.run(["prefetch", "--max-size", "u", "-O", raw_dir, acc], stage="data", sample=acc,
                   log_file=log_file, description=f"{acc}: prefetch")
        sra_file = sra_dir / f"{acc}.sra"
        if not sra_file.exists():
            cand = list(sra_dir.glob("*.sra*")) if sra_dir.exists() else []
            if not cand and not runner.DRY_RUN:
                raise PipelineError(f"prefetch produced no .sra file for {acc}", sample=acc, stage="data")
            sra_file = cand[0] if cand else sra_file
        runner.run(["vdb-validate", sra_dir], stage="data", sample=acc, log_file=log_file,
                   description=f"{acc}: validating .sra (vdb-validate)")
        runner.run(["fasterq-dump", "--split-3", "--skip-technical", "-e", str(threads), "-t", work,
                    "-O", work, sra_file], stage="data", sample=acc, log_file=log_file,
                   description=f"{acc}: fasterq-dump")
        for fq in sorted(work.glob("*.fastq")):
            comp = ["pigz", "-p", str(threads)] if runner.which("pigz") else ["gzip"]
            runner.run([*comp, fq], stage="data", sample=acc, log_file=log_file,
                       description=f"{acc}: compressing {fq.name}")
        produced = sorted(work.glob("*.fastq.gz"))
    if runner.DRY_RUN:
        return finals[0], (finals[1] if paired else None)
    names = {p.name: p for p in produced}
    extra_mates = sorted(n for n in names if re.match(rf"{acc}_[3-9]\.fastq\.gz$", n))
    if extra_mates:
        shutil.rmtree(work, ignore_errors=True)
        raise PipelineError(f"{acc} has more than two reads per spot ({len(extra_mates) + 2} or more)",
                            sample=acc, stage="data",
                            cause="index/cell-barcode/UMI reads are stored with the cDNA read (typical of 10x "
                                  "Genomics single-cell libraries), so this is not a standard bulk library",
                            remedy="analyse single-cell data with Cell Ranger / STARsolo / alevin-fry")
    if paired:
        want = [f"{acc}_1.fastq.gz", f"{acc}_2.fastq.gz"]
        if not all(w in names for w in want):
            raise PipelineError(f"{acc}: expected paired files {want}, got {sorted(names)}", sample=acc,
                                stage="data", cause="the run may actually be single-end",
                                remedy="check library layout in the metadata")
    else:
        want = [f"{acc}.fastq.gz"] if f"{acc}.fastq.gz" in names else [f"{acc}_1.fastq.gz"]
        if want[0] not in names:
            raise PipelineError(f"{acc}: no FASTQ produced", sample=acc, stage="data")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    for w, final in zip(want, finals):
        os.replace(names[w], final)
    leftovers = [n for n in names if n not in want]
    if leftovers:
        ui.info(f"{acc}: ignoring unpaired/technical leftovers {leftovers}")
    shutil.rmtree(work, ignore_errors=True)
    return finals[0], (finals[1] if paired else None)
