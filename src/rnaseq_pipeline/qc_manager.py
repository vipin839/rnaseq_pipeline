"""FastQC and MultiQC execution with post-run validation."""
import json
import re
import zipfile
from pathlib import Path

from . import PipelineError, runner, ui
from . import config as C

FQ_EXT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)


def fastqc_basename(fastq):
    name = Path(fastq).name
    return FQ_EXT.sub("", name) + "_fastqc"


def run_fastqc(files, out_dir, temp_dir, threads, log_file, stage):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    todo = [f for f in files if not fastqc_ok(f, out_dir)]
    if len(todo) < len(files):
        ui.skipped(f"FastQC: {len(files) - len(todo)} file(s) already have valid reports")
    if todo:
        runner.run(["fastqc", "--threads", str(min(threads, len(todo))), "--outdir", out_dir,
                    "--dir", temp_dir, "--noextract", *todo],
                   stage=stage, log_file=log_file, description=f"FastQC on {len(todo)} file(s)")
    return [out_dir / (fastqc_basename(f) + ".zip") for f in files]


def read_fastqc_zip(zip_path):
    """Return (summary_text, data_text) or raise."""
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        data = next((n for n in names if n.endswith("/fastqc_data.txt")), None)
        summ = next((n for n in names if n.endswith("/summary.txt")), None)
        if not data or not summ:
            raise ValueError("fastqc_data.txt or summary.txt missing in zip")
        return z.read(summ).decode(), z.read(data).decode()


def fastqc_total_sequences(data_text):
    m = re.search(r"^Total Sequences\t(\d+)", data_text, re.M)
    return int(m.group(1)) if m else None


def fastqc_problem(fastq, out_dir, expected_reads=None):
    """Why the FastQC report of `fastq` is not acceptable, or None. The report must be complete, belong to THIS
    file (FastQC records the file name it analysed) and count the reads the FASTQ validation counted."""
    base = fastqc_basename(fastq)
    html, z = Path(out_dir) / f"{base}.html", Path(out_dir) / f"{base}.zip"
    if not html.exists() or html.stat().st_size == 0 or not z.exists():
        return "report missing"
    try:
        _, data = read_fastqc_zip(z)
    except (zipfile.BadZipFile, ValueError, OSError):
        return "report unreadable"
    if ">>END_MODULE" not in data:
        return "report incomplete"
    m = re.search(r"^Filename\t(.+)$", data, re.M)
    if not m or m.group(1).strip() != Path(fastq).name:
        return f"report belongs to {m.group(1).strip() if m else 'an unknown file'}, not {Path(fastq).name}"
    n = fastqc_total_sequences(data)
    if n is None or n == 0:
        return "report counts no reads"
    if expected_reads is not None and n != expected_reads:
        return f"read count {n} != {expected_reads} from FASTQ validation"
    return None


def fastqc_ok(fastq, out_dir, expected_reads=None):
    return fastqc_problem(fastq, out_dir, expected_reads) is None


def validate_fastqc(files, out_dir, expected_reads, stage):
    """expected_reads: {fastq_path_str: reads}. Raises PipelineError on any invalid report."""
    problems = []
    for f in files:
        why = fastqc_problem(f, out_dir, expected_reads.get(str(f)))
        if why:
            problems.append(f"{Path(f).name}: FastQC {why}")
    if problems:
        raise PipelineError("FastQC output validation failed:\n  " + "\n  ".join(problems), stage=stage,
                            cause="FastQC crashed, ran out of memory, or input changed",
                            remedy="see the FastQC log; re-run the stage")
    ui.ok(f"FastQC output validated ({len(files)} reports; read counts match FASTQ validation)")


MULTIQC_CONFIG = {
    "title": None,
    "show_analysis_paths": False,
    "fn_clean_exts": ["_fastqc", ".fastq.gz", ".fq.gz", ".fastq", ".fq", ".hisat2.summary", ".bam"],
}


def run_multiqc(in_dirs, out_dir, title, log_file, stage, expected_sources=()):
    """expected_sources: input files that MUST appear in the report. MultiQC skips files it cannot parse (and
    overwrites samples with duplicate names) and still exits 0, so a sample could silently disappear."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = dict(MULTIQC_CONFIG, title=title)
    cfg_path = out_dir / "multiqc_config.yaml"
    C.save_yaml(cfg, cfg_path)
    runner.run(["multiqc", "--force", "--no-ansi", "--outdir", out_dir, "--filename", "multiqc_report.html",
                "--config", cfg_path, *in_dirs],
               stage=stage, log_file=log_file, description=f"MultiQC: {title}")
    if runner.DRY_RUN:
        return out_dir / "multiqc_report.html"
    report = out_dir / "multiqc_report.html"
    data_dir = out_dir / "multiqc_report_data"
    if not data_dir.is_dir():
        data_dir = next((d for d in out_dir.iterdir() if d.is_dir() and d.name.endswith("_data")), None)
    if not report.exists() or report.stat().st_size < 1000 or data_dir is None:
        raise PipelineError("MultiQC did not produce a valid report", stage=stage,
                            remedy=f"see {log_file}")
    if expected_sources:
        src = data_dir / "multiqc_sources.txt"
        seen = set()
        if src.exists():
            for line in src.read_text().splitlines()[1:]:
                cells = line.split("\t")
                if len(cells) >= 4:
                    seen.add(str(Path(cells[3]).resolve()))
        missing = [Path(s) for s in expected_sources if str(Path(s).resolve()) not in seen]
        if missing:
            raise PipelineError(f"MultiQC report is missing {len(missing)} of {len(expected_sources)} input(s): "
                                + ", ".join(m.name for m in missing[:6]), stage=stage,
                                cause="MultiQC skips files it cannot parse, and overwrites samples with the same name",
                                remedy=f"see {log_file}; re-run the step (the listed files may be damaged)")
    gs = next((p for p in data_dir.glob("multiqc_general_stats.txt")), None)
    summary = {"title": title, "inputs": [str(d) for d in in_dirs], "report": str(report),
               "data_dir": str(data_dir), "general_stats_rows": 0}
    if gs:
        summary["general_stats_rows"] = max(0, len(gs.read_text().splitlines()) - 1)
    (out_dir / "multiqc_summary.json").write_text(json.dumps(summary, indent=2))
    ui.ok(f"MultiQC report validated: {report}")
    return report
