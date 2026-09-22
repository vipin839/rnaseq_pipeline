"""`rnaseq-pipeline --check`: a health check that actually exercises the installation.

Every check yields PASS / WARNING / FAIL with a remediation. Functional checks run each scientific tool on a
tiny synthetic job (seconds), so a tool that exists but is broken (missing library, wrong Java, corrupt
install) is reported as FAIL rather than passing because the file exists.
"""
import os
import random
import shutil
import subprocess
import tempfile
from pathlib import Path

from . import __version__, dependency_manager, runner, system_check, ui
from . import config as C

PASS, WARN, FAIL = "PASS", "WARNING", "FAIL"


class Report:
    def __init__(self):
        self.rows = []

    def add(self, area, item, status, detail="", fix=""):
        self.rows.append((area, item, status, detail, fix))

    @property
    def overall(self):
        states = {r[2] for r in self.rows}
        return FAIL if FAIL in states else (WARN if WARN in states else PASS)


def _run(cmd, cwd, timeout=120):
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True, timeout=timeout,
                       env=runner.child_env())
    return r.returncode, (r.stdout + r.stderr)[-400:]


def _tiny_dataset(d):
    """2 kb genome, one 2-exon gene, 200 read pairs drawn from it (reverse-stranded)."""
    rng = random.Random(1)
    genome = "".join(rng.choice("ACGT") for _ in range(2000))
    (d / "g.fa").write_text(">c1\n" + "\n".join(genome[i:i + 60] for i in range(0, 2000, 60)) + "\n")
    (d / "a.gtf").write_text(
        'c1\tt\texon\t201\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'c1\tt\texon\t1001\t1600\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    tx = genome[200:700] + genome[1000:1600]
    comp = str.maketrans("ACGT", "TGCA")
    with open(d / "r1.fq", "w") as f1, open(d / "r2.fq", "w") as f2:
        for i in range(200):
            s = rng.randint(0, len(tx) - 250)
            frag = tx[s:s + 250]
            f1.write(f"@r{i}/1\n{frag[-75:].translate(comp)[::-1]}\n+\n{'I' * 75}\n")
            f2.write(f"@r{i}/2\n{frag[:75]}\n+\n{'I' * 75}\n")


FUNCTIONAL = [
    ("hisat2-build", lambda d: ["hisat2-build", "-q", d / "g.fa", d / "idx"]),
    ("hisat2", lambda d: ["hisat2", "-x", d / "idx", "-1", d / "r1.fq", "-2", d / "r2.fq", "-S", d / "a.sam",
                          "--no-unal"]),
    ("samtools", lambda d: ["samtools", "sort", "-o", d / "a.bam", d / "a.sam"]),
    ("samtools index", lambda d: ["samtools", "index", d / "a.bam"]),
    ("featureCounts", lambda d: ["featureCounts", "-p", "--countReadPairs", "-s", "2", "-a", d / "a.gtf",
                                 "-o", d / "fc.txt", d / "a.bam"]),
    ("stringtie", lambda d: ["stringtie", d / "a.bam", "-G", d / "a.gtf", "-e", "-o", d / "st.gtf"]),
    ("fastp", lambda d: ["fastp", "-i", d / "r1.fq", "-I", d / "r2.fq", "-o", d / "t1.fq", "-O", d / "t2.fq",
                         "-j", d / "fp.json", "-h", d / "fp.html"]),
    ("fastqc", lambda d: ["fastqc", "-q", "-o", d, d / "r1.fq"]),
]


def check_functional(rep, fixes=None):
    """fixes: {executable: command} replacing builds that cannot run on this CPU."""
    fixes = fixes or {}
    d = Path(tempfile.mkdtemp(prefix="rnaseq_doctor_"))
    try:
        _tiny_dataset(d)
        broken = set()
        for name, cmd in FUNCTIONAL:
            tool = cmd(d)[0]
            if runner.which(tool) is None:
                continue  # already reported as missing
            if name in ("hisat2", "samtools", "samtools index", "featureCounts", "stringtie") and broken:
                rep.add("Functional test", name, WARN, "skipped (an earlier step failed)")
                continue
            try:
                code, out = _run(cmd(d), d)
            except subprocess.TimeoutExpired:
                code, out = -1, "timed out"
            if code == 0:
                rep.add("Functional test", name, PASS, "ran a tiny job successfully")
            else:
                broken.add(name)
                meaning = runner.describe_exit(code)
                if "SIGILL" in meaning:
                    fix = fixes.get(tool) or ("install a build of this program for older CPUs (main menu 4); "
                                              "reinstalling the same build will not help")
                    rep.add("Functional test", name, FAIL, f"exit {code}: {meaning}", fix)
                    continue
                last = out.strip().splitlines()[-1][:90] if out.strip() else ""
                rep.add("Functional test", name, FAIL, f"exit {code}: {meaning or last}",
                        "the tool is installed but does not work; reinstall the tools environment (main menu 4)")
        fc = d / "fc.txt.summary"
        if fc.exists():
            assigned = next((int(l.split()[1]) for l in fc.read_text().splitlines() if l.startswith("Assigned")), 0)
            rep.add("Functional test", "end-to-end counts", PASS if assigned >= 150 else FAIL,
                    f"{assigned}/200 read pairs assigned to the test gene",
                    "" if assigned >= 150 else "alignment or counting gives wrong results; reinstall tools")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def check_r(rep, rscript):
    if not rscript:
        rep.add("R", "Rscript", FAIL, "not found", "install the R environment (main menu 4)")
        return
    code, out = _run([rscript, "-e", "suppressPackageStartupMessages({library(DESeq2); library(ggplot2); "
                                      "library(pheatmap); library(jsonlite)}); "
                                      "dds <- DESeqDataSetFromMatrix(matrix(c(10L,20L,30L,40L,5L,6L,7L,8L), 2, "
                                      "dimnames=list(c('g1','g2'), paste0('s',1:4))), "
                                      "data.frame(c=factor(c('a','a','b','b')), row.names=paste0('s',1:4)), ~c); "
                                      "cat('R-OK', ncol(dds))"], Path.cwd(), timeout=300)
    if code == 0 and "R-OK 4" in out:
        rep.add("R", "DESeq2 + plotting packages", PASS, "loaded and built a test dataset")
    else:
        rep.add("R", "DESeq2 + plotting packages", FAIL, out.strip().splitlines()[-1][:100] if out.strip() else "",
                "install/repair the R environment (main menu 4)")


def check_network(rep):
    """Reachable = the server answered at all (any HTTP status below 500), not only with 200."""
    import urllib.error
    import urllib.request
    for name, url in (("ENA", "https://www.ebi.ac.uk/ena/portal/api/"),
                      ("NCBI", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?retmode=json")):
        try:
            with urllib.request.urlopen(url, timeout=15) as r:
                r.read(256)
            rep.add("Network", name, PASS, "reachable")
        except urllib.error.HTTPError as e:
            ok = e.code < 500
            rep.add("Network", name, PASS if ok else WARN, f"answered HTTP {e.code}",
                    "" if ok else "the service reported a server error; try again later")
        except (OSError, ValueError) as e:  # DNS failure, refused connection, timeout
            rep.add("Network", name, WARN, f"not reachable ({type(e).__name__})",
                    "needed only for public data and reference downloads; check the internet connection")


def run(cfg, envs, projects_dir, config_sources, quick=False):
    rep = Report()
    info = system_check.collect(projects_dir, r_bin=envs.rscript(), conda_bin=envs.conda_bin())
    crit, warns = system_check.assess(info, cfg=cfg)
    rep.add("System", "OS", PASS if info["os"] == "Linux" else FAIL, info["distribution"],
            "" if info["os"] == "Linux" else "Linux is required")
    rep.add("System", "Python", PASS if info["python_ok"] else FAIL, info["python"], "Python 3.10+ is required")
    rep.add("System", "CPU", PASS if info["cpu_cores"] >= 2 else WARN, f"{info['cpu_cores']} cores")
    lacking = system_check.missing_x86_64_v3()
    rep.add("System", "CPU instruction set", PASS,
            f"older than x86-64-v3 (no {', '.join(lacking)}); tools are checked for builds that run here" if lacking
            else "x86-64-v3 or newer (AVX2/BMI2)")
    rep.add("System", "RAM", PASS if info["ram_total_gb"] >= 8 else WARN, f"{info['ram_total_gb']} GB",
            "8 GB+ recommended; mammalian genomes need ~8 GB for alignment")
    disk_state = FAIL if any("insufficient disk" in c for c in crit) else (
        WARN if info["disk_free_gb"] < cfg.get("storage", {}).get("min_free_gb_warn", 50) else PASS)
    rep.add("System", "Disk (projects)", disk_state, f"{info['disk_free_gb']} GB free at {info['project_path']}")
    fs_bad = info["filesystem"] in ("vfat", "msdos", "fat")
    rep.add("System", "Filesystem", FAIL if fs_bad else (WARN if info["fs_note"] else PASS),
            f"{info['filesystem']}" + (f" — {info['fs_note']}" if info["fs_note"] else ""),
            "use an ext4/xfs location for projects" if fs_bad or info["fs_note"] else "")
    rep.add("System", "Write permission", PASS if info["project_writable"] else FAIL, info["project_path"],
            "choose a writable projects directory (--projects-dir or config projects_dir)")
    rep.add("System", "Temporary space", PASS if info["temp_writable"] else FAIL,
            f"{info['temp_dir']} ({info['temp_free_gb']} GB)")
    errs = C.validate(cfg, info["cpu_cores"])
    rep.add("Configuration", "values", FAIL if errs else PASS, "; ".join(errs)[:120] if errs else
            f"valid ({len(config_sources)} file(s): {', '.join(Path(s).name for s in config_sources)})",
            "fix the listed values" if errs else "")
    email = (cfg.get("ncbi") or {}).get("email")
    rep.add("Configuration", "NCBI email", PASS if email else WARN, "set" if email else "not set",
            "" if email else f"add ncbi: email: to {C.USER_CONFIG} (NCBI asks tools to identify themselves)")
    store = Path(os.path.expanduser(cfg.get("reference_store", "~/rnaseq_references")))
    try:
        store.mkdir(parents=True, exist_ok=True)
        ok = os.access(store, os.W_OK)
    except OSError:
        ok = False
    rep.add("Configuration", "reference store", PASS if ok else FAIL, str(store),
            "" if ok else "set reference_store to a writable directory")
    conda = envs.conda_bin()
    rep.add("Environment", "conda/mamba", PASS if conda else WARN, str(conda or "not found"),
            "" if conda else "only needed to install/repair tools from the menu")
    tools = dependency_manager.detect_tools()
    fixes = {}
    for t in tools:
        state = {"AVAILABLE": PASS, "OPTIONAL": WARN}.get(t["status"], FAIL)
        detail = f"{t['version'] or '-'} ({t['status'].lower()})"
        fix = "" if state == PASS else ("optional: " + t["purpose"] if state == WARN else "install/repair via main menu 4")
        if t.get("cpu_incompatible"):
            detail = t["note"]
            cmd = dependency_manager.cpu_fix_command(envs, t)
            fix = (f"install a build that runs on this CPU: {cmd}  (or main menu 4)" if cmd else
                   "no build for this CPU is known; reinstalling the same build will not help")
            if cmd:
                fixes[t["exe"]] = fix
        rep.add("Tools", t["key"], state, detail, fix)
    r_info, r_pkgs = dependency_manager.detect_r(envs.rscript())
    rep.add("R", "R version", PASS if r_info["status"] == "AVAILABLE" else FAIL, r_info.get("version") or "missing",
            "" if r_info["status"] == "AVAILABLE" else "R 4.x is required for DESeq2")
    for p in r_pkgs:
        if p["status"] != "AVAILABLE":
            rep.add("R", p["package"], FAIL if p["required"] else WARN, "missing", "install via main menu 4")
    if not quick:
        check_functional(rep, fixes)
        check_r(rep, envs.rscript())
        check_network(rep)
    return rep


def show(rep):
    ui.header(f"RNA-seq PIPELINE {__version__} — HEALTH CHECK")
    area = None
    for a, item, status, detail, fix in rep.rows:
        if a != area:
            ui.section(a.upper())
            area = a
        tag = {"PASS": "OK", "WARNING": "WARNING", "FAIL": "FAILED"}[status]
        ui.status(tag, f"{item:<26} {detail}")
        if fix and status != PASS:
            print(f"{'':12}-> {fix}")
    n = {s: sum(1 for r in rep.rows if r[2] == s) for s in (PASS, WARN, FAIL)}
    print()
    ui.rule("=")
    print(f"HEALTH: {rep.overall}   ({n[PASS]} passed, {n[WARN]} warnings, {n[FAIL]} failed)")
    ui.rule("=")
