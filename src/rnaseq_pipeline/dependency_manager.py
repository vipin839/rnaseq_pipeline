"""Detection, reporting and (user-approved) installation of required software."""
import re
from pathlib import Path

from . import runner, ui

# key, executable, version command, min version, conda package, env, required, purpose
TOOLS = [
    ("fastqc", "fastqc", ["--version"], "0.11.9", "fastqc", "tools", True,
     "Per-file read quality reports (FastQC)"),
    ("multiqc", "multiqc", ["--version"], "1.14", "multiqc", "tools", True,
     "Aggregates FastQC/fastp/HISAT2/featureCounts reports into one HTML"),
    ("fastp", "fastp", ["--version"], "0.23.0", "fastp", "tools", True,
     "Adapter and quality trimming (only if the quality gate recommends it)"),
    ("hisat2", "hisat2", ["--version"], "2.2.0", "hisat2", "tools", True,
     "Splice-aware alignment of RNA-seq reads to the genome"),
    ("hisat2-build", "hisat2-build", ["--version"], "2.2.0", "hisat2", "tools", True,
     "Builds the HISAT2 genome index"),
    ("samtools", "samtools", ["--version"], "1.15", "samtools", "tools", True,
     "BAM sorting, indexing and integrity/QC statistics"),
    ("stringtie", "stringtie", ["--version"], "2.1.0", "stringtie", "tools", True,
     "Reference-guided transcript-level quantification (StringTie2)"),
    ("featureCounts", "featureCounts", ["-v"], "2.0.0", "subread", "tools", True,
     "Gene-level read counting -> count matrix for DESeq2"),
    ("prefetch", "prefetch", ["--version"], "3.0.0", "sra-tools", "tools", True,
     "Downloads .sra files from NCBI SRA"),
    ("fasterq-dump", "fasterq-dump", ["--version"], "3.0.0", "sra-tools", "tools", True,
     "Converts .sra files to FASTQ"),
    ("vdb-validate", "vdb-validate", ["--version"], "3.0.0", "sra-tools", "tools", True,
     "Verifies integrity of downloaded .sra files"),
    ("infer_experiment.py", "infer_experiment.py", ["--version"], "4.0.0", "rseqc", "tools", False,
     "RSeQC: infers library strandedness from alignments"),
    ("pigz", "pigz", ["--version"], "2.4", "pigz", "tools", False,
     "Parallel gzip (faster compression of FASTQ files)"),
    ("curl", "curl", ["--version"], "7.0", "curl", "tools", True,
     "Downloads references and ENA FASTQ files"),
]

# Tools whose current Bioconda builds are compiled for x86-64-v3 (AVX2/BMI2) and die with SIGILL on older CPUs,
# and the newest version whose build does not use those instructions (checked by disassembling the binaries:
# stringtie 3.0.0-3.0.3 all contain AVX2/BMI2 instructions; 2.2.1 and 2.2.3 contain none).
OLD_CPU_VERSIONS = {"stringtie": "2.2.3"}

R_PACKAGES = [
    ("DESeq2", "bioconductor-deseq2", True, "Differential expression statistics"),
    ("ggplot2", "r-ggplot2", True, "Plots"),
    ("pheatmap", "r-pheatmap", True, "Heatmaps"),
    ("EnhancedVolcano", "bioconductor-enhancedvolcano", False, "Volcano plots (ggplot2 fallback if absent)"),
    ("apeglm", "bioconductor-apeglm", False, "LFC shrinkage for MA plots / ranking"),
    ("data.table", "r-data.table", True, "Fast table I/O"),
    ("dplyr", "r-dplyr", True, "Data manipulation"),
    ("tibble", "r-tibble", True, "Data frames"),
    ("readr", "r-readr", True, "Table I/O"),
    ("RColorBrewer", "r-rcolorbrewer", True, "Colour palettes"),
    ("ggrepel", "r-ggrepel", True, "Non-overlapping plot labels"),
    ("jsonlite", "r-jsonlite", True, "Reads the parameter file passed from Python"),
]

VERSION_RE = re.compile(r"(\d+\.\d+(?:\.\d+)?)")


def vtuple(v):
    return tuple(int(x) for x in v.split(".")) if v else ()


def parse_version(text, key):
    if not text:
        return None
    if key == "hisat2" or key == "hisat2-build":
        m = re.search(r"version\s+(\d+\.\d+\.\d+)", text)
        return m.group(1) if m else None
    for line in text.splitlines():
        m = VERSION_RE.search(line)
        if m:
            return m.group(1)
    return None


def detect_tools():
    rows = []
    for key, exe, vcmd, minv, pkg, env, required, purpose in TOOLS:
        path = runner.which(exe)
        rec = {"key": key, "exe": exe, "path": str(path) if path else None, "version": None,
               "min_version": minv, "package": pkg, "env": env, "required": required, "purpose": purpose}
        if not path:
            rec["status"] = "MISSING" if required else "OPTIONAL"
        else:
            code, out = runner.tool_run([path, *vcmd])
            v = parse_version(out, key)
            rec["version"] = v
            meaning = runner.describe_exit(code)
            if "SIGILL" in (meaning or ""):
                rec["status"] = "INCOMPATIBLE"
                rec["note"] = meaning
                rec["cpu_incompatible"] = True
                if key in OLD_CPU_VERSIONS:
                    rec["install_spec"] = f"{pkg}={OLD_CPU_VERSIONS[key]}"
            elif v is None:
                rec["status"] = "INCOMPATIBLE"
                rec["note"] = ("installed but does not run: " + meaning if meaning else
                               "installed but version could not be determined (tool may be broken)")
            elif vtuple(v) < vtuple(minv):
                rec["status"] = "OUTDATED"
            else:
                rec["status"] = "AVAILABLE"
        rows.append(rec)
    return rows


def detect_r(rscript):
    if not rscript:
        return {"rscript": None, "version": None, "status": "MISSING"}, [
            {"package": p, "conda": c, "required": r, "purpose": d, "version": None,
             "status": "MISSING" if r else "OPTIONAL"} for p, c, r, d in R_PACKAGES]
    out = runner.tool_output([rscript, "--version"])
    rv = parse_version(out, "R")
    code = ("pk <- c(%s); for (p in pk) { v <- tryCatch(as.character(packageVersion(p)), error=function(e) 'NA'); "
            "cat(p, v, '\\n') }" % ",".join(f"'{p}'" for p, *_ in R_PACKAGES))
    res = runner.tool_output([rscript, "-e", code], timeout=180) or ""
    found = {}
    for line in res.splitlines():
        parts = line.split()
        if len(parts) == 2:
            found[parts[0]] = None if parts[1] == "NA" else parts[1]
    pkgs = []
    for p, c, r, d in R_PACKAGES:
        v = found.get(p)
        pkgs.append({"package": p, "conda": c, "required": r, "purpose": d, "version": v,
                     "status": "AVAILABLE" if v else ("MISSING" if r else "OPTIONAL")})
    status = "AVAILABLE" if rv and vtuple(rv) >= (4, 0) else ("OUTDATED" if rv else "INCOMPATIBLE")
    return {"rscript": str(rscript), "version": rv, "status": status}, pkgs


def report(tools, r_info, r_pkgs):
    ui.section("SOFTWARE DEPENDENCIES")
    ui.table([(t["key"], t["status"], t["version"] or "-", t["min_version"], t["purpose"]) for t in tools],
             ["Tool", "Status", "Version", "Min", "Purpose"], max_col=60)
    print()
    ui.table([("R", r_info["status"], r_info["version"] or "-", "4.0", "DESeq2 and downstream analysis")]
             + [(p["package"], p["status"], p["version"] or "-", "", p["purpose"]) for p in r_pkgs],
             ["R component", "Status", "Version", "Min", "Purpose"], max_col=60)
    groups = {}
    for t in tools:
        groups.setdefault(t["status"], []).append(t["key"])
    for p in r_pkgs:
        groups.setdefault(p["status"], []).append("R:" + p["package"])
    if r_info["status"] != "AVAILABLE":
        groups.setdefault(r_info["status"], []).append("R")
    print()
    for s in ("AVAILABLE", "MISSING", "OUTDATED", "INCOMPATIBLE", "OPTIONAL"):
        print(f"  {s:<13} {', '.join(groups.get(s, [])) or '-'}")
    return groups


def blocking(tools, r_info, r_pkgs, need_r=True):
    bad = [t["key"] for t in tools if t["required"] and t["status"] in ("MISSING", "OUTDATED", "INCOMPATIBLE")]
    if need_r:
        if r_info["status"] != "AVAILABLE":
            bad.append("R")
        bad += ["R:" + p["package"] for p in r_pkgs if p["required"] and p["status"] != "AVAILABLE"]
    return bad


def write_versions(path, tools, r_info, r_pkgs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["software\tversion\tpath\tstatus"]
    for t in tools:
        lines.append(f"{t['key']}\t{t['version'] or ''}\t{t['path'] or ''}\t{t['status']}")
    lines.append(f"R\t{r_info.get('version') or ''}\t{r_info.get('rscript') or ''}\t{r_info.get('status')}")
    for p in r_pkgs:
        lines.append(f"R:{p['package']}\t{p['version'] or ''}\t\t{p['status']}")
    path.write_text("\n".join(lines) + "\n")
    return path


def cpu_fix_command(envs, tool):
    """The exact command that replaces a CPU-incompatible build, or None if no compatible build is known."""
    if not tool.get("install_spec"):
        return None
    exe = envs.conda_bin()
    name = Path(str(exe)).name if exe else "mamba"
    return f'{name} install -n {envs.names[tool["env"]]} -c conda-forge -c bioconda "{tool["install_spec"]}"'


def tool_versions_map(tools):
    return {t["key"]: t["version"] for t in tools}


def interactive_install(envs, tools, r_info, r_pkgs, log_file):
    """Explain what is missing, show exact commands, ask permission, install, re-validate."""
    missing_tools = [t for t in tools if t["status"] in ("MISSING", "OUTDATED", "INCOMPATIBLE", "OPTIONAL")]
    missing_r = [p for p in r_pkgs if p["status"] != "AVAILABLE"]
    r_missing = r_info["status"] != "AVAILABLE"
    if not missing_tools and not missing_r and not r_missing:
        ui.ok("All dependencies are available.")
        return False
    exe = envs.conda_bin()
    if not exe:
        ui.error("conda/mamba not found — automatic installation is not possible.")
        ui.info("Install Miniforge (no administrator rights needed):")
        print("    curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/"
              "Miniforge3-Linux-x86_64.sh\n    bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3")
        ui.info("Then restart the pipeline.")
        return False
    plans = []
    if missing_tools:
        if envs.prefix("tools") is None:
            plans.append(("Create tools environment '%s'" % envs.names["tools"], envs.create_commands("tools")))
        else:
            pkgs = sorted({t.get("install_spec") or t["package"] for t in missing_tools})
            plans.append(("Install/upgrade in '%s': %s" % (envs.names["tools"], ", ".join(pkgs)),
                          envs.install_commands("tools", pkgs)))
    if missing_r or r_missing:
        if envs.prefix("r") is None:
            plans.append(("Create R environment '%s'" % envs.names["r"], envs.create_commands("r")))
        else:
            pkgs = sorted({p["conda"] for p in missing_r})
            plans.append(("Install in '%s': %s" % (envs.names["r"], ", ".join(pkgs)),
                          envs.install_commands("r", pkgs)))
    ui.section("PROPOSED INSTALLATION (user space, no sudo)")
    for t in missing_tools:
        print(f"  - {t['key']:<20} {t['purpose']}")
        if t.get("cpu_incompatible"):
            print(f"      {t['note']}.")
            if t.get("install_spec"):
                print(f"      Installs {t['install_spec']}, the newest build that runs on this CPU. The version is "
                      "recorded in the project's software versions.")
    for p in missing_r:
        print(f"  - R:{p['package']:<18} {p['purpose']}")
    print("\nCommands that would be executed:")
    for desc, cmd in plans:
        print(f"  # {desc}\n  {runner.fmt(cmd)}")
    if not ui.ask_yes_no("\nProceed with installation?", default=False):
        ui.skipped("installation declined by user")
        return False
    for desc, cmd in plans:
        ui.running(desc)
        envs.run_install(cmd, log_file)
    envs.activate()
    ui.ok("Installation finished; re-validating...")
    return True
