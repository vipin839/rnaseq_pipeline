"""Reference database manager: selection, download, validation, compatibility checks, HISAT2 index."""
import gzip
import hashlib
import os
import re
import shutil
import signal
from datetime import datetime
from pathlib import Path

from . import PipelineError, UserAbort, net, runner, ui
from . import config as C

IUPAC = b"ACGTNacgtnRYKMSWBDHVrykmswbdhv-*"
GTF_STRANDS = {"+", "-", "."}


# ------------------------------ checksums ------------------------------

def bsd_sum(path):
    """BSD 16-bit checksum as printed by `sum` (used in Ensembl CHECKSUMS)."""
    s = 0
    size = 0
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            size += len(chunk)
            for b in chunk:
                s = ((s >> 1) + ((s & 1) << 15) + b) & 0xFFFF
    return s, (size + 1023) // 1024


def bsd_sum_fast(path):
    """Same as bsd_sum but uses the coreutils `sum` binary when available (much faster)."""
    exe = shutil.which("sum")
    if exe:
        out = runner.tool_output([exe, str(path)], timeout=3600)
        if out:
            parts = out.split()
            return int(parts[0]), int(parts[1])
    return bsd_sum(path)


def expected_checksum(pkg, part, filename):
    """Return ('md5', hex) or ('bsd_sum', (sum, blocks)) or None."""
    ck = pkg.get("checksum", {})
    url = pkg[part].get("checksum_url") or ck.get("url")
    if not url:
        return None
    try:
        text = net.get_text(url)
    except PipelineError:
        ui.warn(f"could not fetch checksum list {url}; file integrity will rely on gzip test only")
        return None
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        name = parts[-1].lstrip("./")
        if name.endswith(filename):
            if ck.get("type") == "bsd_sum" and len(parts) >= 3:
                return ("bsd_sum", (int(parts[0]), int(parts[1])))
            if ck.get("type") == "md5":
                return ("md5", parts[0])
    return None


def verify_checksum(path, expected):
    if expected is None:
        return None
    kind, val = expected
    got = net.md5(path) if kind == "md5" else bsd_sum_fast(path)
    if kind == "bsd_sum":
        return got[0] == val[0]
    return got == val


# ------------------------------ validation ------------------------------

def validate_fasta(path, fai_path=None):
    """Stream FASTA; return {names: [..], lengths: {..}, total_bp}. Raises PipelineError on format problems."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise PipelineError(f"genome FASTA missing or empty: {path}", stage="reference")
    names, lengths, seen = [], {}, set()
    cur, n = None, 0
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rb") as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip(b"\r\n")
            if line.startswith(b">"):
                if cur is not None:
                    lengths[cur] = n
                name = line[1:].split(None, 1)[0].decode(errors="replace") if len(line) > 1 else ""
                if not name:
                    raise PipelineError(f"FASTA header without a name at line {lineno}", stage="reference")
                if name in seen:
                    raise PipelineError(f"duplicate sequence name in FASTA: {name}", stage="reference")
                seen.add(name)
                names.append(name)
                cur, n = name, 0
            else:
                if cur is None:
                    if line.strip():
                        raise PipelineError("FASTA does not start with a '>' header", stage="reference",
                                            cause=f"{path} may not be a FASTA file")
                    continue
                bad = line.translate(None, IUPAC)
                if bad:
                    raise PipelineError(f"invalid characters in FASTA sequence at line {lineno}",
                                        stage="reference", cause=repr(bytes(sorted(set(bad))[:5])))
                n += len(line)
    if cur is not None:
        lengths[cur] = n
    if not names:
        raise PipelineError("FASTA contains no sequences", stage="reference")
    zero = [k for k, v in lengths.items() if v == 0]
    if zero:
        raise PipelineError(f"FASTA has empty sequences: {zero[:5]}", stage="reference")
    return {"names": names, "lengths": lengths, "total_bp": sum(lengths.values()), "n_sequences": len(names)}


ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def validate_gtf(path, attribute="gene_id", feature_type="exon"):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise PipelineError(f"annotation GTF missing or empty: {path}", stage="reference")
    opener = gzip.open if path.name.endswith(".gz") else open
    seq_max_end, genes, transcripts = {}, set(), set()
    exons, features, header = 0, {}, []
    missing_attr = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for lineno, line in enumerate(f, 1):
            if line.startswith("#"):
                if len(header) < 20:
                    header.append(line.strip())
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise PipelineError(f"GTF line {lineno} has {len(cols)} columns (expected 9)", stage="reference",
                                    cause="file is not GTF (GFF3 uses key=value attributes and is not supported "
                                          "directly — convert with gffread -T)")
            seq, _, ftype, start, end, _, strand, _, attrs = cols
            try:
                s, e = int(start), int(end)
            except ValueError:
                raise PipelineError(f"GTF line {lineno}: non-integer coordinates", stage="reference") from None
            if s < 1 or e < s:
                raise PipelineError(f"GTF line {lineno}: invalid interval {s}-{e}", stage="reference")
            if strand not in GTF_STRANDS:
                raise PipelineError(f"GTF line {lineno}: invalid strand {strand!r}", stage="reference")
            if "=" in attrs and '"' not in attrs:
                raise PipelineError("annotation looks like GFF3, not GTF", stage="reference",
                                    remedy="convert with: gffread annotation.gff3 -T -o annotation.gtf")
            features[ftype] = features.get(ftype, 0) + 1
            if e > seq_max_end.get(seq, 0):
                seq_max_end[seq] = e
            if ftype == feature_type:
                exons += 1
                a = dict(ATTR_RE.findall(attrs))
                if attribute not in a:
                    missing_attr += 1
                    continue
                genes.add(a[attribute])
                if "transcript_id" in a:
                    transcripts.add(a["transcript_id"])
    if exons == 0:
        raise PipelineError(f"GTF contains no '{feature_type}' features", stage="reference",
                            cause="featureCounts counts this feature type "
                                  "(bacterial NCBI GTFs annotate genes as CDS, not exon)",
                            remedy=f"set featurecounts_parameters.feature_type to a type present in the GTF "
                                   f"(found: {', '.join(sorted(features)) or 'none'})")
    if missing_attr > exons * 0.01:
        raise PipelineError(f"{missing_attr} {feature_type} lines lack {attribute}/transcript_id", stage="reference",
                            remedy="choose a standard GTF (GENCODE/Ensembl/UCSC)")
    build = None
    for h in header:
        m = re.search(r"(GRC[hm]\d+|hg\d+|mm\d+|R64[-\d]*)", h)
        if m and ("genome-build" in h or "description" in h or "assembly" in h.lower()):
            build = m.group(1)
            break
    return {"seqnames": sorted(seq_max_end), "seq_max_end": seq_max_end, "genes": len(genes),
            "transcripts": len(transcripts), "exons": exons, "feature_types": features,
            "header_build": build, "attribute": attribute, "feature_type": feature_type}


def check_compatibility(fa, gtf, genome_assembly=None, annotation_assembly=None):
    """Return (errors, warnings)."""
    errors, warns = [], []
    fa_names = set(fa["names"])
    gtf_names = set(gtf["seqnames"])
    shared = fa_names & gtf_names
    if not shared:
        hint = ""
        if any(n.startswith("chr") for n in fa_names) != any(n.startswith("chr") for n in gtf_names):
            hint = " (one file uses 'chr1' style names and the other '1' style — different providers)"
        errors.append("no chromosome names are shared between genome FASTA and GTF" + hint)
    else:
        missing = gtf_names - fa_names
        if missing:
            frac = len(missing) / len(gtf_names)
            msg = f"{len(missing)} GTF sequence(s) not in FASTA (e.g. {sorted(missing)[:3]})"
            (errors if frac > 0.05 else warns).append(msg)
        over = [s for s in shared if gtf["seq_max_end"][s] > fa["lengths"][s]]
        if over:
            errors.append(f"annotation coordinates exceed chromosome length on {len(over)} sequence(s) "
                          f"(e.g. {over[0]}: GTF end {gtf['seq_max_end'][over[0]]:,} > length "
                          f"{fa['lengths'][over[0]]:,}) — genome and annotation are from DIFFERENT assemblies")
    if genome_assembly and annotation_assembly and _norm(genome_assembly) != _norm(annotation_assembly):
        errors.append(f"declared assemblies differ: genome {genome_assembly} vs annotation {annotation_assembly}")
    hb = gtf.get("header_build")
    if hb and genome_assembly and _norm(hb) != _norm(genome_assembly):
        errors.append(f"GTF header says build {hb}, but genome assembly is {genome_assembly}")
    return errors, warns


def _norm(a):
    a = str(a).lower().replace(" ", "")
    aliases = {"hg38": "grch38", "hg19": "grch37", "mm39": "grcm39", "mm10": "grcm38"}
    a = re.sub(r"\.p\d+$", "", a)
    return aliases.get(a, a)


# ------------------------------ store / paths ------------------------------

class RefPaths:
    def __init__(self, base):
        self.base = Path(base)
        self.genome = self.base / "genome" / "genome.fa"
        self.gtf = self.base / "annotation" / "annotation.gtf"
        self.transcripts = self.base / "transcriptome" / "transcripts.fa.gz"
        self.splice_sites = self.base / "annotation" / "splice_sites.txt"
        self.exons = self.base / "annotation" / "exons.txt"
        self.bed12 = self.base / "annotation" / "annotation.bed12"
        self.index_dir = self.base / "index"
        self.index_prefix = self.index_dir / "genome"
        self.checksums = self.base / "checksums"
        self.logs = self.base / "logs"
        self.manifest = self.base / "reference_manifest.yaml"


def store_dir(cfg, ref):
    root = Path(os.path.expanduser(cfg.get("reference_store", "~/rnaseq_references")))
    return root / ref["id"]


def custom_id(fasta, gtf):
    h = hashlib.sha1(f"{Path(fasta).resolve()}|{Path(gtf).resolve()}".encode()).hexdigest()[:8]
    return f"custom_{h}"


def packages_for(catalog, organism):
    return {k: v for k, v in catalog["packages"].items() if v["organism"] == organism}


def describe(ref):
    return f"{ref.get('organism')} / {ref.get('assembly')} / {ref.get('annotation_release')} ({ref.get('source')})"


# ------------------------------ acquisition ------------------------------

def _gunzip_to(src, dest, log_file):
    tmp = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    runner.run(["gzip", "-dc", src], stage="reference", stdout_file=tmp, log_file=log_file,
               description=f"Decompressing {Path(src).name}")
    if not runner.DRY_RUN:
        os.replace(tmp, dest)


def acquire(ref, paths, log_file):
    """Download (catalog) or link (custom) genome + GTF into the store."""
    for d in (paths.genome.parent, paths.gtf.parent, paths.checksums, paths.logs, paths.transcripts.parent):
        d.mkdir(parents=True, exist_ok=True)
    record = {}
    if ref["kind"] == "catalog":
        pkg = ref["package"]
        for part, dest in (("genome", paths.genome), ("annotation", paths.gtf)):
            url = pkg[part]["url"]
            gz = paths.checksums.parent / "downloads" / Path(url).name
            if dest.exists():
                ui.skipped(f"{part}: already prepared ({dest})")
                record[part] = {"url": url, "path": str(dest)}
                continue
            exp = expected_checksum(pkg, part, Path(url).name)
            net.download(url, gz, stage="reference", log_file=log_file)
            if not runner.DRY_RUN:
                ok = verify_checksum(gz, exp)
                if ok is False:
                    bad = gz.with_name(gz.name + ".badchecksum")
                    os.replace(gz, bad)
                    raise PipelineError(f"{part} download failed checksum verification", stage="reference",
                                        cause=f"file kept as {bad}", remedy="re-run reference preparation")
                (paths.checksums / f"{part}.checksum.txt").write_text(
                    f"file: {gz.name}\nprovider_checksum: {exp}\nverified: {ok}\nsha256: {net_sha256(gz)}\n")
                ui.ok(f"{part} checksum " + ("verified against provider" if ok else "not published; SHA-256 recorded"))
            _gunzip_to(gz, dest, log_file)
            record[part] = {"url": url, "path": str(dest), "checksum": str(exp)}
        if pkg.get("transcripts"):
            record["transcripts"] = {"url": pkg["transcripts"]["url"], "downloaded": False,
                                     "note": "optional; not required for HISAT2/featureCounts in v1"}
    else:
        for part, src, dest in (("genome", ref["fasta"], paths.genome), ("annotation", ref["gtf"], paths.gtf)):
            src = Path(src)
            if dest.exists() or dest.is_symlink():
                record[part] = {"path": str(dest), "source_file": str(src)}
                continue
            if src.name.endswith(".gz"):
                _gunzip_to(src, dest, log_file)
            else:
                dest.symlink_to(src.resolve())
            record[part] = {"path": str(dest), "source_file": str(src)}
            if not runner.DRY_RUN:
                (paths.checksums / f"{part}.checksum.txt").write_text(f"sha256: {net_sha256(src)}\n")
        if ref.get("external_index") and not index_files(paths.index_prefix):
            paths.index_dir.mkdir(parents=True, exist_ok=True)
            for f in index_files(ref["external_index"]) or []:
                suffix = f.name[len(Path(ref["external_index"]).name):]
                (paths.index_dir / f"genome{suffix}").symlink_to(f.resolve())
            record["index"] = {"external_prefix": ref["external_index"]}
    return record


def net_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------ preparation ------------------------------

def index_files(prefix):
    prefix = Path(prefix)
    small = [prefix.parent / f"{prefix.name}.{i}.ht2" for i in range(1, 9)]
    large = [prefix.parent / f"{prefix.name}.{i}.ht2l" for i in range(1, 9)]
    if all(p.exists() for p in small):
        return small
    if all(p.exists() for p in large):
        return large
    return None


def index_valid(prefix, fasta_names):
    files = index_files(prefix)
    if not files or any(f.stat().st_size == 0 for f in files):
        return False, "index files missing or empty"
    out = runner.tool_output(["hisat2-inspect", "-n", str(prefix)], timeout=600)
    if out is None:
        return False, "hisat2-inspect failed"
    names = {l.split()[0] for l in out.splitlines() if l.strip()}
    if fasta_names is not None and names != set(fasta_names):
        return False, f"index sequences ({len(names)}) differ from genome FASTA ({len(fasta_names)})"
    return True, f"{len(names)} sequences"


BUILD_RETRY_SIGNALS = (signal.SIGSEGV, signal.SIGBUS, signal.SIGABRT)


def build_index(paths, fa_info, threads, mem_gb, use_ss, log_file, extra_args=()):
    """Build HISAT2 index into a temp dir, validate, then move into place."""
    runner.run(["hisat2_extract_splice_sites.py", paths.gtf], stage="reference", stdout_file=paths.splice_sites,
               log_file=log_file, description="Extracting known splice sites from GTF")
    runner.run(["hisat2_extract_exons.py", paths.gtf], stage="reference", stdout_file=paths.exons,
               log_file=log_file, description="Extracting exons from GTF")
    if not runner.DRY_RUN and paths.splice_sites.stat().st_size == 0:
        ui.warn("no splice sites in this annotation (normal for bacteria / single-exon genomes)")
        if use_ss:
            # hisat2-build aborts on an empty --ss file; a plain index is correct when there are no introns
            ui.info("building a plain HISAT2 index (no --ss/--exon)")
            use_ss = False
    tmp = paths.index_dir.with_name("index.building")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    cmd = ["hisat2-build", "-p", str(threads)]
    if use_ss:
        cmd += ["--ss", paths.splice_sites, "--exon", paths.exons]
    cmd += [str(x) for x in extra_args]
    cmd += [paths.genome, tmp / "genome"]
    res = runner.run(cmd, stage="reference", log_file=log_file, check=False,
                     description="Building HISAT2 index (this can take a long time for large genomes)")
    if runner.DRY_RUN:
        return
    sig = runner.signal_of(res.returncode)
    if sig in BUILD_RETRY_SIGNALS and threads > 1:
        # hisat2-build 2.2.x occasionally crashes in its multithreaded phase (~0.6% of builds here; 0 of 395
        # single-threaded builds crashed) and the index it builds does not depend on the thread count, so one
        # single-threaded retry is safe. Out-of-memory kills and SIGILL are not retried: they would recur.
        ui.warn(f"hisat2-build crashed ({signal.Signals(sig).name}) — an intermittent crash in HISAT2's multithreaded "
                "index build, not a problem with your files. The partial index was discarded; retrying once with a "
                "single thread (slower; the index is identical)")
        shutil.rmtree(tmp)
        tmp.mkdir(parents=True)
        cmd[cmd.index("-p") + 1] = "1"
        res = runner.run(cmd, stage="reference", log_file=log_file, check=False,
                         description="Building HISAT2 index again (single thread)")
    if res.returncode != 0:
        meaning = runner.describe_exit(res.returncode)
        raise PipelineError(f"hisat2-build failed (exit code {res.returncode})" + (f" — {meaning}" if meaning else ""),
                            stage="reference", cause=(runner.tail_file(log_file) + f"\n  Full log: {log_file}") if log_file else None,
                            remedy=("run 'rnaseq-pipeline --check': it names a build of the program that runs on "
                                    "this CPU" if "SIGILL" in meaning else
                                    "free memory (see the RAM estimate above) or disk space, then resume; the "
                                    "partial index is never used"))
    ok, why = index_valid(tmp / "genome", fa_info["names"])
    if not ok:
        raise PipelineError(f"built HISAT2 index failed validation: {why}", stage="reference",
                            remedy=f"see {log_file}; check RAM/disk")
    if paths.index_dir.exists():
        old = paths.index_dir.with_name(f"index.old.{datetime.now():%Y%m%d%H%M%S}")
        os.replace(paths.index_dir, old)
        ui.info(f"previous index moved aside to {old}")
    os.replace(tmp, paths.index_dir)


def index_ram_needed_gb(genome_bp, use_ss):
    return genome_bp * (60 if use_ss else 2.8) / 1e9 + 0.5


def gtf_to_bed12(gtf, bed, feature_type="exon"):
    """Convert GTF features to BED12 (for RSeQC infer_experiment). For annotations without exons
    (bacteria) the counting feature (CDS) is used instead, keyed by gene."""
    tx = {}
    with open(gtf, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) != 9 or c[2] != feature_type:
                continue
            m = re.search(r'transcript_id "([^"]+)"', c[8]) or re.search(r'gene_id "([^"]+)"', c[8])
            if not m:
                continue
            t = tx.setdefault(m.group(1), [c[0], c[6], []])
            t[2].append((int(c[3]) - 1, int(c[4])))
    tmp = Path(str(bed) + ".part")
    with open(tmp, "w") as out:
        for tid, (chrom, strand, ex) in tx.items():
            ex.sort()
            start, end = ex[0][0], max(e for _, e in ex)
            sizes = ",".join(str(e - s) for s, e in ex)
            starts = ",".join(str(s - start) for s, _ in ex)
            out.write(f"{chrom}\t{start}\t{end}\t{tid}\t0\t{strand if strand in '+-' else '+'}\t{start}\t{end}"
                      f"\t0\t{len(ex)}\t{sizes},\t{starts},\n")
    os.replace(tmp, bed)
    return len(tx)


def write_manifest(path, data):
    C.save_yaml(data, path)


def link_into_project(project, paths):
    """Symlink store files into <project>/reference so the project is self-describing (no copies)."""
    pairs = [(paths.genome, project.path("reference", "genome", "genome.fa")),
             (paths.gtf, project.path("reference", "annotation", "annotation.gtf")),
             (paths.splice_sites, project.path("reference", "annotation", "splice_sites.txt")),
             (paths.exons, project.path("reference", "annotation", "exons.txt")),
             (paths.bed12, project.path("reference", "annotation", "annotation.bed12"))]
    for f in index_files(paths.index_prefix) or []:
        pairs.append((f, project.path("reference", "index", f.name)))
    for src, dst in pairs:
        if not Path(src).exists():
            continue
        if dst.is_symlink() or dst.exists():
            if dst.resolve() == Path(src).resolve():
                continue
            if not dst.is_symlink():
                raise PipelineError(f"{dst} is a real file, not a link; refusing to replace it",
                                    remedy="move it aside manually if it is no longer needed")
            dst.unlink()
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.symlink_to(Path(src).resolve())


# ------------------------------ NCBI assembly search (any organism) ------------------------------

DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome"
NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


def ftp_dir(accession, assembly_name):
    """https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/265/GCF_000013265.1_ASM1326v1"""
    prefix, digits = accession.split("_", 1)
    num = digits.split(".")[0]
    parts = [num[i:i + 3] for i in range(0, 9, 3)]
    safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", assembly_name)
    return f"{NCBI_FTP}/{prefix}/{parts[0]}/{parts[1]}/{parts[2]}/{accession}_{safe_name}"


def search_ncbi(query, limit=10, reference_only=True):
    """Search NCBI Datasets for assemblies of an organism (name or taxid). Returns display records."""
    import json as _json
    import urllib.parse
    q = urllib.parse.quote(str(query).strip())
    url = (f"{DATASETS_API}/taxon/{q}/dataset_report?page_size={limit}"
           f"&filters.assembly_version=current" + ("&filters.reference_only=true" if reference_only else ""))
    text = net.get_text(url)
    try:
        data = _json.loads(text)
    except ValueError:
        raise PipelineError(f"unexpected response from NCBI Datasets for {query!r}") from None
    out = []
    for r in data.get("reports", []):
        info = r.get("assembly_info", {})
        ann = r.get("annotation_info") or {}
        stats = r.get("assembly_stats", {}) or {}
        acc, name = r.get("accession", ""), info.get("assembly_name", "")
        if not acc or not name:
            continue
        out.append({
            "accession": acc, "assembly_name": name,
            "organism": r.get("organism", {}).get("organism_name", ""),
            "strain": (r.get("organism", {}).get("infraspecific_names") or {}).get("strain", ""),
            "level": info.get("assembly_level", ""), "annotation": ann.get("name", ""),
            "annotation_date": ann.get("release_date", ""),
            "genome_bp": int(stats.get("total_sequence_length") or 0),
            "refseq_category": info.get("refseq_category", ""),
            "ftp": ftp_dir(acc, name),
        })
    return out


def ncbi_package(rec):
    """Turn a search result into a catalog-style package (same shape as config/reference_catalog.yaml)."""
    ftp = rec.get("ftp") or ftp_dir(rec["accession"], rec["assembly_name"])
    base = f"{ftp}/{rec['accession']}_{re.sub(r'[^A-Za-z0-9._-]', '_', rec['assembly_name'])}"
    gb = max(0.001, rec["genome_bp"] / 1e9 * 0.3)  # gzipped FASTA is roughly 0.3 bytes per base
    return {
        "organism": rec["organism"], "source": "refseq" if rec["accession"].startswith("GCF") else "genbank",
        "assembly": rec["assembly_name"], "annotation_release": rec["annotation"] or rec["accession"],
        "genome": {"url": f"{base}_genomic.fna.gz", "approx_size_gb": round(gb, 3)},
        "annotation": {"url": f"{base}_genomic.gtf.gz", "approx_size_gb": round(gb * 0.06, 3)},
        "checksum": {"type": "md5", "url": f"{ftp}/md5checksums.txt"},
        "index_approx_size_gb": round(max(0.02, rec["genome_bp"] / 1e9 * 1.5), 2),
        "genome_size_bp": rec["genome_bp"],
    }


def annotation_available(rec):
    """HEAD the GTF: not every assembly has one."""
    pkg = ncbi_package(rec)
    out = runner.tool_output(["curl", "-sIL", "-o", "/dev/null", "-w", "%{http_code}", "--max-time", "25",
                              pkg["annotation"]["url"]], timeout=60)
    return (out or "").strip().endswith("200")


def feature_type_counts(path):
    """Fast scan of GTF column 3 -> {feature_type: count} (used to detect bacterial CDS-only annotation)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    counts = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t", 4)
            if len(cols) > 2:
                counts[cols[2]] = counts.get(cols[2], 0) + 1
    return counts


def suggest_feature_type(counts, configured):
    """Return (suggested_type, reason) when the configured type is clearly wrong for this annotation."""
    have = counts.get(configured, 0)
    genes = counts.get("gene", 0)
    # a usable counting feature occurs at least about once per gene (eukaryotes: many exons per gene)
    if have and (not genes or have >= 0.9 * genes):
        return None, None
    for alt in ("CDS", "exon", "transcript"):
        if alt != configured and counts.get(alt, 0) >= max(1, 0.5 * genes):
            return alt, (f"annotation has {have:,} '{configured}' but {counts[alt]:,} '{alt}' features "
                         f"for {genes:,} genes (typical of bacterial/NCBI annotation)")
    return None, None


# ------------------------------ full preparation (used with or without a project) ------------------------------

def prepare(ref, cfg, threads, mem_gb, ram_available_gb, log_file, on_feature_type=None, confirm_ram=None,
            dry_run=False):
    """Acquire, validate and index a reference in the shared store.

    on_feature_type(alt, why, counts) -> bool   asked when the annotation has no usable counting feature
    confirm_ram(need_gb, avail_gb)     -> bool   asked when the index build may not fit in memory
    Returns (manifest, info) where info carries what a project needs (index prefix, gtf, bed12, ...).
    """
    paths = RefPaths(ref["store"])
    record = acquire(ref, paths, log_file)
    if dry_run:
        build_index(paths, {"names": []}, threads, mem_gb, False, log_file)
        return None, None

    ui.running("Validating genome FASTA (streaming)")
    fa = validate_fasta(paths.genome)
    ui.ok(f"FASTA valid: {fa['n_sequences']} sequences, {fa['total_bp']:,} bp")

    ui.running("Validating GTF")
    ftype = cfg["featurecounts_parameters"]["feature_type"]
    counts = feature_type_counts(paths.gtf)
    alt, why = suggest_feature_type(counts, ftype)
    switched = False
    if alt:
        ui.warn(why)
        ui.info(f"counting '{ftype}' would miss most genes in this annotation")
        if on_feature_type is None or on_feature_type(alt, why, counts):
            ftype, switched = alt, True
    gtf = validate_gtf(paths.gtf, cfg["featurecounts_parameters"]["attribute"], ftype)
    ui.ok(f"GTF valid: {gtf['genes']:,} genes, {gtf['transcripts']:,} transcripts, {gtf['exons']:,} {ftype} features")
    if "exon" not in gtf["feature_types"]:
        ui.info("no 'exon' features (typical for bacterial annotation): HISAT2 will align without splice sites "
                "and StringTie2 transcript quantification is not meaningful")

    errors, warns = check_compatibility(fa, gtf, ref.get("assembly"), ref.get("annotation_assembly",
                                                                              ref.get("assembly")))
    for w in warns:
        ui.warn(w)
    if errors:
        for e in errors:
            ui.error(e)
        if not ref.get("override_compatibility"):
            raise PipelineError("genome and annotation are not compatible", stage="reference",
                                remedy="select a matching genome/annotation pair (same assembly and provider)")
        ui.warn("compatibility errors OVERRIDDEN by user (recorded in manifest)")
    ui.ok("genome/annotation compatibility checks passed")

    runner.run(["samtools", "faidx", paths.genome], stage="reference", log_file=log_file,
               description="Indexing FASTA")
    bed_feature = ftype if counts.get("exon", 0) < counts.get("gene", 0) * 0.9 else "exon"
    n_tx = gtf_to_bed12(paths.gtf, paths.bed12, bed_feature)
    ui.ok(f"BED12 written for strandedness inference ({n_tx:,} {bed_feature} features)")

    setting = cfg["hisat2_build"]["use_splice_sites_in_index"]
    use_ss = (index_ram_needed_gb(fa["total_bp"], True) <= ram_available_gb * 0.9) if setting == "auto" \
        else bool(setting)
    ok, why_idx = index_valid(paths.index_prefix, fa["names"])
    previous = C.load_yaml(paths.manifest) if paths.manifest.exists() else {}
    if ok and paths.splice_sites.exists():
        use_ss = bool(previous.get("index", {}).get("splice_sites_in_index", False))
        ui.skipped(f"existing HISAT2 index is valid ({why_idx}); not rebuilding")
    else:
        need = index_ram_needed_gb(fa["total_bp"], use_ss)
        ui.info(f"HISAT2 index build: splice sites in index = {use_ss} (est. RAM {need:.1f} GB; "
                f"available {ram_available_gb} GB)")
        if not use_ss:
            ui.info("known splice sites will be supplied at alignment time (--known-splicesite-infile)")
        if need > ram_available_gb:
            ui.warn("estimated RAM exceeds available memory; the build may fail or swap heavily")
            if confirm_ram is not None and not confirm_ram(need, ram_available_gb):
                raise UserAbort("index build cancelled")
        build_index(paths, fa, threads, mem_gb, use_ss, log_file, cfg["hisat2_build"].get("extra_args", []))
        ok, why_idx = index_valid(paths.index_prefix, fa["names"])
        if not ok:
            raise PipelineError(f"HISAT2 index invalid after build: {why_idx}", stage="reference")
        ui.ok(f"HISAT2 index built and validated ({why_idx})")

    ui.running("Computing SHA-256 of genome and annotation (once; recorded for reproducibility)")
    checksums = {}
    prev = (previous or {}).get("sha256") or {}
    for label, fpath in (("genome", paths.genome), ("annotation", paths.gtf)):
        st = Path(fpath).stat()
        old = prev.get(label) or {}
        if old.get("size") == st.st_size and old.get("mtime_ns") == st.st_mtime_ns:
            checksums[label] = old
        else:
            checksums[label] = {"sha256": net_sha256(fpath), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    manifest = {
        "sha256": checksums,
        "organism": ref.get("organism"), "source": ref.get("source"), "genome_assembly": ref.get("assembly"),
        "annotation_assembly": ref.get("annotation_assembly", ref.get("assembly")),
        "annotation_release": ref.get("annotation_release"), "label": ref.get("label"), "files": record,
        "genome": {"path": str(paths.genome), "sequences": fa["n_sequences"], "total_bp": fa["total_bp"]},
        "annotation": {"path": str(paths.gtf), "genes": gtf["genes"], "transcripts": gtf["transcripts"],
                       "counted_features": gtf["exons"], "feature_type": ftype,
                       "feature_types": gtf["feature_types"], "header_build": gtf["header_build"]},
        "index": {"prefix": str(paths.index_prefix), "status": "valid", "splice_sites_in_index": use_ss,
                  "splice_sites_file": str(paths.splice_sites), "exons_file": str(paths.exons)},
        "compatibility": {"errors": errors, "warnings": warns, "override": bool(ref.get("override_compatibility"))},
        "store": str(paths.base),
    }
    write_manifest(paths.manifest, manifest)
    info = {"paths": paths, "fa": fa, "gtf": gtf, "feature_type": ftype, "feature_type_switched": switched,
            "counts": counts, "use_ss": use_ss,
            "state": {"index_prefix": str(paths.index_prefix), "index_has_splice_sites": use_ss,
                      "no_splice_sites": paths.splice_sites.exists() and paths.splice_sites.stat().st_size == 0,
                      "splice_sites": str(paths.splice_sites), "gtf": str(paths.gtf), "bed12": str(paths.bed12),
                      "genes": gtf["genes"], "genome_bp": fa["total_bp"]}}
    return manifest, info


def scan_store(store_root):
    """Find usable references under a directory: prepared ones (manifest) and plain FASTA+GTF folders."""
    root = Path(os.path.expanduser(str(store_root)))
    out = []
    if not root.is_dir():
        return out
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        man = d / "reference_manifest.yaml"
        if man.exists():
            m = C.load_yaml(man)
            out.append({"kind": "prepared", "dir": d, "label": m.get("label") or d.name, "manifest": m,
                        "indexed": bool(index_files(d / "index" / "genome"))})
            continue
        fa = sorted(list(d.glob("*.fa")) + list(d.glob("*.fasta")) + list(d.glob("*.fna"))
                    + list(d.glob("*.fa.gz")) + list(d.glob("*.fasta.gz")) + list(d.glob("*.fna.gz")))
        gtf = sorted(list(d.glob("*.gtf")) + list(d.glob("*.gtf.gz")))
        if fa and gtf:
            extra = f", {len(fa)} FASTA / {len(gtf)} GTF candidates" if len(fa) > 1 or len(gtf) > 1 else ""
            out.append({"kind": "files", "dir": d,
                        "label": f"{d.name} ({fa[0].name} + {gtf[0].name}{extra})",
                        "fasta": fa[0], "gtf": gtf[0], "fastas": fa, "gtfs": gtf, "indexed": False})
    return out
