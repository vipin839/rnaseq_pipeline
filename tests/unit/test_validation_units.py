"""Unit tests: path/name/accession validation, config, metadata parsing, count matrix, design, checkpoints,
strandedness parsing, reference compatibility."""
import json
from pathlib import Path

import pytest

from rnaseq import (PipelineError, checkpoint, count_matrix, data_manager, design, geo_manager,
                    reference_manager, strandedness)
from rnaseq import config as C
from rnaseq import validators as V
from rnaseq.project import Project


# ---------------- names / paths / accessions ----------------
@pytest.mark.parametrize("bad", ["", " ", "../x", "a/b", "a b", "-x", "x;rm -rf", "$(id)", "a" * 65, ".."])
def test_project_name_rejects(bad):
    with pytest.raises(ValueError):
        V.project_name(bad)


def test_project_name_ok():
    assert V.project_name(" Liver_study-2.v1 ") == "Liver_study-2.v1"


def test_within_blocks_traversal(tmp_path):
    with pytest.raises(ValueError):
        V.within(tmp_path, "../../etc/passwd")
    assert V.within(tmp_path, "a/b") == (tmp_path / "a" / "b").resolve()


@pytest.mark.parametrize("acc,kind", [("SRR1066657", "run"), ("err123456", "run"), ("SRX403442", "experiment"),
                                      ("PRJNA232717", "study"), ("SRP034844", "study"), ("GSE52778", "geo_series"),
                                      ("SAMN02487085", "sample")])
def test_accessions(acc, kind):
    assert V.accession(acc) == (acc.upper(), kind)


@pytest.mark.parametrize("bad", ["SRR12", "XYZ123456", "GSE", "SRR1066657;ls", ""])
def test_bad_accessions(bad):
    with pytest.raises(ValueError):
        V.accession(bad)


def test_accession_list_dedup():
    assert [a for a, _ in V.accession_list("SRR1066657, SRR1066657 SRR1066658")] == ["SRR1066657", "SRR1066658"]


def test_numeric_validators():
    assert V.probability("0.05") == 0.05
    for bad in ("0", "1", "-0.1", "abc", "nan"):
        with pytest.raises(ValueError):
            V.probability(bad)
    with pytest.raises(ValueError):
        V.threads(0, 8)
    with pytest.raises(ValueError):
        V.threads(99, 8)
    assert V.percentage(50) == 50


# ---------------- config ----------------
def test_default_config_valid():
    assert C.validate(C.load(), 8) == []


@pytest.mark.parametrize("key,val", [("alpha", 1.5), ("alpha", 0), ("log2fc_threshold", -1), ("threads", 0),
                                     ("strandedness", "sideways"), ("design_formula", "condition"),
                                     ("design_formula", "~ a*b"), ("quality_gate.adapter_max_percent_trim", 150),
                                     ("fastp_parameters.extra_args", "--foo"), ("top_n_genes", 0)])
def test_invalid_config_values(key, val):
    cfg = C.load()
    C.set_value(cfg, key, val)
    assert C.validate(cfg, 8), f"{key}={val} should be rejected"


# ---------------- sample detection / metadata parsing ----------------
def test_pairing(tmp_path):
    for n in ["A_R1.fastq.gz", "A_R2.fastq.gz", "B_S2_L001_R1_001.fastq.gz", "B_S2_L001_R2_001.fastq.gz",
              "C_1.fq.gz", "C_2.fq.gz"]:
        (tmp_path / n).write_bytes(b"")
    s = data_manager.pair_files(data_manager.find_fastqs([tmp_path]))
    assert [x["id"] for x in s] == ["A", "B_S2_L001", "C"]
    assert all(x["r2"] is not None for x in s)
    assert s[1]["bio_unit"] == "B"


def test_mixed_layout_rejected(tmp_path):
    for n in ["A_R1.fastq.gz", "A_R2.fastq.gz", "B.fastq.gz"]:
        (tmp_path / n).write_bytes(b"")
    with pytest.raises(PipelineError):
        data_manager.pair_files(data_manager.find_fastqs([tmp_path]))


def test_geo_soft_parsing():
    text = ("^SAMPLE = GSM1\n!Sample_title = ctrl_rep1\n!Sample_characteristics_ch1 = treatment: none\n"
            "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX384345\n"
            "^SAMPLE = GSM2\n!Sample_title = dex_rep1\n!Sample_characteristics_ch1 = treatment: dex\n")
    s = geo_manager.parse_soft(text)
    assert s[0]["gsm"] == "GSM1" and s[0]["srx"] == ["SRX384345"] and s[1]["characteristics"]["treatment"] == "dex"


# ---------------- count matrix ----------------
def _write(path, text):
    path.write_text(text)
    return path


def test_count_matrix_ok(tmp_path):
    f = _write(tmp_path / "m.tsv", "Gene_ID\tA\tB\ng1\t1\t2\ng2\t0\t5\n")
    info = count_matrix.validate_file(f, ["A", "B"], 2)
    assert info["library_sizes"] == {"A": 1, "B": 7}


@pytest.mark.parametrize("body,err", [
    ("Gene_ID\tA\tB\ng1\t1\t-2\n", "negative"), ("Gene_ID\tA\tB\ng1\t1\t2.5\n", "non-integer"),
    ("Gene_ID\tA\tB\ng1\t1\tNA\n", "missing"), ("Gene_ID\tA\tB\ng1\t1\t2\ng1\t3\t4\n", "duplicate gene"),
    ("Gene_ID\tA\tA\ng1\t1\t2\n", "duplicate sample"), ("Gene_ID\tA\tLength\ng1\t1\t2\n", "annotation columns"),
    ("Gene_ID\tA\tB\ng1\t0\t2\n", "zero total"), ("Gene_ID\tA\tB\ng1\t1\n", "wrong number")])
def test_count_matrix_rejects(tmp_path, body, err):
    f = _write(tmp_path / "m.tsv", body)
    with pytest.raises(PipelineError) as e:
        count_matrix.validate_file(f)
    assert err in str(e.value)


def test_sample_matching():
    assert count_matrix.check_metadata_match(["A", "B"], ["A", "B"]) == []
    assert count_matrix.check_metadata_match(["A", "B"], ["B", "A"])  # order
    assert "not in metadata" in count_matrix.check_metadata_match(["A", "B", "C"], ["A", "B"])[0]


# ---------------- design ----------------
def table(rows, cols=("condition",)):
    return {"columns": list(cols), "rows": {s: dict(zip(cols, v)) for s, v in rows.items()}}


def test_design_valid():
    t = table({"a": ("C",), "b": ("C",), "c": ("T",), "d": ("T",)})
    assert design.validate(t, "~ condition", "C", [("T", "C")]) == []


def test_design_requires_replicates():
    t = table({"a": ("C",), "b": ("C",), "c": ("T",)})
    assert any("replicates" in p for p in design.validate(t, "~ condition", "C", [("T", "C")]))


def test_design_confounded_batch():
    t = table({"a": ("b1", "C"), "b": ("b1", "C"), "c": ("b2", "T"), "d": ("b2", "T")}, ("batch", "condition"))
    probs = design.validate(t, "~ batch + condition", "C", [("T", "C")])
    assert any("not full rank" in p for p in probs)


def test_design_batch_ok():
    t = table({"a": ("b1", "C"), "b": ("b2", "C"), "c": ("b1", "T"), "d": ("b2", "T"), "e": ("b1", "C"),
               "f": ("b2", "T")}, ("batch", "condition"))
    assert design.validate(t, "~ batch + condition", "C", [("T", "C")]) == []


@pytest.mark.parametrize("formula", ["condition", "~ condition*batch", "~ missing", "~ condition + condition", ""])
def test_design_invalid_formula(formula):
    t = table({"a": ("C",), "b": ("C",), "c": ("T",), "d": ("T",)})
    assert design.validate(t, formula, "C", [])


def test_design_empty_and_bad_level():
    t = table({"a": ("C",), "b": ("",), "c": ("T",), "d": ("1bad",)})
    probs = design.validate(t, "~ condition", "C", [])
    assert any("empty" in p for p in probs) and any("invalid value" in p for p in probs)


def test_design_import_missing_sample(tmp_path):
    t = design.new_table(["A", "B", "C"])
    f = _write(tmp_path / "meta.tsv", "sample\tcondition\nA\tx\nB\ty\n")
    with pytest.raises(PipelineError, match="missing samples"):
        design.import_tsv(t, f)


# ---------------- checkpoints ----------------
def test_checkpoint_tamper_and_cascade(tmp_path):
    p = Project.create("cp", tmp_path)
    cps = checkpoint.Checkpoints(p)
    a = p.path("data", "a.txt")
    a.write_text("hello")
    cps.write("s1", [a], params={"x": 1})
    b = p.path("data", "b.txt")
    b.write_text("world")
    cps.write("s2", [b], depends_on=["s1"])
    assert cps.verify_files("s1") == [] and cps.verify_files("s2") == []
    fp = cps.fingerprint("s1")
    cps.write("s1", [a], params={"x": 1})  # identical re-run -> same fingerprint, downstream stays valid
    assert cps.fingerprint("s1") == fp and cps.verify_files("s2") == []
    a.write_text("HELLO")  # same size, different content
    assert any("content changed" in x or "modified" in x for x in cps.verify_files("s1", deep=True))
    cps.write("s1", [a], params={"x": 1})
    assert any("upstream" in x for x in cps.verify_files("s2"))
    b.unlink()
    assert any("missing" in x for x in cps.verify_files("s2"))


# ---------------- strandedness ----------------
PE = """This is PairEnd Data
Fraction of reads failed to determine: 0.0100
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0200
Fraction of reads explained by "1+-,1-+,2++,2--": 0.9700
"""
UN = PE.replace("0.0200", "0.4900").replace("0.9700", "0.5000")
AMB = PE.replace("0.0200", "0.3000").replace("0.9700", "0.6900")


def test_strandedness_calls():
    assert strandedness.call(strandedness.parse_infer(PE), 0.8, 0.1, 0.5)[0] == "reverse"
    assert strandedness.call(strandedness.parse_infer(UN), 0.8, 0.1, 0.5)[0] == "unstranded"
    assert strandedness.call(strandedness.parse_infer(AMB), 0.8, 0.1, 0.5)[0] is None
    assert strandedness.call(strandedness.parse_infer(""), 0.8, 0.1, 0.5)[0] is None
    c = {"a": ("reverse", "high", ""), "b": ("forward", "high", "")}
    assert strandedness.consensus(c)[0] is None


# ---------------- reference compatibility ----------------
def _fa(path, seqs):
    with open(path, "w") as f:
        for n, s in seqs.items():
            f.write(f">{n}\n{s}\n")
    return path


def _gtf(path, lines):
    path.write_text("".join(f"{c}\tsrc\texon\t{s}\t{e}\t.\t+\t.\tgene_id \"g{i}\"; transcript_id \"t{i}\";\n"
                            for i, (c, s, e) in enumerate(lines)))
    return path


def test_reference_compatible(tmp_path):
    fa = reference_manager.validate_fasta(_fa(tmp_path / "g.fa", {"chr1": "ACGT" * 100}))
    gtf = reference_manager.validate_gtf(_gtf(tmp_path / "a.gtf", [("chr1", 1, 100)]))
    assert reference_manager.check_compatibility(fa, gtf, "X1", "X1") == ([], [])


def test_reference_chr_style_mismatch(tmp_path):
    fa = reference_manager.validate_fasta(_fa(tmp_path / "g.fa", {"chr1": "ACGT" * 100}))
    gtf = reference_manager.validate_gtf(_gtf(tmp_path / "a.gtf", [("1", 1, 100)]))
    errs, _ = reference_manager.check_compatibility(fa, gtf)
    assert errs and "chr1" in errs[0]


def test_reference_assembly_mismatch_by_coordinates(tmp_path):
    # annotation made for a longer chromosome (different assembly) -> coordinates overflow
    fa = reference_manager.validate_fasta(_fa(tmp_path / "g.fa", {"chr1": "ACGT" * 100}))
    gtf = reference_manager.validate_gtf(_gtf(tmp_path / "a.gtf", [("chr1", 300, 900)]))
    errs, _ = reference_manager.check_compatibility(fa, gtf)
    assert any("DIFFERENT assemblies" in e for e in errs)


def test_reference_declared_assembly_mismatch(tmp_path):
    fa = reference_manager.validate_fasta(_fa(tmp_path / "g.fa", {"chr1": "ACGT" * 100}))
    gtf = reference_manager.validate_gtf(_gtf(tmp_path / "a.gtf", [("chr1", 1, 10)]))
    errs, _ = reference_manager.check_compatibility(fa, gtf, "GRCh38", "GRCh37")
    assert any("declared assemblies differ" in e for e in errs)
    assert reference_manager.check_compatibility(fa, gtf, "GRCh38", "hg38")[0] == []


@pytest.mark.parametrize("content,msg", [("", "missing or empty"), ("ACGT\n", "does not start"),
                                          (">a\nACGT\n>a\nAC\n", "duplicate"), (">a\nAC1T\n", "invalid characters")])
def test_bad_fasta(tmp_path, content, msg):
    f = tmp_path / "x.fa"
    f.write_text(content)
    with pytest.raises(PipelineError, match=msg):
        reference_manager.validate_fasta(f)


def test_gff3_rejected(tmp_path):
    f = tmp_path / "x.gtf"
    f.write_text("chr1\tsrc\texon\t1\t10\t.\t+\t.\tID=e1;Parent=t1\n")
    with pytest.raises(PipelineError, match="GFF3"):
        reference_manager.validate_gtf(f)


def test_missing_reference(tmp_path):
    with pytest.raises(PipelineError, match="missing or empty"):
        reference_manager.validate_fasta(tmp_path / "none.fa")
    with pytest.raises(PipelineError, match="missing or empty"):
        reference_manager.validate_gtf(tmp_path / "none.gtf")


def test_bsd_sum_matches_coreutils(tmp_path):
    f = tmp_path / "x.bin"
    f.write_bytes(bytes(range(256)) * 50)
    assert reference_manager.bsd_sum(f) == reference_manager.bsd_sum_fast(f)
