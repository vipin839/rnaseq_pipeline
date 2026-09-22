
from conftest import write_fastq
from rnaseq_pipeline import fastq_validator as FV


def status(rows):
    return {r["mate"]: r["status"] for r in rows}


def test_valid_pair(good_pair):
    rows = FV.validate_sample("S1", *good_pair)
    assert status(rows) == {"R1": "PASS", "R2": "PASS"}
    assert rows[0]["reads"] == 100 and rows[0]["mean_length"] == 40


def test_single_end_uncompressed(tmp_path):
    f = write_fastq(tmp_path / "a.fastq", [("x", "ACGTN", "IIIII")], gz=False)
    rows = FV.validate_sample("A", f)
    assert rows[0]["status"] == "PASS" and rows[0]["reads"] == 1


def test_quality_length_mismatch(tmp_path):
    f = write_fastq(tmp_path / "a.fq.gz", [("x", "ACGT", "III"), ("y", "ACGT", "IIII")])
    r = FV.validate_sample("A", f)[0]
    assert r["status"] == "FAILED" and r["problem"] == "length mismatch" and r["record"] == 1


def test_invalid_base(tmp_path):
    f = write_fastq(tmp_path / "a.fq.gz", [("x", "ACGU", "IIII")])
    assert FV.validate_sample("A", f)[0]["problem"] == "invalid sequence character"


def test_invalid_quality_char(tmp_path):
    f = write_fastq(tmp_path / "a.fq.gz", [("x", "ACGT", "II I")])
    assert FV.validate_sample("A", f)[0]["problem"] == "invalid quality character"


def test_bad_header_and_separator(tmp_path):
    f = tmp_path / "a.fq"
    f.write_text("x\nACGT\n+\nIIII\n")
    assert FV.validate_sample("A", f)[0]["problem"] == "invalid header"
    f.write_text("@x\nACGT\n-\nIIII\n")
    assert FV.validate_sample("A", f)[0]["problem"] == "invalid separator"


def test_truncated_record(tmp_path):
    f = tmp_path / "a.fq"
    f.write_text("@x\nACGT\n+\nIIII\n@y\nACGT\n")
    assert FV.validate_sample("A", f)[0]["problem"] == "truncated file"


def test_truncated_gzip(tmp_path):
    f = write_fastq(tmp_path / "a.fq.gz", [(f"r{i}", "ACGT" * 20, "I" * 80) for i in range(5000)])
    data = f.read_bytes()
    f.write_bytes(data[: len(data) // 2])
    r = FV.validate_sample("A", f)[0]
    assert r["status"] == "FAILED"
    assert r["problem"] in ("gzip integrity failure", "truncated file", "length mismatch", "invalid header",
                            "invalid separator")


def test_corrupted_gzip_bytes(tmp_path):
    f = write_fastq(tmp_path / "a.fq.gz", [(f"r{i}", "ACGT" * 20, "I" * 80) for i in range(5000)])
    data = bytearray(f.read_bytes())
    for i in range(200, 260):
        data[i] ^= 0xFF
    f.write_bytes(bytes(data))
    assert FV.validate_sample("A", f)[0]["status"] == "FAILED"


def test_pair_count_mismatch(tmp_path):
    a = write_fastq(tmp_path / "R1.fq.gz", [(f"r{i}/1", "ACGT", "IIII") for i in range(10)])
    b = write_fastq(tmp_path / "R2.fq.gz", [(f"r{i}/2", "ACGT", "IIII") for i in range(9)])
    rows = FV.validate_sample("P", a, b)
    assert all(r["status"] == "FAILED" for r in rows)
    assert any(r["problem"] == "paired read count mismatch" for r in rows)


def test_mate_id_mismatch(tmp_path):
    a = write_fastq(tmp_path / "R1.fq.gz", [("r1 1:N", "ACGT", "IIII"), ("r2 1:N", "ACGT", "IIII")])
    b = write_fastq(tmp_path / "R2.fq.gz", [("r1 2:N", "ACGT", "IIII"), ("rX 2:N", "ACGT", "IIII")])
    rows = FV.validate_sample("P", a, b)
    assert any(r["problem"] == "mate ID mismatch" for r in rows)


def test_empty_and_missing(tmp_path):
    e = tmp_path / "e.fq.gz"
    e.write_bytes(b"")
    assert FV.validate_sample("A", e)[0]["problem"] == "empty file"
    assert FV.validate_sample("A", tmp_path / "nope.fq")[0]["problem"] == "file not found"


def test_phred64_warning(tmp_path):
    f = write_fastq(tmp_path / "a.fq.gz", [("x", "ACGT", "hhhh")])
    assert FV.validate_sample("A", f)[0]["status"] == "WARNING"


def test_report_and_parallel(tmp_path, good_pair):
    jobs = [("S1", str(good_pair[0]), str(good_pair[1]), "ACGTNacgtn.", 33)] * 1
    rows = FV.validate_many(jobs, workers=2)
    p = FV.write_report(rows, tmp_path / "rep.tsv")
    lines = p.read_text().splitlines()
    assert lines[0].startswith("sample\tfile") and len(lines) == 3
