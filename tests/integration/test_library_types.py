"""F13: library types end-to-end through the real CLI, up to gene counting.

For each synthetic library we know the true strandedness and layout. The pipeline must infer strandedness
correctly and featureCounts must then assign nearly every read — only possible if `-s` was right.
A deliberately wrong strandedness setting must trigger the disagreement prompt and the low-assignment warning.
"""
import json
import os
import subprocess
import sys

import pytest

from conftest import ROOT, have

sys.path.insert(0, str(ROOT / "tests" / "data"))
import make_synthetic  # noqa: E402

pytestmark = pytest.mark.skipif(not have("fastqc", "hisat2", "samtools", "featureCounts", "stringtie", "multiqc",
                                         "infer_experiment.py"), reason="bioinformatics tools not installed")


def run_until_counts(tmp_path, strand, paired, configured="auto", strand_answer=("y",)):
    data = tmp_path / "data"
    make_synthetic.make(data, pairs_per_sample=20000, strand=strand, paired=paired, adapters=False)
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(f"reference_store: {tmp_path / 'store'}\ntrim_adapters: never\nstrandedness: {configured}\n")
    answers = ["1", "Lib", str(tmp_path / "projects"), "5", str(data / "fastq"), "", "1", "y",
               "5", str(data / "genome.fa"), str(data / "annotation.gtf"), "Synthetic", "SYNTH1", "", "v1", "y",
               "n", "1", "y",
               "1", "1", "1", "1",      # data, fastq verification, raw QC, quality gate (trimming decided by config)
               "y",                     # continue with samples flagged REVIEW (<1M reads: expected for test data)
               "1", "1", "1", "1",      # trimming (skipped), reference, alignment, BAM QC
               "1", *strand_answer,     # strandedness
               "1", "1"]                # StringTie2, featureCounts; input then ends -> stops at count matrix
    r = subprocess.run([str(ROOT / "rnaseq_pipeline"), "--config", str(cfg), "--projects-dir",
                        str(tmp_path / "projects")], input="\n".join(answers) + "\n", capture_output=True,
                       text=True, timeout=1800, env=dict(os.environ, PYTHONNOUSERSITE="1"))
    p = tmp_path / "projects" / "RNAseq_Lib"
    assert (p / "checkpoints" / "featurecounts_completed.json").exists(), r.stdout[-3000:]
    state = json.loads((p / "project.json").read_text())
    stats = json.loads((p / "featurecounts" / "featurecounts_stats.json").read_text())
    return state, stats, r.stdout, p


@pytest.mark.parametrize("strand,paired,flag", [("unstranded", True, "0"), ("forward", True, "1"),
                                                ("reverse", False, "2")])
def test_library_type_inferred_and_counted(tmp_path, strand, paired, flag):
    state, stats, out, p = run_until_counts(tmp_path, strand, paired)
    s = state["strandedness"]
    assert s["value"] == strand and s["featurecounts_flag"] == flag and s["source"].startswith("RSeQC")
    assert state["read_type"] == ("paired" if paired else "single")
    for sid, v in stats.items():
        assert v["assigned_pct"] > 95, f"{sid}: only {v['assigned_pct']}% assigned with -s {flag}"
    import csv
    for row in csv.DictReader(open(p / "alignment" / "reports" / "alignment_summary.tsv"), delimiter="\t"):
        assert int(row["primary"]) == int(row["input_reads"]) * (2 if paired else 1)


def test_wrong_strandedness_is_caught(tmp_path):
    """Config says 'forward' for reverse-stranded data; the user keeps the configured value on purpose."""
    state, stats, out, _ = run_until_counts(tmp_path, "reverse", True, configured="forward",
                                            strand_answer=("1",))
    assert "disagrees with inferred 'reverse'" in out
    assert state["strandedness"]["value"] == "forward"
    assert all(v["assigned_pct"] < 10 for v in stats.values())
    assert "low assignment rate" in out
