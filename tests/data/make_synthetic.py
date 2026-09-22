#!/usr/bin/env python3
"""Generate a small, fully controlled RNA-seq dataset with known truth.

Outputs (in OUT_DIR):
  genome.fa, annotation.gtf
  fastq/{Ctrl1,Ctrl2,Ctrl3,Treat1,Treat2,Treat3}_R{1,2}.fastq.gz   (paired-end, reverse-stranded, dUTP-like)
  truth.tsv  (gene_id, true_log2fc, direction)
Known truth: GENE0001-0008 up 4x in Treat, GENE0009-0016 down 4x, others unchanged.
~15% of fragments are shorter than the read length -> adapter read-through (should trigger trimming).
"""
import gzip
import math
import os
import random
import sys
from pathlib import Path

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def make(out, pairs_per_sample=150000, read_len=75, n_genes=60, seed=7, strand="reverse", paired=True,
         adapters=True):
    """strand: 'reverse' (dUTP-like: read 1 antisense), 'forward' (read 1 sense) or 'unstranded' (random).
    paired=False writes single-end files <sample>.fastq.gz containing read 1 only."""
    assert strand in ("reverse", "forward", "unstranded")
    rng = random.Random(seed)
    out = Path(out)
    (out / "fastq").mkdir(parents=True, exist_ok=True)
    chroms = {"chr1": 320000, "chr2": 240000, "chr3": 160000}
    seqs = {c: "".join(rng.choice("ACGT") for _ in range(L)) for c, L in chroms.items()}
    with open(out / "genome.fa", "w") as f:
        for c, s in seqs.items():
            f.write(f">{c}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")

    genes, g = [], 0
    for c, L in chroms.items():
        pos = 2000
        while pos < L - 12000 and g < n_genes:
            g += 1
            gene_strand = "+" if rng.random() < 0.5 else "-"
            n_ex = rng.randint(1, 4)
            exons, p = [], pos
            for _ in range(n_ex):
                el = rng.randint(250, 700)
                exons.append((p, p + el - 1))
                p += el + rng.randint(300, 1500)
            genes.append({"id": f"GENE{g:04d}", "chrom": c, "strand": gene_strand, "exons": exons})
            pos = p + rng.randint(1500, 4000)
    with open(out / "annotation.gtf", "w") as f:
        f.write("#!genome-build SYNTH1\n")
        for ge in genes:
            s, e = ge["exons"][0][0], ge["exons"][-1][1]
            gid, tid = ge["id"], ge["id"].replace("GENE", "TX") + ".1"
            attr_g = f'gene_id "{gid}"; gene_name "Syn{gid[4:]}"; gene_type "protein_coding";'
            f.write(f"{ge['chrom']}\tsynth\tgene\t{s}\t{e}\t.\t{ge['strand']}\t.\t{attr_g}\n")
            attr_t = attr_g + f' transcript_id "{tid}";'
            f.write(f"{ge['chrom']}\tsynth\ttranscript\t{s}\t{e}\t.\t{ge['strand']}\t.\t{attr_t}\n")
            for i, (a, b) in enumerate(ge["exons"], 1):
                f.write(f"{ge['chrom']}\tsynth\texon\t{a}\t{b}\t.\t{ge['strand']}\t.\t{attr_t} exon_number \"{i}\";\n")
        # transcript sequences (sense orientation)
    for ge in genes:
        seq = "".join(seqs[ge["chrom"]][a - 1:b] for a, b in ge["exons"])
        ge["tx"] = seq if ge["strand"] == "+" else rc(seq)

    base = {ge["id"]: math.exp(rng.gauss(3.0, 1.2)) for ge in genes}
    lfc = {}
    for i, ge in enumerate(genes):
        lfc[ge["id"]] = 2.0 if i < 8 else (-2.0 if i < 16 else 0.0)
    with open(out / "truth.tsv", "w") as f:
        f.write("gene_id\ttrue_log2fc\tdirection\n")
        for ge in genes:
            d = "up" if lfc[ge["id"]] > 0 else ("down" if lfc[ge["id"]] < 0 else "none")
            f.write(f"{ge['id']}\t{lfc[ge['id']]}\t{d}\n")

    samples = [("Ctrl1", 0), ("Ctrl2", 0), ("Ctrl3", 0), ("Treat1", 1), ("Treat2", 1), ("Treat3", 1)]
    for si, (name, grp) in enumerate(samples):
        srng = random.Random(seed * 100 + si)
        # negative-binomial-like biological noise: gamma multiplier with dispersion 0.04
        w = {}
        for ge in genes:
            mu = base[ge["id"]] * (2 ** (lfc[ge["id"]] * grp)) * len(ge["tx"])
            shape = 1 / 0.04
            w[ge["id"]] = srng.gammavariate(shape, mu / shape)
        tot = sum(w.values())
        ids = [ge["id"] for ge in genes]
        cum, acc = [], 0.0
        for gid in ids:
            acc += w[gid] / tot
            cum.append(acc)
        by_id = {ge["id"]: ge for ge in genes}
        n1 = f"{name}_R1.fastq.gz" if paired else f"{name}.fastq.gz"
        n2 = f"{name}_R2.fastq.gz" if paired else None
        with gzip.open(out / "fastq" / n1, "wt", compresslevel=3) as f1, \
                (gzip.open(out / "fastq" / n2, "wt", compresslevel=3) if paired else open(os.devnull, "w")) as f2:
            for n in range(pairs_per_sample):
                r = srng.random()
                lo, hi = 0, len(cum) - 1
                while lo < hi:
                    mid = (lo + hi) // 2
                    if cum[mid] < r:
                        lo = mid + 1
                    else:
                        hi = mid
                tx = by_id[ids[lo]]["tx"]
                short = adapters and srng.random() < 0.15
                flen = srng.randint(40, 70) if short else int(srng.gauss(220, 30))
                flen = max(40, min(flen, len(tx)))
                start = srng.randint(0, len(tx) - flen)
                frag = tx[start:start + flen]
                antisense_first = strand == "reverse" or (strand == "unstranded" and srng.random() < 0.5)
                first, second = (rc(frag), frag) if antisense_first else (frag, rc(frag))
                r1 = (first + ADAPTER_R1 + "A" * read_len)[:read_len]
                r2 = (second + ADAPTER_R2 + "A" * read_len)[:read_len]
                r1, q1 = mutate(r1, srng, read_len)
                r2, q2 = mutate(r2, srng, read_len)
                rid = f"@SYN{si}:{n + 1}"
                f1.write(f"{rid} 1:N:0:ACGT\n{r1}\n+\n{q1}\n")
                f2.write(f"{rid} 2:N:0:ACGT\n{r2}\n+\n{q2}\n")
        print(f"  wrote {name} ({pairs_per_sample:,} pairs)", flush=True)
    return out


def mutate(seq, rng, L):
    s = list(seq)
    q = []
    for i in range(L):
        # quality declines slightly towards the 3' end
        qual = max(2, min(40, int(rng.gauss(38 - 10 * (i / L) ** 3, 2))))
        if rng.random() < 10 ** (-qual / 10):
            s[i] = rng.choice("ACGT")
        q.append(chr(qual + 33))
    return "".join(s), "".join(q)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "synthetic"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 150000
    strand = sys.argv[3] if len(sys.argv) > 3 else "reverse"
    paired = (sys.argv[4] if len(sys.argv) > 4 else "paired") == "paired"
    make(out, n, strand=strand, paired=paired)
    print(f"done: {out}")
