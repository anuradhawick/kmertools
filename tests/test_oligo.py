import pykmertools as kt
import pathlib
from Bio import SeqIO

dir = pathlib.Path(__file__).parent


def test_oligo():
    oligo_gen = kt.OligoComputer(4)
    seqs = [
        str(seq.seq)
        for seq in list(SeqIO.parse(dir.joinpath("../test_data/reads.fq"), "fastq"))
    ]
    oligos_generated = [
        list(map(lambda x: round(x, 6), line))
        for line in oligo_gen.vectorise_batch(seqs)
    ]
    oligos_truth = [
        list(map(float, line.strip().split()))
        for line in open(dir.joinpath("../test_data/expected_fa.kmers"))
        .read()
        .splitlines()
    ]
    for g, t in zip(oligos_generated, oligos_truth):
        assert g == t
