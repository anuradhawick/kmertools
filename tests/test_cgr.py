import pykmertools as kt
import pathlib
from Bio import SeqIO

dir = pathlib.Path(__file__).parent


def test_cgr():
    cgr_gen = kt.CgrComputer(1)
    seqs = [
        str(seq.seq)
        for seq in list(SeqIO.parse(dir.joinpath("../test_data/reads.fq"), "fastq"))
    ]
    cgrs_generated = cgr_gen.vectorise_batch(seqs)
    cgrs_truth = [
        [eval(item) for item in line.split(" ")]
        for line in open(dir.joinpath("../test_data/expected_reads.cgr"))
        .read()
        .splitlines()
    ]
    for g, t in zip(cgrs_generated, cgrs_truth):
        assert g == t
