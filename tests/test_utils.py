from pykmertools import utils as ktutils


def test_to_acgt():
    kmer_1 = ktutils.to_acgt(111, 5)
    kmer_2 = ktutils.to_acgt(27, 5)

    assert kmer_1 == "ACGTT"
    assert kmer_2 == "AACGT"


def test_to_numeric():
    kmer_1, kmer_2 = ktutils.to_numeric("ACGTT")

    assert kmer_1 == 111
    assert kmer_2 == 27
