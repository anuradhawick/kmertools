import pykmertools as kt
from pykmertools import utils as ktutils


def test_kmers():
    kmer_gen = kt.KmerGenerator("ACGTCC", 3)
    kmers = list(kmer_gen)
    kmers_acgt = ["ACG", "CGT", "GTC", "TCC"]

    for (fmer, _), acgt_mer in zip(kmers, kmers_acgt):
        assert ktutils.to_acgt(fmer, len(acgt_mer)) == acgt_mer


def test_wide_kmers():
    sequence = "ACGT" * 24
    ksize = 70
    fmer, rmer = next(iter(kt.KmerGenerator(sequence, ksize)))

    assert type(fmer) is int
    assert type(rmer) is int
    assert fmer.bit_length() > 128
    assert ktutils.to_acgt(fmer, ksize) == sequence[:ksize]
    assert (fmer, rmer) == ktutils.to_numeric(sequence[:ksize])
