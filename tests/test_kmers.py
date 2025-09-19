import pykmertools as kt
from  pykmertools import utils as ktutils


def test_kmers():
    kmer_gen = kt.KmerGenerator("ACGTCC", 3)
    kmers = list(kmer_gen)
    kmers_acgt = ["ACG", "CGT", "GTC", "TCC"]

    for (fmer, _), acgt_mer in zip(kmers, kmers_acgt):
        assert ktutils.to_acgt(fmer, len(acgt_mer)) == acgt_mer
