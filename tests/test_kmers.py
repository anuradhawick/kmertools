import pykmertools as kt


def test_kmers():
    kmer_gen = kt.KmerGenerator("ACGTCC", 3)
    kmers = list(kmer_gen)
    kmers_acgt = ["ACG", "CGT", "GTC", "TCC"]

    for (fmer, _), acgt_mer in zip(kmers, kmers_acgt):
        assert kmer_gen.to_acgt(fmer) == acgt_mer
