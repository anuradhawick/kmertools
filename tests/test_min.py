import pykmertools as kt
import pathlib

dir = pathlib.Path(__file__).parent


def test_min():
    min_gen = kt.MinimiserGenerator(
        "ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA",
        31,
        7,
    )
    mins = [
        "ACGATAT",
        "ACGCCTA",
        "AGAGCTA",
        "AAATCCC",
        "AATCCCG",
        "AATCGAT",
        "AAAGCGC",
    ]

    for (kmer, _, _), min in zip(min_gen, mins):
        assert min_gen.to_acgt(kmer) == min


def test_wide_minimisers():
    sequence = "ACGT" * 24
    min_gen = kt.MinimiserGenerator(sequence, 80, 70)
    minimiser, start, end = next(iter(min_gen))

    assert type(minimiser) is int
    assert minimiser.bit_length() > 128
    assert (start, end) == (0, len(sequence))
    assert min_gen.to_acgt(minimiser) == kt.utils.to_acgt(minimiser, 70)
