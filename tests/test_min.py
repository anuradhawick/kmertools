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
