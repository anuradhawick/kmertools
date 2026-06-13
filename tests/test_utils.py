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


def test_wide_numeric_roundtrip():
    kmer = ("ACGT" * 25)[:99]
    forward, reverse = ktutils.to_numeric(kmer)

    assert type(forward) is int
    assert type(reverse) is int
    assert forward.bit_length() > 128
    assert ktutils.to_acgt(forward, len(kmer)) == kmer


def test_wide_numeric_length_limit():
    try:
        ktutils.to_numeric("A" * 1025)
    except ValueError as error:
        assert "between 1 and 1024" in str(error)
    else:
        raise AssertionError("expected a ValueError for a 1025-base k-mer")
