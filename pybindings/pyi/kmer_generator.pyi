# pykmertools/kmer_generator.pyi

from typing import Iterator, Tuple

class KmerGenerator:
    """
    An iterator object to generate k-mers as (forward, reverse) numeric kmer tuples.
    """

    def __init__(self, seq: str, ksize: int) -> None:
        """
        Initialise the KmerGenerator.

        Args:
            seq (str): The DNA sequence to generate k-mers from.
            ksize (int): The size of k-mers to generate.
        """
        ...

    def __iter__(self) -> Iterator[Tuple[int, int]]:
        """
        Return an iterator that yields (forward, reverse) numeric kmer tuples.

        Returns:
            Iterator[Tuple[int, int]]: An iterator over k-mer tuples (forward and reverse strands).
        """
        ...

    def to_acgt(self, kmer: int) -> str:
        """
        Initialise the KmerGenerator.

        Args:
            kmer (int): value of the k-mer.

        Returns:
            str: ACGT alphabetic representation of the kmer.
        """
        ...
