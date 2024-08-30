# pykmertools/minimiser_generator.pyi

from typing import Iterator, Tuple

class MinimiserGenerator:
    """
    An iterator object to iterate minimisers as (kmer, start, end) numeric minimiser tuples.
    """

    def __init__(self, seq: str, wsize: int, msize: int) -> None:
        """
        Initialise the MinimiserGenerator.

        Args:
            seq (str): The DNA sequence to generate minimisers from.
            wsize (int): size of the window.
            msize (int): size of the minimiser.
        """
        ...

    def __iter__(self) -> Iterator[Tuple[int, int, int]]:
        """
        Return an iterator that yields (kmer, start, end) numeric minimiser tuples.

        Returns:
            Iterator[Tuple[int, int, int]]: An iterator over minimiser tuples as the tuple (minimiser, start, end).
        """
        ...

    def to_acgt(self, mmer: int) -> str:
        """
        Initialise the KmerGenerator.

        Args:
            mmer (int): value of the minimiser.

        Returns:
            str: ACGT alphabetic representation of the kmer.
        """
        ...
