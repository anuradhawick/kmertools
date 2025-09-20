"""
Pykmertools: kmertools python wrapper

Modules:
    OligoComputer      - computing oligonucleotide frequency vectors
                         from DNA sequences
    CgrComputer        - computing chaos game representations
                          for DNA sequences
    KmerGenerator      - an iterator object to generate k-mers
                         as (forward, reverse) numeric kmer tuples
    MinimiserGenerator - an iterator object to iterate minimisers
                         as (kmer, start, end) numeric minimiser tuples
"""

from typing import List, Tuple, Dict, Iterator, Protocol
from typing_extensions import runtime_checkable

Point = Tuple[float, float]

@runtime_checkable
class _UtilsModule(Protocol):
    @staticmethod
    def to_acgt(kmer: int, ksize: int) -> str:
        """
        Convert numeric kmer to string kmer

        Args:
            kmer (int): value of the k-mer.
            ksize (int): size of the k-mer.

        Returns:
            str: ACGT alphabetic representation of the kmer.
        """
        ...

    @staticmethod
    def to_numeric(kmer: str) -> Tuple[int, int]:
        """
        Convert string kmer to numeric kmer.

        Args:
            kmer (str): ACGT alphabetic representation of the kmer.

        Returns:
            Tuple[int, int]: A tuple containing forward and reverse k-mers

        Raises:
            ValueError: kmer length is different.
        """
        ...

utils: _UtilsModule

class CgrComputer:
    """
    Computing chaos game representations (CGR) for DNA sequences.
    """

    def __init__(self, vecsize: int) -> None:
        """
        Initialise the CGR counter.

        Args:
            vecsize (int): Size of the vector to initialise the CGR map.
        """
        ...

    def vectorise_one(self, seq: str) -> List[Point]:
        """
        Generate the CGR for a single sequence.

        Args:
            seq (str): The sequence as a string.

        Returns:
            List[Point]: A list of points representing the CGR of the sequence.

        Raises:
            ValueError: If the sequence contains an invalid nucleotide.
        """
        ...

    def vectorise_batch(self, seqs: List[str]) -> List[List[Point]]:
        """
        Generate the CGRs for a batch of sequences.

        Args:
            seqs (List[str]): A list of sequences.

        Returns:
            List[List[Point]]: A list of lists, each containing points representing the CGR of a sequence.

        Raises:
            ValueError: If any sequence contains an invalid nucleotide.
        """
        ...

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

    def kmer_pos_maps(self) -> Tuple[List[int], Dict[int, int], int]:
        """
        Get the k-mer position maps for the KmerGenerator.

        Returns:
            Tuple[List[int], Dict[int, int], int]: A tuple containing:
                - A list mapping of size 4^ksize each index corresponds to minimum complement
                  k-mer in integer representation.
                - A dictionary mapping minimum complement k-mer integer representation to its
                  index in above list.
                - The total number of possible minimum complement k-mers.
        """
        ...

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

class OligoComputer:
    """
    Computing oligonucleotide frequency vectors from DNA sequences.
    """

    def __init__(self, k: int) -> None:
        """
        Initialise the OligoComputer.

        Args:
            k (int): The size of the oligonucleotides to compute.
        """
        ...

    def vectorise_one(
        self, seq: str, norm: bool = True, mins: bool = True
    ) -> List[float]:
        """
        Generate the CGR for a single sequence.

        Args:
            seq (str): The sequence as a string.
            norm (bool): enable normalisation by counts.
            mins (bool): count minimum complement k-mers only.

        Returns:
            List[float]: A list of floats representing the oligonuclotide frequency vector of the sequence.

        Raises:
            ValueError: If the sequence contains an invalid nucleotide.
        """
        ...

    def vectorise_batch(
        self, seqs: List[str], norm: bool = True, mins: bool = True
    ) -> List[List[float]]:
        """
        Generate the CGRs for a batch of sequences.

        Args:
            seqs (List[str]): A list of sequences.
            norm (bool): enable normalisation by counts.
            mins (bool): count minimum complement k-mers only.

        Returns:
            List[List[Point]]: A list of lists, each representing the oligonuclotide frequency vector of the sequence.

        Raises:
            ValueError: If any sequence contains an invalid nucleotide.
        """
        ...

    def get_header(self, mins: bool = True) -> List[str]:
        """
        Generate the header for oligo nucleotide vector.

        Args:
            mins (bool): count minimum complement k-mers only.

        Returns:
            List[str]: A list of strings representing the header for the oligo nucleotide vector.
        """
        ...

__all__ = [
    "CgrComputer",
    "KmerGenerator",
    "MinimiserGenerator",
    "OligoComputer",
    "utils",
]
