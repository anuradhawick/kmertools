from typing import List

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

    def vectorise_one(self, seq: str, norm: bool = True) -> List[float]:
        """
        Generate the CGR for a single sequence.

        Args:
            seq (str): The sequence as a string.
            norm (bool): enable normalisation by counts.

        Returns:
            List[float]: A list of floats representing the oligonuclotide frequency vector of the sequence.

        Raises:
            ValueError: If the sequence contains an invalid nucleotide.
        """
        ...

    def vectorise_batch(self, seqs: List[str], norm: bool = True) -> List[List[float]]:
        """
        Generate the CGRs for a batch of sequences.

        Args:
            seqs (List[str]): A list of sequences.
            norm (bool): enable normalisation by counts.

        Returns:
            List[List[Point]]: A list of lists, each representing the oligonuclotide frequency vector of the sequence.

        Raises:
            ValueError: If any sequence contains an invalid nucleotide.
        """
        ...
