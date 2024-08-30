from typing import List, Tuple, Dict

Point = Tuple[float, float]

class CgrComputer:
    """
    Computing chaos game representations (CGR) for DNA sequences.
    """

    cgr_center: Point
    cgr_map: Dict[int, Point]

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
