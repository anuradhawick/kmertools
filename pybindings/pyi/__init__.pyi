# pykmertools/__init__.pyi

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

from .cgr_computer import CgrComputer
from .kmer_generator import KmerGenerator
from .minimiser_generator import MinimiserGenerator
from .oligo_computer import OligoComputer

__all__ = ["CgrComputer", "KmerGenerator", "MinimiserGenerator", "OligoComputer"]
