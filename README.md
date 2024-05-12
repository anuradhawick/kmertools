# kmertools: DNA Vectorisation Tool

![GitHub License](https://img.shields.io/github/license/anuradhawick/kmertools)
[![Cargo tests](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml)
[![codecov](https://codecov.io/gh/anuradhawick/kmertools/graph/badge.svg?token=IDGRE54SSQ)](https://codecov.io/gh/anuradhawick/kmertools)

## Overview

`kmertools` is a k-mer based feature extraction tool designed to support metagenomics and other bioinformatics analytics. This tool leverages k-mer analysis to vectorize DNA sequences, facilitating the use of these vectors in various AI/ML applications.

## Features

- **Oligonucleotide Frequency Vectors:** Generate frequency vectors for oligonucleotides.
- **Chaos Game Representation (CGR):** Compute CGR vectors for DNA sequences.
- **Coverage Histograms:** Create coverage histograms to analyze the depth of sequencing reads.
- **Minimiser Binning:** Efficiently bin sequences using minimisers to reduce data complexity.

## Installation

You can install `kmertools` directly from the source by cloning the repository and using Rust's package manager `cargo`.

```bash
git clone https://github.com/your-repository/kmertools.git
cd kmertools
cargo build --release
```

## Usage

### General Syntax

The tool uses a command-line interface with the following general syntax:
```
kmertools <command> [options]
```

### Commands

- **Comp**: Processes the input file and outputs vectors.
- **Cov**: Generates coverage histograms based on the reads.
- **Min**: Bins reads using minimisers.


## License

TBD (GPL 3)

## Author

Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
