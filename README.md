# kmertools: DNA Vectorisation Tool

![GitHub License](https://img.shields.io/github/license/anuradhawick/kmertools)
[![Cargo tests](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml)
[![codecov](https://codecov.io/gh/anuradhawick/kmertools/graph/badge.svg?token=IDGRE54SSQ)](https://codecov.io/gh/anuradhawick/kmertools)

<div align="center">
<pre>
$$\   $$\                                   $$$$$$$$\                     $$\           
$$ | $$  |                                  \__$$  __|                    $$ |          
$$ |$$  / $$$$$$\$$$$\   $$$$$$\   $$$$$$\     $$ |    $$$$$$\   $$$$$$\  $$ | $$$$$$$\ 
$$$$$  /  $$  _$$  _$$\ $$  __$$\ $$  __$$\    $$ |   $$  __$$\ $$  __$$\ $$ |$$  _____|
$$  $$<   $$ / $$ / $$ |$$$$$$$$ |$$ |  \__|   $$ |   $$ /  $$ |$$ /  $$ |$$ |\$$$$$$\  
$$ |\$$\  $$ | $$ | $$ |$$   ____|$$ |         $$ |   $$ |  $$ |$$ |  $$ |$$ | \____$$\ 
$$ | \$$\ $$ | $$ | $$ |\$$$$$$$\ $$ |         $$ |   \$$$$$$  |\$$$$$$  |$$ |$$$$$$$  |
\__|  \__|\__| \__| \__| \_______|\__|         \__|    \______/  \______/ \__|\_______/ 
</pre>
</div>
                                                         
## Overview

`kmertools` is a k-mer based feature extraction tool designed to support metagenomics and other bioinformatics analytics. This tool leverages k-mer analysis to vectorize DNA sequences, facilitating the use of these vectors in various AI/ML applications.

## Features

- **Oligonucleotide Frequency Vectors:** Generate frequency vectors for oligonucleotides.
- **Minimiser Binning:** Efficiently bin sequences using minimisers to reduce data complexity.
- PLANNED **Chaos Game Representation (CGR):** Compute CGR vectors for DNA sequences.
- PLANNED **Coverage Histograms:** Create coverage histograms to analyze the depth of sequencing reads.

## Installation

You can install `kmertools` directly from the source by cloning the repository and using Rust's package manager `cargo`.

```bash
git clone https://github.com/your-repository/kmertools.git
cd kmertools
cargo build --release
```

Now add the binary to path (you may modify `~/.bashrc` or `~/.zshrc`)

```sh
# to add to current terminal
export PATH=$PATH:$(pwd)/target/release/
# to save to ~/.bashrc
echo "export PATH=\$PATH:$(pwd)/target/release/" >> ~/.bashrc
# to save to ~/.zshrc
echo "export PATH=\$PATH:$(pwd)/target/release/" >> ~/.zshrc
```

## Help

The tool uses a command-line interface with the following command.

```sh
kmertools --help
```

OR

```sh
kmertools <COMMAND> --help
```

You will be greeted with a help window as follows.

```sh
kmertools: DNA vectorisation

Usage: kmertools <COMMAND>

Commands:
  comp  Generate sequence composition based features
  cov   Generates coverage histogram based on the reads
  min   Bin reads using minimisers
  help  Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help (see more with '--help')
  -V, --version  Print version
```

## Commands

### Composition computations

```sh
Generate sequence composition based features

Usage: kmertools comp <COMMAND>

Commands:
  oligo  Generate oligonucleotide frequency vectors
  help   Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

### Coverage computations (Under development)

```sh
Generates coverage histogram based on the reads

Usage: kmertools cov [OPTIONS] --input <INPUT> --output <OUTPUT>

Options:
  -i, --input <INPUT>
          Input file path

  -a, --alt-input <ALT_INPUT>
          Input file path, for k-mer counting

  -o, --output <OUTPUT>
          Output vectors path

  -k, --k-size <K_SIZE>
          K size for the coverage histogram
          
          [default: 15]

  -p, --preset <PRESET>
          Output type to write
          
          [default: spc]

          Possible values:
          - csv: Comma separated format
          - tsv: Tab separated format
          - spc: Space separated format

  -s, --bin-size <BIN_SIZE>
          Bin size for the coverage histogram
          
          [default: 16]

  -c, --bin-count <BIN_COUNT>
          Number of bins for the coverage histogram
          
          [default: 16]

      --counts
          Disable normalisation and output raw counts

  -t, --threads <THREADS>
          Thread count for computations 0=auto
          
          [default: 0]

  -h, --help
          Print help (see a summary with '-h')
```

### Minimiser computations

```sh
Bin reads using minimisers

Usage: kmertools min [OPTIONS] --input <INPUT> --output <OUTPUT>

Options:
  -i, --input <INPUT>
          Input file path

  -o, --output <OUTPUT>
          Output vectors path

  -m, --m-size <M_SIZE>
          Minimiser size
          
          [default: 10]

  -w, --w-size <W_SIZE>
          Window size
          
          0 - emits one minimiser per sequence (useful for sequencing reads)
          w_size must be longer than m_size
          
          [default: 0]

  -p, --preset <PRESET>
          Output type to write
          
          [default: s2m]

          Possible values:
          - s2m: Conver sequences into minimiser representation
          - m2s: Group sequences by minimiser

  -t, --threads <THREADS>
          Thread count for computations 0=auto
          
          [default: 0]

  -h, --help
          Print help (see a summary with '-h')
```

### K-mer counting

```sh
Count k-mers

Usage: kmertools ctr [OPTIONS] --input <INPUT> --output <OUTPUT> --k-size <K_SIZE>

Options:
  -i, --input <INPUT>      Input file path
  -o, --output <OUTPUT>    Output vectors path
  -k, --k-size <K_SIZE>    k size for counting
  -m, --memory <MEMORY>    Max memory in GB [default: 6]
  -t, --threads <THREADS>  Thread count for computations 0=auto [default: 0]
  -h, --help               Print help
```

## Authors

* Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
* Vijini Mallawaarachchi [https://vijinimallawaarachchi.com](https://vijinimallawaarachchi.com)

## Support and contributions

Please get in touch via author websites or GitHub issues. Thanks!