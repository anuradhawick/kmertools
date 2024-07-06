# kmertools: DNA Vectorisation Tool

![GitHub License](https://img.shields.io/github/license/anuradhawick/kmertools)
[![Cargo tests](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml)
[![Clippy check](https://github.com/anuradhawick/kmertools/actions/workflows/clippy_check.yml/badge.svg)](https://github.com/anuradhawick/kmertools/actions/workflows/clippy_check.yml)
![Conda](https://img.shields.io/conda/v/bioconda/kmertools)
![Conda](https://img.shields.io/conda/dn/bioconda/kmertools)
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

**NEW:** `kmertools` is not available on bioconda at [https://anaconda.org/bioconda/kmertools](https://anaconda.org/bioconda/kmertools).

## Features

- **Oligonucleotide Frequency Vectors:** Generate frequency vectors for oligonucleotides.
- **Minimiser Binning:** Efficiently bin sequences using minimisers to reduce data complexity.
- **Chaos Game Representation (CGR):** Compute CGR vectors for DNA sequences based on k-mers or whole sequence transformation.
- **Coverage Histograms:** Create coverage histograms to analyze the depth of sequencing reads.

## Installation

### Option 1: from bioconda (recommended)

You can install `kmertools` from Bioconda at https://anaconda.org/bioconda/kmertools. Make sure you have [conda](https://docs.conda.io/en/latest/) installed.

```bash
# create conda environment and install kmertools
conda create -n kmertools -c bioconda kmertools

# activate environment
conda activate kmertools
```

### Option 2: from sources

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
source ~/.bashrc

# to save to ~/.zshrc for Mac
echo "export PATH=\$PATH:$(pwd)/target/release/" >> ~/.zshrc
source ~/.zshrc
```

## Test the installation

After setting up, run the following command to print out the `kmertools` help message.

```bash
kmertools --help
```

## Help

Please read our comprehensive [Wiki](https://github.com/anuradhawick/kmertools/wiki).

## Authors

* Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
* Vijini Mallawaarachchi [https://vijinimallawaarachchi.com](https://vijinimallawaarachchi.com)

## Support and contributions

Please get in touch via author websites or GitHub issues. Thanks!
