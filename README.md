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

Please read our comprehensive [Wiki](https://github.com/anuradhawick/kmertools/wiki).

## Authors

* Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
* Vijini Mallawaarachchi [https://vijinimallawaarachchi.com](https://vijinimallawaarachchi.com)

## Support and contributions

Please get in touch via author websites or GitHub issues. Thanks!