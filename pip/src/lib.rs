use clap::Parser;
use kmertools::args::{cli, Cli};
use pybindings::{
    cgr::CgrComputer,
    kmer::{register_utils_module, KmerGenerator},
    min::MinimiserGenerator,
    oligo::OligoComputer,
};
use pyo3::prelude::*;

#[pyfunction]
// TODO: remove after https://github.com/PyO3/maturin/issues/368 is resolved
fn run_cli(_py: Python) -> PyResult<()> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    let parsed_args = Cli::parse_from(&args);
    cli(parsed_args);
    Ok(())
}

/// Pykmertools: kmertools python wrapper
/// Modules:
///     OligoComputer      - computing oligonucleotide frequency vectors
///                          from DNA sequences
///     CgrComputer        - computing chaos game representations
///                           for DNA sequences
///     KmerGenerator      - an iterator object to generate k-mers
///                          as (forward, reverse) numeric kmer tuples
///     MinimiserGenerator - an iterator object to iterate minimisers
///                          as (kmer, start, end) numeric minimiser tuples
///     utils              - utility functions
#[pymodule]
fn pykmertools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<OligoComputer>()?;
    m.add_class::<CgrComputer>()?;
    m.add_class::<KmerGenerator>()?;
    m.add_class::<MinimiserGenerator>()?;
    m.add_function(wrap_pyfunction!(run_cli, m)?)?;
    register_utils_module(m)?;
    Ok(())
}
