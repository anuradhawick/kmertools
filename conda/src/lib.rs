use pybindings::{
    cgr::CgrComputer,
    kmer::{register_utils_module, KmerGenerator},
    min::MinimiserGenerator,
    oligo::OligoComputer,
};
use pyo3::prelude::*;

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
    register_utils_module(m)?;
    Ok(())
}
