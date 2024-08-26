mod cgr;
mod cov;
mod oligo;
use cgr::CgrComputer;
use oligo::OligoComputer;
use pyo3::prelude::*;

/// Pykmertools: kmertools python wrapper
/// Modules:
///     OligoComputer - computing oligonucleotide frequency vectors
///                     from DNA sequences
///     CgrComputer   - computing chaos game representations
///                     for DNA sequences
#[pymodule]
fn pykmertools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<OligoComputer>()?;
    m.add_class::<CgrComputer>()?;
    Ok(())
}
