mod cgr;
mod oligo;
use cgr::CgrComputer;
use oligo::OligoComputer;
use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
fn pykmertools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<OligoComputer>()?;
    m.add_class::<CgrComputer>()?;
    Ok(())
}
