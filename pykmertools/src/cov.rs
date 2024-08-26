use pyo3::prelude::*;

#[pyclass]
pub struct CovComputer {}

#[pymethods]
impl CovComputer {
    #[new]
    #[pyo3(signature = ())]
    fn new() -> Self {
        Self {}
    }
}
