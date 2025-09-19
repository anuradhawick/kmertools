use kmer::{kmer::KmerGenerator as RsKmerGenerator, kmer_to_numeric, numeric_to_kmer, Kmer};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::{collections::HashMap, mem::transmute, sync::Arc};

/// Computer for generating k-mers
#[pyclass]
pub struct KmerGenerator {
    _data: Arc<[u8]>,
    _kg: RsKmerGenerator<'static>,
    ksize: usize,
}

#[pymethods]
impl KmerGenerator {
    /// Initialise the kmer iterator
    /// Attributes:
    ///     seq (str): string from which to extract k-mers
    ///     ksize (int): size of the k-mers to count
    #[new]
    #[pyo3(signature = (seq, ksize))]
    pub fn new(seq: String, ksize: usize) -> Self {
        let _data: Arc<[u8]> = Arc::from(seq.into_boxed_str().into_boxed_bytes());
        let static_str: &'static [u8] = unsafe { transmute(Arc::as_ref(&_data)) };
        let _kg = RsKmerGenerator::new(static_str, ksize);
        Self { _kg, _data, ksize }
    }

    /// Get the k-mer position maps for the KmerGenerator
    #[pyo3(signature = ())]
    pub fn kmer_pos_maps(&self) -> (Vec<usize>, HashMap<usize, u64>, usize) {
        RsKmerGenerator::kmer_pos_maps(self.ksize)
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(Kmer, Kmer)> {
        slf._kg.next()
    }
}

/// Translate numeric k-mer to ACGT
/// Attributes:
///     kmer (int): value of the k-mer
///     ksize (int): size of the k-mer
#[pyfunction]
pub fn to_acgt(kmer: u64, ksize: usize) -> String {
    numeric_to_kmer(kmer, ksize)
}

/// Translate ACGT kmer to numeric pair
/// Attributes:
///     kmer (str): k-mer string
#[pyfunction]
pub fn to_numeric(kmer: String) -> PyResult<(u64, u64)> {
    if kmer.len() > 32 {
        return Err(PyValueError::new_err(format!(
            "Invalid k-mer length: {}, must be <= 32",
            kmer.len()
        )));
    }
    Ok(kmer_to_numeric(&kmer))
}

pub fn register_utils_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let utils_module = PyModule::new(parent_module.py(), "utils")?;
    let _ = utils_module.add_function(wrap_pyfunction!(to_acgt, &utils_module)?);
    let _ = utils_module.add_function(wrap_pyfunction!(to_numeric, &utils_module)?);
    parent_module.add_submodule(&utils_module)?;
    Ok(())
}
