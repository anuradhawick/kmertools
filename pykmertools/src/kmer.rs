use std::{mem::transmute, sync::Arc};

use kmer::{kmer::KmerGenerator as RsKmerGenerator, numeric_to_kmer, Kmer};
use pyo3::prelude::*;

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

    /// Translate numeric k-mer to ACGT
    /// Attributes:
    ///     kmer (int): value of the k-mer
    #[pyo3(signature = (kmer))]
    pub fn to_acgt(&self, kmer: u64) -> String {
        numeric_to_kmer(kmer, self.ksize)
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(Kmer, Kmer)> {
        slf._kg.next()
    }
}
