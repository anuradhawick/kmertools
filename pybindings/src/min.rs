use std::{mem::transmute, sync::Arc};

use kmer::{minimiser::MinimiserGenerator as RsMinimiserGenerator, numeric_to_kmer, Kmer};
use pyo3::prelude::*;

/// Computer for generating k-mers
#[pyclass]
pub struct MinimiserGenerator {
    _data: Arc<[u8]>,
    _mg: RsMinimiserGenerator<'static>,
    msize: usize,
}

#[pymethods]
impl MinimiserGenerator {
    /// Initialise the kmer iterator
    /// Attributes:
    ///     seq (str): string from which to extract k-mers
    ///     wsize (int): size of the window
    ///     msize (int): size of the minimiser
    #[new]
    #[pyo3(signature = (seq, wsize, msize))]
    pub fn new(seq: String, wsize: usize, msize: usize) -> Self {
        let _data: Arc<[u8]> = Arc::from(seq.into_boxed_str().into_boxed_bytes());
        let static_str: &'static [u8] = unsafe { transmute(Arc::as_ref(&_data)) };
        let _mg = RsMinimiserGenerator::new(static_str, wsize, msize);
        Self { _mg, _data, msize }
    }

    /// Translate numeric minimiser to ACGT
    /// Attributes:
    ///     mmer (int): value of the minimiser
    #[pyo3(signature = (mmer))]
    pub fn to_acgt(&self, mmer: u64) -> String {
        numeric_to_kmer(mmer, self.msize)
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Translate numeric k-mer to ACGT
    /// Returns:
    ///     Tuple[int, int, int]: minimiser, start pos, end pos
    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(Kmer, usize, usize)> {
        slf._mg.next()
    }
}
