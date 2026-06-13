use std::{mem::transmute, sync::Arc};

use crate::kmer::{word_from_py, word_to_py};
use kmer::{
    minimiser::MinimiserGenerator as RsMinimiserGenerator, numeric_to_kmer, KmerU1024, KmerU2048,
    KmerU256, KmerU512,
};
use pyo3::{exceptions::PyValueError, prelude::*};

enum MinimiserGeneratorInner {
    U64(Box<RsMinimiserGenerator<'static, u64>>),
    U128(Box<RsMinimiserGenerator<'static, u128>>),
    U256(Box<RsMinimiserGenerator<'static, KmerU256>>),
    U512(Box<RsMinimiserGenerator<'static, KmerU512>>),
    U1024(Box<RsMinimiserGenerator<'static, KmerU1024>>),
    U2048(Box<RsMinimiserGenerator<'static, KmerU2048>>),
}

impl MinimiserGeneratorInner {
    fn new(seq: &'static [u8], wsize: usize, msize: usize) -> PyResult<Self> {
        if wsize < msize {
            return Err(PyValueError::new_err(
                "window size must be at least minimiser size",
            ));
        }

        match msize {
            1..=32 => Ok(Self::U64(Box::new(RsMinimiserGenerator::new(
                seq, wsize, msize,
            )))),
            33..=64 => Ok(Self::U128(Box::new(RsMinimiserGenerator::new(
                seq, wsize, msize,
            )))),
            65..=128 => Ok(Self::U256(Box::new(RsMinimiserGenerator::new(
                seq, wsize, msize,
            )))),
            129..=256 => Ok(Self::U512(Box::new(RsMinimiserGenerator::new(
                seq, wsize, msize,
            )))),
            257..=512 => Ok(Self::U1024(Box::new(RsMinimiserGenerator::new(
                seq, wsize, msize,
            )))),
            513..=1024 => Ok(Self::U2048(Box::new(RsMinimiserGenerator::new(
                seq, wsize, msize,
            )))),
            _ => Err(PyValueError::new_err(
                "minimiser size must be between 1 and 1024 bases",
            )),
        }
    }

    fn next(&mut self, py: Python<'_>) -> PyResult<Option<(Py<PyAny>, usize, usize)>> {
        macro_rules! next_minimiser {
            ($generator:expr) => {
                $generator
                    .next()
                    .map(|(minimiser, start, end)| Ok((word_to_py(py, minimiser)?, start, end)))
                    .transpose()
            };
        }

        match self {
            Self::U64(generator) => next_minimiser!(generator),
            Self::U128(generator) => next_minimiser!(generator),
            Self::U256(generator) => next_minimiser!(generator),
            Self::U512(generator) => next_minimiser!(generator),
            Self::U1024(generator) => next_minimiser!(generator),
            Self::U2048(generator) => next_minimiser!(generator),
        }
    }
}

/// Computer for generating minimisers
#[pyclass]
pub struct MinimiserGenerator {
    generator: MinimiserGeneratorInner,
    _data: Arc<[u8]>,
    msize: usize,
}

#[pymethods]
impl MinimiserGenerator {
    /// Initialise the minimiser iterator
    /// Attributes:
    ///     seq (str): string from which to extract minimisers
    ///     wsize (int): size of the window
    ///     msize (int): size of the minimiser
    #[new]
    #[pyo3(signature = (seq, wsize, msize))]
    pub fn new(seq: String, wsize: usize, msize: usize) -> PyResult<Self> {
        let data: Arc<[u8]> = Arc::from(seq.into_boxed_str().into_boxed_bytes());
        let static_str: &'static [u8] = unsafe { transmute(Arc::as_ref(&data)) };
        let generator = MinimiserGeneratorInner::new(static_str, wsize, msize)?;
        Ok(Self {
            generator,
            _data: data,
            msize,
        })
    }

    /// Translate numeric minimiser to ACGT
    /// Attributes:
    ///     mmer (int): value of the minimiser
    #[pyo3(signature = (mmer))]
    pub fn to_acgt(&self, mmer: &Bound<'_, PyAny>) -> PyResult<String> {
        macro_rules! convert {
            ($word:ty) => {
                Ok(numeric_to_kmer(word_from_py::<$word>(mmer)?, self.msize))
            };
        }

        match self.msize {
            1..=32 => convert!(u64),
            33..=64 => convert!(u128),
            65..=128 => convert!(KmerU256),
            129..=256 => convert!(KmerU512),
            257..=512 => convert!(KmerU1024),
            513..=1024 => convert!(KmerU2048),
            _ => Err(PyValueError::new_err(
                "minimiser size must be between 1 and 1024 bases",
            )),
        }
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Returns:
    ///     Tuple[int, int, int]: minimiser, start pos, end pos
    pub fn __next__(
        mut slf: PyRefMut<'_, Self>,
        py: Python<'_>,
    ) -> PyResult<Option<(Py<PyAny>, usize, usize)>> {
        slf.generator.next(py)
    }
}
