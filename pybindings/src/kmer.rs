use kmer::{
    kmer::KmerGenerator as RsKmerGenerator, kmer_to_numeric, numeric_to_kmer, KmerU1024, KmerU2048,
    KmerU256, KmerU512, KmerWord,
};
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyAnyMethods, PyBytes, PyBytesMethods, PyInt},
    IntoPyObjectExt,
};
use std::{collections::HashMap, mem::transmute, sync::Arc};

pub(crate) trait PyKmerWord: KmerWord {
    fn into_py_int(self, py: Python<'_>) -> PyResult<Py<PyAny>>;
    fn from_py_int(value: &Bound<'_, PyAny>) -> PyResult<Self>;
}

macro_rules! impl_native_py_kmer_word {
    ($word:ty) => {
        impl PyKmerWord for $word {
            #[inline]
            fn into_py_int(self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                self.into_py_any(py)
            }

            #[inline]
            fn from_py_int(value: &Bound<'_, PyAny>) -> PyResult<Self> {
                value.extract()
            }
        }
    };
}

impl_native_py_kmer_word!(u64);
impl_native_py_kmer_word!(u128);

macro_rules! impl_wide_py_kmer_word {
    ($word:ty, $bytes:expr) => {
        impl PyKmerWord for $word {
            fn into_py_int(self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                let bytes = self
                    .as_limbs()
                    .iter()
                    .flat_map(|limb| limb.to_le_bytes())
                    .collect::<Vec<_>>();
                let bytes = PyBytes::new(py, &bytes);
                py.get_type::<PyInt>()
                    .call_method1("from_bytes", (bytes, "little"))
                    .map(Bound::unbind)
            }

            fn from_py_int(value: &Bound<'_, PyAny>) -> PyResult<Self> {
                let bytes = value
                    .call_method1("to_bytes", ($bytes, "little"))?
                    .cast_into::<PyBytes>()?;
                Self::try_from_le_slice(bytes.as_bytes()).ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "invalid non-negative k-mer integer for {} bits",
                        Self::BITS
                    ))
                })
            }
        }
    };
}

impl_wide_py_kmer_word!(KmerU256, 32);
impl_wide_py_kmer_word!(KmerU512, 64);
impl_wide_py_kmer_word!(KmerU1024, 128);
impl_wide_py_kmer_word!(KmerU2048, 256);

// Use direct native conversion for u64/u128 and packed bytes for wider words.
pub(crate) fn word_to_py<K: PyKmerWord>(py: Python<'_>, value: K) -> PyResult<Py<PyAny>> {
    value.into_py_int(py)
}

// Extract native words directly and wider words from packed little-endian bytes.
pub(crate) fn word_from_py<K: PyKmerWord>(value: &Bound<'_, PyAny>) -> PyResult<K> {
    K::from_py_int(value)
}

enum KmerGeneratorInner {
    U64(Box<RsKmerGenerator<'static, u64>>),
    U128(Box<RsKmerGenerator<'static, u128>>),
    U256(Box<RsKmerGenerator<'static, KmerU256>>),
    U512(Box<RsKmerGenerator<'static, KmerU512>>),
    U1024(Box<RsKmerGenerator<'static, KmerU1024>>),
    U2048(Box<RsKmerGenerator<'static, KmerU2048>>),
}

impl KmerGeneratorInner {
    fn new(seq: &'static [u8], ksize: usize) -> PyResult<Self> {
        match ksize {
            1..=32 => Ok(Self::U64(Box::new(RsKmerGenerator::new(seq, ksize)))),
            33..=64 => Ok(Self::U128(Box::new(RsKmerGenerator::new(seq, ksize)))),
            65..=128 => Ok(Self::U256(Box::new(RsKmerGenerator::new(seq, ksize)))),
            129..=256 => Ok(Self::U512(Box::new(RsKmerGenerator::new(seq, ksize)))),
            257..=512 => Ok(Self::U1024(Box::new(RsKmerGenerator::new(seq, ksize)))),
            513..=1024 => Ok(Self::U2048(Box::new(RsKmerGenerator::new(seq, ksize)))),
            _ => Err(PyValueError::new_err(
                "k-mer size must be between 1 and 1024 bases",
            )),
        }
    }

    fn next(&mut self, py: Python<'_>) -> PyResult<Option<(Py<PyAny>, Py<PyAny>)>> {
        macro_rules! next_pair {
            ($generator:expr) => {
                $generator
                    .next()
                    .map(|(forward, reverse)| {
                        Ok((word_to_py(py, forward)?, word_to_py(py, reverse)?))
                    })
                    .transpose()
            };
        }

        match self {
            Self::U64(generator) => next_pair!(generator),
            Self::U128(generator) => next_pair!(generator),
            Self::U256(generator) => next_pair!(generator),
            Self::U512(generator) => next_pair!(generator),
            Self::U1024(generator) => next_pair!(generator),
            Self::U2048(generator) => next_pair!(generator),
        }
    }
}

/// Computer for generating k-mers
#[pyclass]
pub struct KmerGenerator {
    generator: KmerGeneratorInner,
    _data: Arc<[u8]>,
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
    pub fn new(seq: String, ksize: usize) -> PyResult<Self> {
        let data: Arc<[u8]> = Arc::from(seq.into_boxed_str().into_boxed_bytes());
        let static_str: &'static [u8] = unsafe { transmute(Arc::as_ref(&data)) };
        let generator = KmerGeneratorInner::new(static_str, ksize)?;
        Ok(Self {
            generator,
            _data: data,
            ksize,
        })
    }

    /// Get the k-mer position maps for the KmerGenerator
    #[pyo3(signature = ())]
    pub fn kmer_pos_maps(&self) -> PyResult<(Vec<usize>, HashMap<usize, u64>, usize)> {
        if self.ksize > 31 {
            return Err(PyValueError::new_err(
                "k-mer position maps require ksize <= 31",
            ));
        }
        Ok(RsKmerGenerator::<u64>::kmer_pos_maps(self.ksize))
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(
        mut slf: PyRefMut<'_, Self>,
        py: Python<'_>,
    ) -> PyResult<Option<(Py<PyAny>, Py<PyAny>)>> {
        slf.generator.next(py)
    }
}

fn to_acgt_with<K: PyKmerWord>(kmer: &Bound<'_, PyAny>, ksize: usize) -> PyResult<String> {
    Ok(numeric_to_kmer(word_from_py::<K>(kmer)?, ksize))
}

/// Translate numeric k-mer to ACGT
/// Attributes:
///     kmer (int): value of the k-mer
///     ksize (int): size of the k-mer
#[pyfunction]
pub fn to_acgt(kmer: &Bound<'_, PyAny>, ksize: usize) -> PyResult<String> {
    match ksize {
        1..=32 => to_acgt_with::<u64>(kmer, ksize),
        33..=64 => to_acgt_with::<u128>(kmer, ksize),
        65..=128 => to_acgt_with::<KmerU256>(kmer, ksize),
        129..=256 => to_acgt_with::<KmerU512>(kmer, ksize),
        257..=512 => to_acgt_with::<KmerU1024>(kmer, ksize),
        513..=1024 => to_acgt_with::<KmerU2048>(kmer, ksize),
        _ => Err(PyValueError::new_err(
            "k-mer size must be between 1 and 1024 bases",
        )),
    }
}

fn to_numeric_with<K: PyKmerWord>(py: Python<'_>, kmer: &str) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
    let (forward, reverse) = kmer_to_numeric::<K>(kmer);
    Ok((word_to_py(py, forward)?, word_to_py(py, reverse)?))
}

/// Translate ACGT kmer to numeric pair
/// Attributes:
///     kmer (str): k-mer string
#[pyfunction]
pub fn to_numeric(py: Python<'_>, kmer: String) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
    match kmer.len() {
        1..=32 => to_numeric_with::<u64>(py, &kmer),
        33..=64 => to_numeric_with::<u128>(py, &kmer),
        65..=128 => to_numeric_with::<KmerU256>(py, &kmer),
        129..=256 => to_numeric_with::<KmerU512>(py, &kmer),
        257..=512 => to_numeric_with::<KmerU1024>(py, &kmer),
        513..=1024 => to_numeric_with::<KmerU2048>(py, &kmer),
        _ => Err(PyValueError::new_err(
            "k-mer length must be between 1 and 1024 bases",
        )),
    }
}

pub fn register_utils_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let utils_module = PyModule::new(parent_module.py(), "utils")?;
    utils_module.add_function(wrap_pyfunction!(to_acgt, &utils_module)?)?;
    utils_module.add_function(wrap_pyfunction!(to_numeric, &utils_module)?)?;
    parent_module.add_submodule(&utils_module)?;
    Ok(())
}
