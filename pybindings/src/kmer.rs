use kmer::{
    kmer::KmerGenerator as RsKmerGenerator, kmer_to_numeric, numeric_to_kmer, KmerU1024, KmerU2048,
    KmerU256, KmerU512, KmerWord,
};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyAnyMethods};
use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    mem::transmute,
    str::FromStr,
    sync::Arc,
};

pub(crate) trait PyKmerWord: KmerWord + Display + FromStr
where
    <Self as FromStr>::Err: Debug,
{
}

impl<K> PyKmerWord for K
where
    K: KmerWord + Display + FromStr,
    <K as FromStr>::Err: Debug,
{
}

// Convert any KmerWord to Python's arbitrary-precision int through exact decimal text.
pub(crate) fn word_to_py<K: Display>(py: Python<'_>, value: K) -> PyResult<Py<PyAny>> {
    py.import("builtins")?
        .getattr("int")?
        .call1((value.to_string(),))
        .map(Bound::unbind)
}

// Convert a Python int to the selected fixed-width KmerWord through exact decimal text.
pub(crate) fn word_from_py<K: PyKmerWord>(value: &Bound<'_, PyAny>) -> PyResult<K>
where
    <K as FromStr>::Err: Debug,
{
    let value_string = value.str()?;
    value_string.to_str()?.parse().map_err(|error| {
        PyValueError::new_err(format!(
            "invalid non-negative k-mer integer for {} bits: {error:?}",
            K::BITS
        ))
    })
}

enum KmerGeneratorInner {
    U64(RsKmerGenerator<'static, u64>),
    U128(RsKmerGenerator<'static, u128>),
    U256(RsKmerGenerator<'static, KmerU256>),
    U512(RsKmerGenerator<'static, KmerU512>),
    U1024(RsKmerGenerator<'static, KmerU1024>),
    U2048(RsKmerGenerator<'static, KmerU2048>),
}

impl KmerGeneratorInner {
    fn new(seq: &'static [u8], ksize: usize) -> PyResult<Self> {
        match ksize {
            1..=32 => Ok(Self::U64(RsKmerGenerator::new(seq, ksize))),
            33..=64 => Ok(Self::U128(RsKmerGenerator::new(seq, ksize))),
            65..=128 => Ok(Self::U256(RsKmerGenerator::new(seq, ksize))),
            129..=256 => Ok(Self::U512(RsKmerGenerator::new(seq, ksize))),
            257..=512 => Ok(Self::U1024(RsKmerGenerator::new(seq, ksize))),
            513..=1024 => Ok(Self::U2048(RsKmerGenerator::new(seq, ksize))),
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

fn to_acgt_with<K: PyKmerWord>(kmer: &Bound<'_, PyAny>, ksize: usize) -> PyResult<String>
where
    <K as FromStr>::Err: Debug,
{
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

fn to_numeric_with<K: PyKmerWord>(py: Python<'_>, kmer: &str) -> PyResult<(Py<PyAny>, Py<PyAny>)>
where
    <K as FromStr>::Err: Debug,
{
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
