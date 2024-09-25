use composition::cgr::cgr_maps;
use pyo3::{exceptions::PyValueError, prelude::*};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;

type Point = (f64, f64);

/// Computer for generating chaos game representation (cgr)
#[pyclass]
pub struct CgrComputer {
    cgr_center: Point,
    cgr_map: HashMap<u8, Point>,
}

#[pymethods]
impl CgrComputer {
    /// Initialise the cgr counter
    /// Attributes:
    ///     ksize (int): size of the k-mers to count
    #[new]
    #[pyo3(signature = (vecsize))]
    fn new(vecsize: usize) -> Self {
        let (cgr_center, cgr_map) = cgr_maps(vecsize as f64);

        Self {
            cgr_center,
            cgr_map,
        }
    }

    /// Generate the cgr
    /// Attributes:
    ///     seq (str): sequence as a string
    #[pyo3(signature = (seq))]
    fn vectorise_one(&self, seq: String) -> PyResult<Vec<Point>> {
        let mut cgr = Vec::with_capacity(seq.len());
        let mut cgr_marker = self.cgr_center;

        for s in seq.as_bytes().iter() {
            if let Some(&cgr_corner) = self.cgr_map.get(s) {
                cgr_marker = (
                    (cgr_corner.0 + cgr_marker.0) / 2.0,
                    (cgr_corner.1 + cgr_marker.1) / 2.0,
                );
                cgr.push(cgr_marker);
            } else {
                return Err(PyValueError::new_err("Bad nucleotide, unable to proceed"));
            }
        }

        Ok(cgr)
    }

    /// Generate the cgrs
    /// Attributes:
    ///     seq (list[str]): list of sequences
    #[pyo3(signature = (seqs))]
    fn vectorise_batch(&self, seqs: Vec<String>) -> PyResult<Vec<Vec<Point>>> {
        seqs.into_par_iter()
            .map(|seq| self.vectorise_one(seq))
            .collect()
    }
}
