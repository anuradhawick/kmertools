use pyo3::{exceptions::PyValueError, prelude::*};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;

type Point = (f64, f64);

fn cgr_maps(vecsize: f64) -> (Point, HashMap<u8, Point>) {
    let cgr_a: Point = (0.0, 0.0);
    let cgr_t: Point = (vecsize, 0.0);
    let cgr_g: Point = (vecsize, vecsize);
    let cgr_c: Point = (0.0, vecsize);
    let cgr_center: Point = (vecsize / 2.0, vecsize / 2.0);

    let cgr_dict: HashMap<u8, Point> = [
        (b'A', cgr_a), // Adenine
        (b'T', cgr_t), // Thymine
        (b'G', cgr_g), // Guanine
        (b'C', cgr_c), // Cytosine
        (b'U', cgr_t), // Uracil (demethylated form of thymine)
        (b'a', cgr_a), // Adenine
        (b't', cgr_t), // Thymine
        (b'g', cgr_g), // Guanine
        (b'c', cgr_c), // Cytosine
        (b'u', cgr_t), // Uracil/Thymine
    ]
    .iter()
    .cloned()
    .collect();

    (cgr_center, cgr_dict)
}

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
    ///     norm (bool): enable normalisation by counts
    #[pyo3(signature = (seqs))]
    fn vectorise_batch(&self, seqs: Vec<String>) -> PyResult<Vec<Vec<Point>>> {
        seqs.into_par_iter()
            .map(|seq| self.vectorise_one(seq))
            .collect()
    }
}
