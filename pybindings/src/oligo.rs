use kmer::{kmer::KmerGenerator, numeric_to_kmer};
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;

/// Computer for generating oligonucleotide frequency vectors
#[pyclass]
pub struct OligoComputer {
    ksize: usize,
    kcount: usize,
    pos_map: Vec<usize>,
    pos_kmer: HashMap<usize, u64>,
}

#[pymethods]
impl OligoComputer {
    /// Initialise the kmer counter
    /// Attributes:
    ///     ksize (int): size of the k-mers to count
    #[new]
    #[pyo3(signature = (ksize))]
    fn new(ksize: usize) -> Self {
        let (min_mer_pos_map, pos_min_mer_map, kcount) = KmerGenerator::kmer_pos_maps(ksize);

        Self {
            ksize,
            kcount,
            pos_map: min_mer_pos_map,
            pos_kmer: pos_min_mer_map,
        }
    }

    /// Generate the oligo nucletide vector
    /// Attributes:
    ///     seq (str): sequence as a string
    ///     norm (bool): enable normalisation by counts
    ///     mins (bool): count minimum complement k-mers only
    #[pyo3(signature = (seq, norm=true, mins=true))]
    fn vectorise_one(&self, seq: String, norm: bool, mins: bool) -> Vec<f64> {
        let kcount = if mins {
            self.kcount
        } else {
            4_u64.pow(self.ksize as u32) as usize
        };
        let mut vec = vec![0_f64; kcount];
        let mut total = 0_f64;

        for (fmer, rmer) in KmerGenerator::new(seq.as_bytes(), self.ksize) {
            if mins {
                let min_mer = u64::min(fmer, rmer);
                unsafe {
                    // we already know the size of the vector and
                    // min_mer is absolutely smaller than that
                    let &min_mer_pos = self.pos_map.get_unchecked(min_mer as usize);
                    *vec.get_unchecked_mut(min_mer_pos) += 1_f64;
                    total += 1_f64;
                }
            } else {
                unsafe {
                    *vec.get_unchecked_mut(fmer as usize) += 1_f64;
                    total += 2_f64;
                }
            }
        }
        if norm {
            vec.iter_mut().for_each(|el| *el /= f64::max(1_f64, total));
        }
        vec
    }

    /// Generate the oligo nucletide vector
    /// Attributes:
    ///     seq (list[str]): list of sequences
    ///     norm (bool): enable normalisation by counts
    ///     mins (bool): count minimum complement k-mers only
    #[pyo3(signature = (seqs, norm=true, mins=true))]
    fn vectorise_batch(&self, seqs: Vec<String>, norm: bool, mins: bool) -> Vec<Vec<f64>> {
        seqs.into_par_iter()
            .map(|seq| self.vectorise_one(seq, norm, mins))
            .collect()
    }

    /// Generate the header for oligo nucletide vector
    #[pyo3(signature = (mins=true))]
    fn get_header(&self, mins: bool) -> Vec<String> {
        if mins {
            let mut kmers = vec![String::new(); self.kcount];
            for (&pos, &kmer) in self.pos_kmer.iter() {
                kmers[pos] = numeric_to_kmer(kmer, self.ksize);
            }
            kmers
        } else {
            let mut kmers = vec![String::new(); 4_u64.pow(self.ksize as u32) as usize];
            for kmer in 0..(4_u64.pow(self.ksize as u32)) {
                kmers[kmer as usize] = numeric_to_kmer(kmer, self.ksize);
            }
            kmers
        }
    }
}
