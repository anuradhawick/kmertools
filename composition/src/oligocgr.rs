use kmer::{kmer::KmerGenerator, numeric_to_kmer};
use ktio::seq::{SeqFormat, Sequence, Sequences};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufWriter, Write},
};

const GB_4: usize = 4 * (1 << 30);
type Point = (f64, f64);

// Code and test adopted from https://github.com/skatila/pycgr (as of 2:55 am Friday, 7 June 2024 Coordinated Universal Time (UTC))
// Git tag 922ebd4bec482ae6522453c14d3f9dc5d1c99995
// Under full compliance of GPL-3.0 license (https://github.com/skatila/pycgr/blob/922ebd4bec482ae6522453c14d3f9dc5d1c99995/LICENSE)
pub struct OligoCgrComputer {
    in_path: String,
    out_path: String,
    threads: usize,
    norm: bool,
    ksize: usize,
    memory: usize,
    cgr_center: Point,
    cgr_map: HashMap<u8, Point>,
    kmers: Vec<String>,
    pos_map: Vec<usize>,
    kcount: usize,
}

impl OligoCgrComputer {
    pub fn new(in_path: String, out_path: String, ksize: usize, vecsize: usize) -> Self {
        let (cgr_center, cgr_map) = OligoCgrComputer::cgr_maps(vecsize as f64);
        let (min_mer_pos_map, pos_min_mer_map, kcount) = KmerGenerator::kmer_pos_maps(ksize);
        let mut kmers = vec![String::new(); kcount];

        for (&pos, &kmer) in pos_min_mer_map.iter() {
            kmers[pos] = numeric_to_kmer(kmer, ksize);
        }

        Self {
            in_path,
            out_path,
            ksize,
            threads: rayon::current_num_threads(),
            norm: true,
            memory: GB_4,
            cgr_center,
            cgr_map,
            kmers,
            pos_map: min_mer_pos_map,
            kcount,
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn set_norm(&mut self, norm: bool) {
        self.norm = norm;
    }

    pub fn vectorise(&self) -> Result<(), String> {
        let mut reader = ktio::seq::get_reader(&self.in_path).unwrap();
        let buffer = reader
            .fill_buf()
            .map_err(|_| String::from("Invalid stream"))?;
        let format = if buffer[0] == b'>' {
            SeqFormat::Fasta
        } else {
            SeqFormat::Fastq
        };
        let records = Sequences::new(format, reader).unwrap();
        let file = File::create(&self.out_path)
            .map_err(|_| format!("Unable to write to file: {}", self.out_path))?;
        let mut out_buffer = BufWriter::new(file);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();

        pool.install(|| {
            rayon::scope(|_| {
                let mut buffer = Vec::with_capacity(1000);
                let mut total = 0_usize;

                // Define a closure to handle buffer processing
                let mut process_buffer = |buffer: &Vec<Sequence>| {
                    let result = buffer
                        .par_iter()
                        .map(|seq| {
                            let kvec = self.vectorise_one(&seq.seq).unwrap();
                            let kvec_str: Vec<String> = kvec
                                .iter()
                                .map(|val| format!("({},{},{})", val.0 .0, val.0 .1, val.1))
                                .collect();
                            format!("{}\n", kvec_str.join(" "))
                        })
                        .collect::<Vec<String>>()
                        .join("");
                    out_buffer.write_all(result.as_bytes()).unwrap();
                };

                for record in records {
                    total += record.seq.len();
                    buffer.push(record);

                    if total >= self.memory {
                        process_buffer(&buffer);
                        buffer.clear();
                        total = 0;
                    }
                }

                if !buffer.is_empty() {
                    process_buffer(&buffer);
                }
            });
        });
        Ok(())
    }

    fn vectorise_one(&self, seq: &[u8]) -> Result<Vec<(Point, f64)>, String> {
        let mut cgr = Vec::with_capacity(seq.len());
        let freqs = self.seq_to_kmer(seq);

        for (kmer, freq) in self.kmers.iter().zip(freqs.iter()) {
            let mut cgr_marker = self.cgr_center;
            for s in kmer.as_bytes() {
                if let Some(&cgr_corner) = self.cgr_map.get(s) {
                    cgr_marker = (
                        (cgr_corner.0 + cgr_marker.0) / 2.0,
                        (cgr_corner.1 + cgr_marker.1) / 2.0,
                    );
                } else {
                    return Err("Bad nucleotide, unable to proceed".to_string());
                }
            }
            cgr.push((cgr_marker, *freq));
        }

        Ok(cgr)
    }

    fn seq_to_kmer(&self, seq: &[u8]) -> Vec<f64> {
        let mut vec = vec![0_f64; self.kcount];
        let mut total = 0_f64;

        for (fmer, rmer) in KmerGenerator::new(seq, self.ksize) {
            let min_mer = u64::min(fmer, rmer);
            unsafe {
                // we already know the size of the vector and
                // min_mer is absolutely smaller than that
                let &min_mer_pos = self.pos_map.get_unchecked(min_mer as usize);
                *vec.get_unchecked_mut(min_mer_pos) += 1_f64;
                total += 1_f64;
            }
        }
        if self.norm {
            vec.iter_mut().for_each(|el| *el /= f64::max(1_f64, total));
        }
        vec
    }

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
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn oligo_cgr_vec_norm_test() {
        let cgr = OligoCgrComputer::new(PATH_FQ.to_owned(), "".to_owned(), 4, 16);
        let res = cgr
            .vectorise_one("aaaatgatgaaatagagagactttattaa".as_bytes())
            .unwrap();
        assert_eq!(res[0].0 .0, 0.5);
        assert_eq!(res[0].0 .1, 0.5);
        assert_eq!(res[0].1, 1.0 / (29 - 4 + 1) as f64);
    }

    #[test]
    fn oligo_cgr_vec_unnorm_test() {
        let mut cgr = OligoCgrComputer::new(PATH_FQ.to_owned(), "".to_owned(), 4, 16);
        cgr.set_norm(false);
        let res = cgr
            .vectorise_one("aaaatgatgaaatagagagactttattaa".as_bytes())
            .unwrap();
        assert_eq!(res[0].0 .0, 0.5);
        assert_eq!(res[0].0 .1, 0.5);
        assert_eq!(res[0].1, 1.0);
    }

    #[test]
    fn oligo_cgr_complete_unnorm_test() {
        let mut cgr = OligoCgrComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/reads.k4.cgr".to_owned(),
            4,
            16,
        );
        cgr.set_threads(4);
        cgr.set_norm(false);
        cgr.vectorise().unwrap();

        assert_eq!(
            fs::read("../test_data/expected_reads.k4.cgr").unwrap(),
            fs::read("../test_data/reads.k4.cgr").unwrap()
        )
    }
}
