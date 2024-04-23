use kmer::KmerGenerator;
use ktio::{
    mmap::MMWriter,
    seq::{get_reader, SeqFormat, Sequences},
};
use scc::HashMap as SccMap;
use std::{
    cmp::min,
    collections::HashMap,
    sync::{Arc, Mutex},
};

const NUMBER_SIZE: usize = 8;
const MB_100: usize = 100 * (1 << 20);

pub struct CovComputer {
    in_path: String,
    in_path_kmer: String,
    out_path: String,
    ksize: usize,
    kmer_counts: SccMap<kmer::Kmer, f64>,
    threads: usize,
    norm: bool,
    delim: String,
    bin_size: usize,
    bin_count: usize,
}

impl CovComputer {
    pub fn new(
        in_path: String,
        out_path: String,
        ksize: usize,
        bin_size: usize,
        bin_count: usize,
    ) -> Self {
        Self {
            in_path: in_path.clone(),
            in_path_kmer: in_path,
            out_path,
            ksize,
            threads: rayon::current_num_threads(),
            norm: true,
            delim: " ".to_owned(),
            kmer_counts: SccMap::new(),
            bin_size,
            bin_count,
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn set_norm(&mut self, norm: bool) {
        self.norm = norm;
    }

    pub fn set_delim(&mut self, delim: String) {
        self.delim = delim;
    }

    pub fn set_kmer_path(&mut self, path: String) {
        self.in_path_kmer = path;
    }

    pub fn vectorise(&self) -> Result<(), String> {
        self.build_table()?;
        self.vectorise_mmap()?;
        Ok(())
    }

    fn build_table(&self) -> Result<(), String> {
        // TODO: for k-mer sizes smaller than 16; we could simply use an array of size 4^k
        let format = SeqFormat::get(&self.in_path_kmer).ok_or("Unsupported input format")?;
        let reader = get_reader(&self.in_path_kmer)?;
        let records = Sequences::new(format, reader)?;
        let records_arc = Arc::new(Mutex::new(records));

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();

        pool.scope(|scope| {
            for _ in 0..self.threads {
                let records_arc_clone = Arc::clone(&records_arc);
                scope.spawn(move |_| {
                    let mut buffer = Vec::with_capacity(1000);
                    let mut total = 0_usize;
                    let mut buffer_map: HashMap<kmer::Kmer, f64> = HashMap::new();
                    loop {
                        {
                            let mut records = records_arc_clone.lock().unwrap();

                            for record in records.by_ref() {
                                total += record.seq.len();
                                buffer.push(record.seq);

                                if total >= MB_100 {
                                    break;
                                }
                            }
                        };

                        if total > 0 {
                            for seq in &buffer {
                                let gen = KmerGenerator::new(seq, self.ksize);
                                for (fmer, rmer) in gen {
                                    buffer_map
                                        .entry(kmer::Kmer::min(fmer, rmer))
                                        .and_modify(|val| *val += 1_f64)
                                        .or_insert(1_f64);
                                }
                            }
                            buffer.clear();
                            total = 0;
                            for (kmer, count) in buffer_map.iter() {
                                self.kmer_counts
                                    .entry(*kmer)
                                    .and_modify(|val| *val += *count)
                                    .or_insert(*count);
                            }
                            buffer_map.shrink_to(0);
                        } else {
                            break;
                        }
                    }
                });
            }
        });

        Ok(())
    }

    fn vectorise_mmap(&self) -> Result<(), String> {
        let per_line_size = self.bin_count * (NUMBER_SIZE + 1);
        // pre-calculate file size
        let estimated_file_size = {
            let format = SeqFormat::get(&self.in_path).unwrap();
            let reader = ktio::seq::get_reader(&self.in_path).unwrap();
            Sequences::seq_count(format, reader)
        } * per_line_size;
        // memmap
        let mut mmap = ktio::mmap::mmap_file_for_writing(&self.out_path, estimated_file_size)?;
        // get reader
        let format = SeqFormat::get(&self.in_path).unwrap();
        let reader = ktio::seq::get_reader(&self.in_path).unwrap();
        let records = Sequences::new(format, reader).unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();
        let records_arc = Arc::new(Mutex::new(records));

        pool.scope(|scope| {
            let mm_slice: MMWriter<u8> = MMWriter::new(&mut mmap[..]);
            for _ in 0..self.threads {
                let records_arc_clone = Arc::clone(&records_arc);
                scope.spawn(move |_| {
                    loop {
                        let record = { records_arc_clone.lock().unwrap().next() };
                        if let Some(record) = record {
                            let kvec = self.vectorise_one(&record.seq);
                            // optimise this with pre-sized string
                            let kvec_str: Vec<String> = kvec
                                .iter()
                                .map(|val| format!("{:.*}", NUMBER_SIZE - 2, val))
                                .collect();
                            let kvec_str = format!("{}\n", kvec_str.join(&self.delim));
                            let start_pos = kvec_str.len() * record.n;
                            unsafe {
                                mm_slice.write_at(kvec_str.as_bytes(), start_pos);
                            }
                        } else {
                            // end of iteration
                            break;
                        }
                    }
                });
            }
        });

        Ok(())
    }

    fn vectorise_one(&self, seq: &str) -> Vec<f64> {
        let mut vec = vec![0_f64; self.bin_count];
        let mut total = 0_f64;

        for (fmer, rmer) in KmerGenerator::new(seq, self.ksize) {
            let min_mer = u64::min(fmer, rmer);
            let count = self
                .kmer_counts
                .get(&min_mer)
                .map_or(0_f64, |val| *val.get());
            total += count;
            let pos = (count / self.bin_size as f64).floor() as usize;
            let pos = min(self.bin_count - 1, pos);
            unsafe {
                *vec.get_unchecked_mut(pos) += count;
            }
        }
        if self.norm {
            vec.iter_mut().for_each(|el| *el /= f64::max(1_f64, total));
        }
        vec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn kmer_counts_test() {
        let cov = CovComputer::new(PATH_FQ.to_owned(), "".to_owned(), 4, 10, 100);
        cov.build_table().unwrap();
        // AAAA/TTTT (0) = 3
        assert_eq!(cov.kmer_counts.get(&0_u64).unwrap().get(), &3_f64);
        // TGAA/AACA (224) = 3
        assert_eq!(cov.kmer_counts.get(&224_u64).unwrap().get(), &3_f64);
    }

    #[test]
    fn kmer_vec_test() {
        let cov = CovComputer::new(PATH_FQ.to_owned(), "".to_owned(), 4, 2, 5);
        cov.build_table().unwrap();
        let vec = cov.vectorise_one("ACGTNTTTT");
        assert_eq!(vec, vec![0.25, 0.75, 0.0, 0.0, 0.0,]);
    }

    #[test]
    fn kmer_vec_full_test() {
        let cov = CovComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_fq_cov.kmers".to_owned(),
            4,
            2,
            5,
        );
        cov.vectorise().unwrap();
        assert_eq!(
            fs::read("../test_data/expected_fq_cov.kmers").unwrap(),
            fs::read("../test_data/computed_fq_cov.kmers").unwrap()
        )
    }
}
