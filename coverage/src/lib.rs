use std::sync::{Arc, Mutex};

use dashmap::DashMap;
use kmer::KmerGenerator;
use ktio::seq::{get_reader, SeqFormat, Sequences};

const MB_100: usize = 100 * (2 << 20);

pub struct CovComputer<'a> {
    in_path: &'a str,
    in_path_kmer: &'a str,
    out_path: &'a str,
    ksize: usize,
    kmer_counts: DashMap<kmer::Kmer, u32>,
    threads: usize,
    norm: bool,
    delim: String,
}

impl<'a> CovComputer<'a> {
    pub fn new(in_path: &'a str, out_path: &'a str, ksize: usize) -> Self {
        Self {
            in_path,
            in_path_kmer: in_path,
            out_path,
            ksize,
            threads: rayon::current_num_threads(),
            norm: true,
            delim: " ".to_owned(),
            kmer_counts: DashMap::new(),
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

    pub fn set_kmer_path(&mut self, path: &'a str) {
        self.in_path_kmer = path;
    }

    pub fn vectorise(&self) -> Result<(), String> {
        if self.in_path == "-" {
            return self.vectorise_batch();
        }
        self.vectorise_mmap()
    }

    fn build_table(&self) -> Result<(), String> {
        let format = SeqFormat::get(self.in_path_kmer).ok_or("Unsupported input format")?;
        let reader = get_reader(self.in_path_kmer)?;
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
                                    self.kmer_counts
                                        .entry(kmer::Kmer::min(fmer, rmer))
                                        .and_modify(|val| *val += 1)
                                        .or_insert(1);
                                }
                            }
                            buffer.clear();
                            total = 0;
                        } else {
                            break;
                        }
                    }
                });
            }
        });

        // println!("{:#?}", self.kmer_counts);

        Ok(())
    }

    fn vectorise_batch(&self) -> Result<(), String> {
        todo!()
    }

    fn vectorise_mmap(&self) -> Result<(), String> {
        todo!()
    }

    fn vectorise_one(&self, seq: &str) -> Vec<f64> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn kmer_counts_test() {
        let cov = CovComputer::new(PATH_FQ, "", 4);
        cov.build_table().unwrap();
        // AAAA/TTTT (0) = 3
        assert_eq!(*cov.kmer_counts.get(&0_u64).unwrap(), 3_u32);
        // TGAA/AACA (224) = 3
        assert_eq!(*cov.kmer_counts.get(&224_u64).unwrap(), 3_u32);
    }
}
