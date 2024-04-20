use kmer::KmerGenerator;
use ktio::mmap::MMWriter;
use ktio::seq::{SeqFormat, Sequences};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;
use std::sync::Mutex;

const NUMBER_SIZE: usize = 8;
const GB_4: usize = 2 << 30;

pub struct CompositionComputer<'a> {
    in_path: &'a str,
    out_path: &'a str,
    ksize: usize,
    kcount: usize,
    threads: usize,
    pos_map: Vec<usize>,
    norm: bool,
    delim: String,
}

impl<'a> CompositionComputer<'a> {
    pub fn new(in_path: &'a str, out_path: &'a str, ksize: usize) -> Self {
        let (pos_map, kcount) = KmerGenerator::kmer_to_vec_pos_map(ksize);
        Self {
            in_path,
            out_path,
            ksize,
            kcount,
            pos_map,
            threads: rayon::current_num_threads(),
            norm: true,
            delim: " ".to_owned(),
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

    pub fn vectorise(&self) -> Result<(), String> {
        if self.in_path == "-" {
            return self.vectorise_batch();
        }
        self.vectorise_mmap()
    }

    fn vectorise_batch(&self) -> Result<(), String> {
        let reader = ktio::seq::get_reader(self.in_path).unwrap();
        let format = SeqFormat::get("stdin.fa").unwrap();
        let records = Sequences::new(format, reader).unwrap();
        let file = File::create(self.out_path)
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

                for record in records {
                    total += record.seq.len();
                    buffer.push(record);

                    if total >= GB_4 {
                        let result = buffer
                            .par_iter()
                            .map(|seq| {
                                let kvec = self.vectorise_one(&seq.seq);
                                // optimise this with pre-sized string
                                let kvec_str: Vec<String> = kvec
                                    .iter()
                                    .map(|val| format!("{:.*}", NUMBER_SIZE - 2, val))
                                    .collect();
                                format!("{}\n", kvec_str.join(&self.delim))
                            })
                            .collect::<Vec<String>>()
                            .join("");
                        out_buffer.write_all(result.as_bytes()).unwrap();
                        buffer.clear();
                        total = 0;
                    }
                }

                if total > 0 {
                    // optimise this with pre-sized string
                    let result = buffer
                        .par_iter()
                        .map(|seq| {
                            let kvec = self.vectorise_one(&seq.seq);
                            let kvec_str: Vec<String> = kvec
                                .iter()
                                .map(|val| format!("{:.*}", NUMBER_SIZE - 2, val))
                                .collect();
                            format!("{}\n", kvec_str.join(&self.delim))
                        })
                        .collect::<Vec<String>>()
                        .join("");
                    out_buffer.write_all(result.as_bytes()).unwrap();
                    buffer.clear();
                }
            });
        });

        Ok(())
    }

    fn vectorise_mmap(&self) -> Result<(), String> {
        let per_line_size = self.kcount * (NUMBER_SIZE + 1);
        // pre-calculate file size
        let estimated_file_size = {
            let format = SeqFormat::get(self.in_path).unwrap();
            let reader = ktio::seq::get_reader(self.in_path).unwrap();
            let records = Sequences::new(format, reader).unwrap();
            records.count()
        } * per_line_size;
        // memmap
        let mut mmap = ktio::mmap::mmap_file_for_writing(self.out_path, estimated_file_size)?;
        // get reader
        let format = SeqFormat::get(self.in_path).unwrap();
        let reader = ktio::seq::get_reader(self.in_path).unwrap();
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
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn seq_count_test() {
        let path = "../test_data/reads.fq";
        let reader = ktio::seq::get_reader(path).unwrap();
        let seqs = Sequences::new(SeqFormat::get(path).unwrap(), reader).unwrap();
        assert_eq!(seqs.count(), 2);
    }

    #[test]
    fn kmer_vec_norm_test() {
        let com = CompositionComputer::new("../test_data/reads.fq", "../test_data/reads.kmers", 4);
        let kvec = com.vectorise_one("AAAANGAGA");
        assert_eq!(kvec[0], 0.5);
    }

    #[test]
    fn kmer_vec_unnorm_test() {
        let mut com =
            CompositionComputer::new("../test_data/reads.fq", "../test_data/reads.kmers", 4);
        com.set_norm(false);
        let kvec = com.vectorise_one("AAAANGAGA");
        assert_eq!(kvec[0], 1.0);
    }

    #[test]
    fn vec_mmap_test() {
        let com =
            CompositionComputer::new("../test_data/reads.fq", "../test_data/computed_fa.kmers", 4);
        let _ = com.vectorise_mmap();
        assert_eq!(
            fs::read("../test_data/expected_fa.kmers").unwrap(),
            fs::read("../test_data/computed_fa.kmers").unwrap()
        )
    }

    #[test]
    fn vec_mmap_threaded_test() {
        for _ in 0..4 {
            let mut com = CompositionComputer::new(
                "../test_data/reads.fq",
                "../test_data/computed_fa.kmers",
                4,
            );
            com.set_threads(8);
            let _ = com.vectorise_mmap();
            assert_eq!(
                fs::read("../test_data/expected_fa.kmers").unwrap(),
                fs::read("../test_data/computed_fa.kmers").unwrap()
            )
        }
    }
}
