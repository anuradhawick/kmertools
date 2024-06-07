use kmer::kmer::KmerGenerator;
use kmer::numeric_to_kmer;
use ktio::mmap::MMWriter;
use ktio::seq::{SeqFormat, Sequence, Sequences};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::sync::Arc;
use std::sync::Mutex;

const NUMBER_SIZE: usize = 8;
const GB_4: usize = 4 * (1 << 30);

pub struct OligoComputer {
    in_path: String,
    out_path: String,
    ksize: usize,
    kcount: usize,
    threads: usize,
    pos_map: Vec<usize>,
    pos_kmer: HashMap<usize, u64>,
    norm: bool,
    delim: String,
    memory: usize,
    header: bool,
}

impl OligoComputer {
    pub fn new(in_path: String, out_path: String, ksize: usize) -> Self {
        let (min_mer_pos_map, pos_min_mer_map, kcount) = KmerGenerator::kmer_pos_maps(ksize);
        Self {
            in_path,
            out_path,
            ksize,
            kcount,
            pos_map: min_mer_pos_map,
            pos_kmer: pos_min_mer_map,
            threads: rayon::current_num_threads(),
            norm: true,
            delim: " ".to_owned(),
            memory: GB_4,
            header: false,
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

    pub fn set_max_memory(&mut self, memory: usize) {
        self.memory = memory;
    }

    pub fn set_header(&mut self, header: bool) {
        self.header = header;
    }

    fn get_header(&self) -> Vec<String> {
        let mut kmers = vec![String::new(); self.kcount];
        for (&pos, &kmer) in self.pos_kmer.iter() {
            kmers[pos] = numeric_to_kmer(kmer, self.ksize);
        }
        kmers
    }

    // this function cannot be fully tested becaue we cannot have stdin at test time
    // TODO remove stdin if needed
    #[cfg(not(tarpaulin_include))]
    pub fn vectorise(&self) -> Result<(), String> {
        if self.in_path == "-" || !self.norm {
            return self.vectorise_batch();
        }
        self.vectorise_mmap()
    }

    fn vectorise_batch(&self) -> Result<(), String> {
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

        if self.header {
            let header = self.get_header().join(&self.delim) + "\n";
            out_buffer.write_all(header.as_bytes()).unwrap();
        }

        pool.install(|| {
            rayon::scope(|_| {
                let mut buffer = Vec::with_capacity(1000);
                let mut total = 0_usize;

                // Define a closure to handle buffer processing
                let mut process_buffer = |buffer: &Vec<Sequence>| {
                    let result = buffer
                        .par_iter()
                        .map(|seq| {
                            let kvec = self.vectorise_one(&seq.seq);
                            let kvec_str: Vec<String> = kvec
                                .iter()
                                .map(|val| {
                                    if self.norm {
                                        format!("{:.*}", NUMBER_SIZE - 2, val)
                                    } else {
                                        format!("{}", val)
                                    }
                                })
                                .collect();
                            format!("{}\n", kvec_str.join(&self.delim))
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

    fn vectorise_mmap(&self) -> Result<(), String> {
        // only works for normalised (we need fixed length outputs)
        assert!(self.norm);
        let per_line_size = self.kcount * (NUMBER_SIZE + 1);
        // pre-calculate file size
        let mut estimated_file_size = {
            let format = SeqFormat::get(&self.in_path).unwrap();
            let reader = ktio::seq::get_reader(&self.in_path).unwrap();
            Sequences::seq_stats(format, reader).seq_count
        } * per_line_size;
        let mut header = String::new();
        if self.header {
            header = self.get_header().join(&self.delim) + "\n";
            estimated_file_size += header.len();
        }
        // memmap
        let mut mmap = ktio::mmap::mmap_file_for_writing(&self.out_path, estimated_file_size)?;
        // get reader
        let format = SeqFormat::get(&self.in_path).unwrap();
        let reader = ktio::seq::get_reader(&self.in_path).unwrap();
        let records = Sequences::new(format, reader).unwrap();
        let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();
        let records_arc = Arc::new(Mutex::new(records));

        pool.scope(|scope| {
            let mm_slice: MMWriter<u8> = MMWriter::new(&mut mmap[..]);
            if self.header {
                unsafe {
                    mm_slice.write_at(header.as_bytes(), 0);
                }
            }
            for _ in 0..self.threads {
                let records_arc_clone = Arc::clone(&records_arc);
                let header_len = header.len();
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
                                mm_slice.write_at(kvec_str.as_bytes(), start_pos + header_len);
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

    fn vectorise_one(&self, seq: &[u8]) -> Vec<f64> {
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
    use super::*;
    use std::fs;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn kmer_vec_norm_test() {
        let com = OligoComputer::new(PATH_FQ.to_owned(), "../test_data/reads.kmers".to_owned(), 4);
        let kvec = com.vectorise_one(b"AAAANGAGA");
        assert_eq!(kvec[0], 0.5);
    }

    #[test]
    fn kmer_vec_unnorm_test() {
        let mut com =
            OligoComputer::new(PATH_FQ.to_owned(), "../test_data/reads.kmers".to_owned(), 4);
        com.set_norm(false);
        let kvec = com.vectorise_one(b"AAAANGAGA");
        assert_eq!(kvec[0], 1.0);
        assert_eq!(kvec.iter().fold(0.0, |acc, v| acc + v), 2.0);
    }

    #[test]
    fn vec_mmap_test() {
        let com = OligoComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_fa.kmers".to_owned(),
            4,
        );
        let _ = com.vectorise_mmap();
        assert_eq!(
            fs::read("../test_data/expected_fa.kmers").unwrap(),
            fs::read("../test_data/computed_fa.kmers").unwrap()
        )
    }

    #[test]
    fn vec_mmap_threaded_test() {
        for _ in 0..4 {
            let mut com = OligoComputer::new(
                PATH_FQ.to_owned(),
                "../test_data/computed_fa_mmap.kmers".to_owned(),
                4,
            );
            com.set_threads(8);
            let _ = com.vectorise_mmap();
            assert_eq!(
                fs::read("../test_data/expected_fa.kmers").unwrap(),
                fs::read("../test_data/computed_fa_mmap.kmers").unwrap()
            );
        }
    }

    #[test]
    fn vec_batch_threaded_test() {
        for _ in 0..4 {
            let mut com = OligoComputer::new(
                PATH_FQ.to_owned(),
                "../test_data/computed_fa_batch.kmers".to_owned(),
                4,
            );
            com.set_max_memory(100);
            com.set_threads(8);
            let _ = com.vectorise_batch();
            assert_eq!(
                fs::read("../test_data/expected_fa.kmers").unwrap(),
                fs::read("../test_data/computed_fa_batch.kmers").unwrap()
            );
            com.set_max_memory(GB_4);
            com.set_threads(8);
            let _ = com.vectorise_batch();
            assert_eq!(
                fs::read("../test_data/expected_fa.kmers").unwrap(),
                fs::read("../test_data/computed_fa_batch.kmers").unwrap()
            );
        }
    }

    #[test]
    fn vec_batch_threaded_unnorm_test() {
        for _ in 0..4 {
            let mut com = OligoComputer::new(
                PATH_FQ.to_owned(),
                "../test_data/computed_fa_batch_unnorm.kmers".to_owned(),
                4,
            );
            com.set_threads(8);
            com.set_norm(false);
            let _ = com.vectorise_batch();
            assert_eq!(
                fs::read("../test_data/expected_fa_batch_unnorm.kmers").unwrap(),
                fs::read("../test_data/computed_fa_batch_unnorm.kmers").unwrap()
            );
        }
    }

    #[test]
    fn get_header_test() {
        let com = OligoComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_fa_batch_unnorm.kmers".to_owned(),
            4,
        );
        let header = com.get_header();
        assert_eq!(header[0], "AAAA");
        assert_eq!(header[135], "TTAA");
    }

    #[test]
    fn vec_batch_with_header_test() {
        let mut com = OligoComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_fa_batch_header.kmers".to_owned(),
            4,
        );
        com.set_header(true);
        let _ = com.vectorise_batch();
        assert_eq!(
            fs::read("../test_data/computed_fa_batch_header.kmers").unwrap(),
            fs::read("../test_data/expected_fa_header.kmers").unwrap()
        );
    }

    #[test]
    fn vec_mmap_with_header_test() {
        let mut com = OligoComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_fa_mmap_header.kmers".to_owned(),
            4,
        );
        com.set_header(true);
        let _ = com.vectorise_mmap();
        assert_eq!(
            fs::read("../test_data/computed_fa_mmap_header.kmers").unwrap(),
            fs::read("../test_data/expected_fa_header.kmers").unwrap()
        );
    }
}
