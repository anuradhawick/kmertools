use kmer::KmerGenerator;
use ktio::mmap::MMWriter;
use ktio::seq::{SeqFormat, Sequences};
use std::sync::Arc;
use std::sync::Mutex;

const NUMBER_SIZE: usize = 8;

pub struct CompositionComputer<'a> {
    in_path: &'a str,
    out_path: &'a str,
    ksize: usize,
    kcount: usize,
    scount: usize,
    threads: usize,
    pos_map: Vec<usize>,
    norm: bool,
}

impl<'a> CompositionComputer<'a> {
    pub fn new(
        in_path: &'a str,
        out_path: &'a str,
        ksize: usize,
        norm: bool,
        threads: usize,
    ) -> Self {
        let (pos_map, kcount) = KmerGenerator::kmer_to_vec_pos_map(ksize);
        let mut comp = CompositionComputer {
            in_path,
            out_path,
            ksize,
            kcount,
            scount: 0,
            pos_map,
            threads,
            norm,
        };
        comp.scount = comp.count_seqs();
        comp
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

    fn count_seqs(&mut self) -> usize {
        let reader = ktio::seq::get_reader(self.in_path).unwrap();
        let seqs = Sequences::new(SeqFormat::get(self.in_path).unwrap(), reader).unwrap();
        seqs.count()
    }

    fn vectorise_mmap(&self) -> Result<(), String> {
        let per_line_size = self.kcount * (NUMBER_SIZE + 1);
        let estimated_file_size = self.scount * per_line_size;
        let mut mmap = ktio::mmap::mmap_file_for_writing(self.out_path, estimated_file_size)?;
        let reader = ktio::seq::get_reader(self.in_path).unwrap();
        let format = SeqFormat::get(self.in_path).unwrap();
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
                            let kvec_str: Vec<String> = kvec
                                .iter()
                                .map(|val| format!("{:.*}", NUMBER_SIZE - 2, val))
                                .collect();
                            let kvec_str = format!("{}\n", kvec_str.join(" "));
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
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn seq_count_test() {
        let com = CompositionComputer::new(
            "../test_data/reads.fq",
            "../test_data/reads.kmers",
            4,
            true,
            8,
        );
        assert_eq!(com.scount, 2);
    }

    #[test]
    fn kmer_vec_test() {
        let com = CompositionComputer::new(
            "../test_data/reads.fq",
            "../test_data/reads.kmers",
            4,
            true,
            8,
        );
        let kvec = com.vectorise_one("AAAANGAGA");
        assert_eq!(kvec[0], 1.0);
    }

    #[test]
    fn kmer_vec_norm_test() {
        let com = CompositionComputer::new(
            "../test_data/reads.fq",
            "../test_data/reads.kmers",
            4,
            true,
            8,
        );
        let kvec = com.vectorise_one("AAAANGAGA");
        assert_eq!(kvec[0], 0.5);
    }

    #[test]
    fn vec_mmap_test() {
        let com = CompositionComputer::new(
            "../test_data/reads.fq",
            "../test_data/computed_fa.kmers",
            4,
            true,
            8,
        );
        let _ = com.vectorise_mmap();
        assert_eq!(com.scount, 2);
        assert_eq!(
            fs::read("../test_data/expected_fa.kmers").unwrap(),
            fs::read("../test_data/computed_fa.kmers").unwrap()
        )
    }
}
