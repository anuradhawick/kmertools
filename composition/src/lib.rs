use kmer::KmerGenerator;
use memmap2::{MmapMut, MmapOptions};
use std::fs;
use std::fs::OpenOptions;
use std::io::BufReader;
use std::io::{Error, Write};
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};
use threadpool::ThreadPool;

fn vectorise_one(
    kcount: usize,
    ksize: usize,
    pos_map: &[usize],
    seq: &str,
    norm: bool,
) -> Vec<f64> {
    let mut vec = vec![0_f64; kcount];
    let mut total = 0_f64;

    for (fmer, rmer) in KmerGenerator::new(seq, ksize) {
        let min_mer = u64::min(fmer, rmer);
        unsafe {
            // we already know the size of the vector and
            // min_mer is absolutely smaller than that
            let &min_mer_pos = pos_map.get_unchecked(min_mer as usize);
            *vec.get_unchecked_mut(min_mer_pos) += 1_f64;
            total += 1_f64;
        }
    }
    if norm {
        vec.iter_mut().for_each(|el| *el /= f64::max(1_f64, total));
    }
    vec
}

pub struct CompositionComputer<'a> {
    in_path: &'a str,
    out_path: &'a str,
    ksize: usize,
    kcount: usize,
    scount: usize,
    pos_map: Vec<usize>,
}

impl<'a> CompositionComputer<'a> {
    pub fn new(in_path: &'a str, out_path: &'a str, ksize: usize) -> Self {
        let (pos_map, kcount) = KmerGenerator::kmer_to_vec_pos_map(ksize);
        let mut comp = CompositionComputer {
            in_path,
            out_path,
            ksize,
            kcount,
            scount: 0,
            pos_map,
        };
        comp.scount = comp.count_seqs();
        comp
    }

    fn count_seqs(&self) -> usize {
        let file = fs::File::open(self.in_path).unwrap();
        let reader = BufReader::new(file);
        let seqs = io::Sequences::new("fq", reader).unwrap();
        seqs.count()
    }

    fn create_writable_mmap(&self) -> Result<MmapMut, Error> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .truncate(true)
            .create(true)
            .open(self.out_path)?;

        let per_line_size = self.kcount * (8 + 1);
        let estimated_file_size = self.scount * per_line_size;
        file.set_len(estimated_file_size as u64)?;
        let mmap = unsafe { MmapOptions::new().len(estimated_file_size).map_mut(&file)? };

        Ok(mmap)
    }

    fn vectorise_mmap(&self) -> Result<(), String> {
        let mut mmap = self.create_writable_mmap().unwrap();
        let file = fs::File::open(self.in_path).unwrap();
        let reader = BufReader::new(file);
        let records = Arc::new(Mutex::new(io::Sequences::new("fq", reader).unwrap()));
        let pool: ThreadPool = ThreadPool::new(8);
        let pos_map = self.pos_map.clone();
        let kcount = Arc::new(self.kcount);
        let ksize = Arc::new(self.ksize);

        for _ in 0..8 {
            let records_clone = Arc::clone(&records);
            let kcount_clone = Arc::clone(&kcount);
            let ksize_clone = Arc::clone(&ksize);
            let pos_map_clone = pos_map.clone();
            pool.execute(move || {
                while let Some(record) = records_clone.lock().unwrap().next() {
                    let kvec = vectorise_one(
                        *kcount_clone,
                        *ksize_clone,
                        &pos_map_clone,
                        &record.seq,
                        true,
                    );
                    println!("processed {:?} {:#?}", record.id, &kvec.as_slice()[0..5]);
                }
            });
            todo!("Implement file writing")
        }
        pool.join();

        Ok(())
    }

    fn vectorise_batch(&self) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq_count_test() {
        let com = CompositionComputer::new("../test_data/reads.fq", "../test_data/reads.kmers", 4);
        assert_eq!(com.count_seqs(), 2);
    }

    #[test]
    fn kmer_vec_test() {
        let com = CompositionComputer::new("../test_data/reads.fq", "../test_data/reads.kmers", 4);
        let kvec = vectorise_one(com.kcount, com.ksize, &com.pos_map, "AAAANGAGA", false);
        assert_eq!(kvec[0], 1.0);
    }

    #[test]
    fn kmer_vec_norm_test() {
        let com = CompositionComputer::new("../test_data/reads.fq", "../test_data/reads.kmers", 4);
        let kvec = vectorise_one(com.kcount, com.ksize, &com.pos_map, "AAAANGAGA", true);
        assert_eq!(kvec[0], 0.5);
    }

    #[test]
    fn vec_mmap_test() {
        let com = CompositionComputer::new("../test_data/reads.fq", "../test_data/reads.kmers", 4);
        let kvec = com.vectorise_mmap();
        // assert_eq!(kvec[0], 0.5);
    }
}
