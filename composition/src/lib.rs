use kmer::KmerGenerator;
use memmap2::{MmapMut, MmapOptions};
use std::fs;
use std::fs::OpenOptions;
use std::io::BufReader;
use std::io::Error;
use std::sync::Mutex;

const NUMBER_SIZE: usize = 8;
const THREADS: usize = 12;

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
        let seqs = io::Sequences::new(io::SeqFormat::get(self.in_path).unwrap(), reader).unwrap();
        seqs.count()
    }

    fn create_writable_mmap(&self) -> Result<MmapMut, Error> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .truncate(true)
            .create(true)
            .open(self.out_path)?;

        let per_line_size = self.kcount * (NUMBER_SIZE + 1);
        let estimated_file_size = self.scount * per_line_size;
        file.set_len(estimated_file_size as u64)?;
        let mmap = unsafe { MmapOptions::new().len(estimated_file_size).map_mut(&file)? };

        Ok(mmap)
    }

    fn vectorise_mmap(&self) -> Result<(), String> {
        let mut mmap = self.create_writable_mmap().unwrap();
        let file = fs::File::open(self.in_path).unwrap();
        let reader = BufReader::new(file);
        let lock = Mutex::new(0);
        let mut records =
            io::Sequences::new(io::SeqFormat::get(self.in_path).unwrap(), reader).unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(THREADS)
            .build()
            .unwrap();
        let mut record_index: usize = 0;

        for _ in 0..THREADS {
            pool.install(|| loop {
                let (record, record_index) = {
                    let _guard = lock.lock().unwrap();
                    let record = records.next();
                    let thread_record_index = record_index;
                    record_index += 1;
                    (record, thread_record_index)
                };

                if let Some(record) = record {
                    let kvec =
                        vectorise_one(self.kcount, self.ksize, &self.pos_map, &record.seq, true);
                    let kvec_str: Vec<String> = kvec
                        .iter()
                        .map(|val| format!("{:.*}", NUMBER_SIZE - 2, val))
                        .collect();
                    let kvec_str = format!("{}\n", kvec_str.join(" "));
                    let start_pos = kvec_str.len() * record_index;
                    let pos = &mut mmap[start_pos..start_pos + kvec_str.bytes().len()];
                    pos.copy_from_slice(kvec_str.as_bytes());
                    // println!("{}", record_index);
                } else {
                    break;
                }
            });
        }

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
        println!("{:?}", com.scount);
        let _ = com.vectorise_mmap();
    }
}
