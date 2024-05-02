use kmer::{kmer::KmerGenerator, Kmer};
use ktio::seq::*;
use scc::HashMap as SccMap;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    sync::{mpsc, Arc, Mutex},
};

const MAX_MEM: u64 = 5000 * (1_u64 << 20);

struct DSKCounter {
    in_path: String,
    out_path: String,
    ksize: f64,
    total_seq_count: f64,
    total_seq_length: f64,
    threads: usize,
    // algorithm
    v: u64,
    n_iters: u64,
    n_parts: u64,
    // temp files
    temp_file_names: Vec<String>,
    temp_files: SccMap<Kmer, BufWriter<File>>,
}

impl DSKCounter {
    pub fn new(in_path: String, out_path: String, ksize: usize) -> Self {
        let format = SeqFormat::get(&in_path).unwrap();
        let reader = get_reader(&in_path).unwrap();
        let stats = Sequences::seq_stats(format, reader);
        let v = (stats.total_length + stats.seq_count - ksize * stats.seq_count) as u64;
        Self {
            in_path,
            ksize: ksize as f64,
            total_seq_count: stats.seq_count as f64,
            total_seq_length: stats.total_length as f64,
            v,
            n_iters: 1,
            threads: rayon::current_num_threads(),
            out_path,
            n_parts: Self::parts(v as f64, ksize as f64, 1.0, MAX_MEM as f64),
            temp_file_names: Vec::new(),
            temp_files: SccMap::new(),
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    fn parts(v: f64, k: f64, _n_iters: f64, m: f64) -> u64 {
        let log_2k = (2.0 * k).log2().ceil();
        // 32 bits may be common in hashmaps, but scc may be different
        // investigate further
        let numerator = v * 2.0_f64.powf(log_2k) + 32_f64;
        let denominator = 0.7_f64 * m;

        (numerator / denominator).ceil() as u64
    }

    fn part_to_fname(&self, part: &str) -> String {
        format!("{}.{}.{}", self.out_path, part, "tmp")
    }

    fn write_to_buckets(&self) -> Result<(), String> {
        let format = SeqFormat::get(&self.in_path).ok_or("Unsupported input format")?;
        let outfile = File::create(&self.out_path).unwrap();
        let mut out_writer = BufWriter::new(outfile);

        println!("No parts: {}", self.n_parts);

        let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();

        pool.scope(|scope| {
            let format = SeqFormat::get(&self.in_path).unwrap();
            let reader = ktio::seq::get_reader(&self.in_path).unwrap();
            let records = Sequences::new(format, reader).unwrap();
            let records_arc = Arc::new(Mutex::new(records));

            for _ in 0..1 {
                let records_arc_clone = Arc::clone(&records_arc);
                scope.spawn(move |_| {
                    let mut table = HashMap::new();
                    let mut r = 0;
                    loop {
                        let record = { records_arc_clone.lock().unwrap().next() };

                        if let Some(record) = record {
                            for (fmer, rmer) in KmerGenerator::new(&record.seq, self.ksize as usize)
                            {
                                let min_mer = Kmer::min(fmer, rmer);
                                table.entry(min_mer).and_modify(|v| *v += 1).or_insert(1);
                            }
                            r += 1;
                            // break;
                            if table.len() > 10_000_000 {
                                println!(
                                    "Flush a hash table of size: {} counted reads {}",
                                    table.len(),
                                    r
                                );
                                break;
                            }
                        } else {
                            // end of iteration
                            break;
                        }
                    }
                });
            }
        });

        // for part in 0..1 {
        //     // let fname = self.part_to_fname(&part.to_string());
        //     // let file = File::create(&fname).unwrap();
        //     // let mut writer = BufWriter::with_capacity(1 << 20, file);
        //     let reader = ktio::seq::get_reader(&self.in_path).unwrap();
        //     let records = Sequences::new(format, reader).unwrap();

        //     for record in records {
        //         for (fmer, rmer) in KmerGenerator::new(&record.seq, self.ksize as usize) {
        //             // let min_mer = Kmer::min(fmer, rmer);
        //             // let kmer_part = min_mer % self.n_parts;

        //             // if part == kmer_part {
        //             //     // writer.write_all(&min_mer.to_ne_bytes()).unwrap();
        //             // }
        //         }
        //     }

        //     // writer.flush().unwrap();
        //     // drop(writer);

        //     // let file = File::open(&fname).expect("This file must exist");
        //     // let mut reader = BufReader::new(file);
        //     // let mut buffer = [0u8; 8]; // Buffer to hold chunks of 8 bytes
        //     //                            // let counts = SccMap::new();
        //     // let mut counts = HashMap::new();

        //     // while reader.read_exact(&mut buffer).is_ok() {
        //     //     let number = u64::from_le_bytes(buffer);
        //     //     counts
        //     //         .entry(number)
        //     //         .and_modify(|val| *val += 1)
        //     //         .or_insert(1_u64);
        //     // }

        //     // counts.scan(|k, v| {
        //     //     let line = format!("{} {}\n", k, v);
        //     //     out_writer.write_all(line.as_bytes()).unwrap();
        //     //     //
        //     // });

        //     // ktio::fops::delete_file_if_exists(fname).unwrap();
        // }

        Ok(())
    }

    fn aggregate(&self) {
        let outfile = File::create(&self.out_path).unwrap();
        let mut writer = BufWriter::new(outfile);

        // aggregate each file
        for fname in &self.temp_file_names {
            let file = File::open(fname).expect("This file must exist");
            let mut reader = BufReader::new(file);
            let mut buffer = [0u8; 8]; // Buffer to hold chunks of 8 bytes
            let counts = SccMap::new();

            while reader.read_exact(&mut buffer).is_ok() {
                let number = u64::from_le_bytes(buffer);
                counts
                    .entry(number)
                    .and_modify(|val| *val += 1)
                    .or_insert(1_u64);
            }

            counts.scan(|k, v| {
                let str = format!("{} {}\n", k, v);
                writer.write_all(str.as_bytes()).unwrap();
                //
            });

            ktio::fops::delete_file_if_exists(fname).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // const PATH_FQ: &str = "../test_data/reads.fq";
    const PATH_FQ: &str = "/home/anuvini/Downloads/reads.fasta";

    #[test]
    fn kmer_counts_test() {
        let mut ctr = DSKCounter::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts.kmers".to_owned(),
            15,
        );
        println!("Partitions: {}", ctr.n_parts);
        println!("Total kmers: {}", ctr.v);
        ctr.write_to_buckets().unwrap();
        // ctr.aggregate();
    }
}
