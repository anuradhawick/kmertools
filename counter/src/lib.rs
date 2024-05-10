use indicatif::{ProgressBar, ProgressStyle};
use kmer::{kmer::KmerGenerator, kmer_minimisers::KmerMinimiserGenerator, minimiser, Kmer};
use ktio::seq::{get_reader, SeqFormat, Sequences};
use rayon::prelude::*;
use scc::HashMap as SccMap;
use std::{
    cmp::{max, min},
    collections::HashMap,
    fs,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc, Mutex,
    },
};

// not very useful to make this a globally accepted type
type SeqArc = Arc<Mutex<Sequences<BufReader<Box<dyn Read + Sync + Send>>>>>;

pub struct CountComputer {
    in_path: String,
    out_path: String,
    ksize: usize,
    threads: usize,
    records: SeqArc,
    chunk: Arc<AtomicU64>,
    n_parts: u64,
    memory_ceil_gb: f64,
}

impl CountComputer {
    pub fn new(in_path: String, out_path: String, ksize: usize) -> Self {
        let format = SeqFormat::get(&in_path).unwrap();
        let reader = ktio::seq::get_reader(&in_path).unwrap();
        let records = Sequences::new(format, reader).unwrap();
        let memory_ceil_gb = 6_f64;
        let reader = get_reader(&in_path).unwrap();
        let format = SeqFormat::get(&in_path).unwrap();
        let stats = Sequences::seq_stats(format, reader);
        let data_size_gb = stats.total_length as f64 / (1 << 30) as f64;
        // assuming 8 bytes per kmer
        let n_parts = max(
            1,
            (8_f64 * data_size_gb / (2_f64 * memory_ceil_gb)).ceil() as u64,
        );

        Self {
            in_path,
            out_path,
            ksize,
            threads: rayon::current_num_threads(),
            records: Arc::new(Mutex::new(records)),
            chunk: Arc::new(AtomicU64::new(0)),
            n_parts,
            memory_ceil_gb,
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn count(&self) -> u64 {
        let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();
        let total_records = Arc::new(AtomicU64::new(0));
        // an estimate of worse case kmer count
        let total_kmers_so_far = Arc::new(AtomicU64::new(0));
        let pbar = ProgressBar::new_spinner();
        let counts_table: Vec<SccMap<Kmer, u32>> = vec![SccMap::new(); self.n_parts as usize];
        let counts_table_arc = Arc::new(counts_table);
        // make pbar for all bases struct wide

        pool.scope(|scope| {
            for _ in 0..self.threads {
                let records_arc_clone = Arc::clone(&self.records);
                let total_records_clone = Arc::clone(&total_records);
                let pbar_clone = pbar.clone();
                let counts_table_arc_clone = Arc::clone(&counts_table_arc);
                let total_kmers_so_far_clone = Arc::clone(&total_kmers_so_far);

                scope.spawn(move |_| {
                    loop {
                        // when limit reached exit without further reads
                        if total_kmers_so_far_clone.load(Ordering::Relaxed)
                            > (1_000_000_000_f64 * self.memory_ceil_gb / 8.0) as u64
                        {
                            break;
                        }
                        let record = { records_arc_clone.lock().unwrap().next() };
                        if let Some(record) = record {
                            total_records_clone.fetch_add(1, Ordering::Acquire);
                            for (fmer, rmer) in KmerGenerator::new(&record.seq, self.ksize) {
                                let min_mer = min(fmer, rmer);
                                unsafe {
                                    counts_table_arc_clone
                                        .get_unchecked((min_mer % self.n_parts) as usize)
                                        .entry(min_mer)
                                        .and_modify(|v| *v += 1)
                                        .or_insert(1);
                                }
                            }

                            total_kmers_so_far_clone
                                .fetch_add(record.seq.len() as u64, Ordering::Relaxed);
                            let recs = total_records_clone.load(Ordering::Acquire);
                            if recs % 10000 == 0 {
                                pbar_clone.set_message(format!(
                                    "Processed no. of sequences from chunk: {}: {}",
                                    self.chunk.load(Ordering::Relaxed),
                                    recs
                                ));
                                pbar_clone.tick();
                            }
                        } else {
                            // end of iteration
                            break;
                        }
                    }
                });
            }
        });

        let chunk = self.chunk.load(Ordering::Acquire);
        let recs = total_records.load(Ordering::Acquire);
        self.chunk.fetch_add(1, Ordering::Acquire);
        pbar.set_message(format!(
            "Processed no. of sequences from chunk: {}: {}",
            chunk, recs,
        ));
        pbar.finish();

        pool.scope(|_| {
            counts_table_arc
                .par_iter()
                .enumerate()
                .for_each(|(part, map)| {
                    let outf = fs::File::create(format!(
                        "{}.part_{}_chunk_{}",
                        self.out_path, part, chunk
                    ))
                    .unwrap();
                    let mut buff = BufWriter::new(outf);
                    map.scan(|k, v| {
                        buff.write_all(format!("{}\t{:?}\n", k, v).as_bytes())
                            .unwrap();
                    });
                })
        });

        total_records.load(Ordering::Acquire)
    }

    pub fn merge(&self) {
        let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();
        let chunks = 27;
        let pbar = ProgressBar::new(self.n_parts * chunks);
        let outf = fs::File::create(format!("{}.counts", self.out_path)).unwrap();
        let mut buff = BufWriter::new(outf);
        pbar.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
            )
            .unwrap()
            .progress_chars("#>-"),
        );

        for part in 0..self.n_parts {
            let completed = Arc::new(AtomicU64::new(0));
            let map: SccMap<Kmer, u32> = SccMap::new();
            let map_arc = Arc::new(map);
            pbar.set_message(format!("Merging partition: {}", part));

            pool.scope(|scope| {
                for chunk in 0..27 {
                    let map_arc_clone = Arc::clone(&map_arc);
                    let pbar_clone = pbar.clone();
                    let completed_clone = Arc::clone(&completed);

                    scope.spawn(move |_| {
                        let file = fs::File::open(format!(
                            "{}.part_{}_chunk_{}",
                            self.out_path, part, chunk
                        ))
                        .unwrap();
                        let buff = BufReader::new(file);
                        for line in buff.lines().map_while(Result::ok) {
                            let mut parts = line.trim().split('\t');
                            let kmer: Kmer = parts.next().unwrap().parse().unwrap();
                            let count: u32 = parts.next().unwrap().parse().unwrap();
                            *map_arc_clone.entry(kmer).or_insert(0) += count;
                        }
                        completed_clone.fetch_add(1, Ordering::Acquire);
                        pbar_clone.inc(1);
                    });
                }
            });

            map_arc.scan(|k, v| {
                buff.write_all(format!("{}\t{:?}\n", k, v).as_bytes())
                    .unwrap();
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // const PATH_FQ: &str = "../test_data/reads.fq";
    const PATH_FQ: &str = "/Users/wic053/Downloads/reads.fasta";

    #[test]
    fn count_test() {
        // for p in 0..10 {
        let ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts".to_owned(),
            15,
        );
        println!("n_parts: {}", ctr.n_parts);
        // println!("n_parts: {}", ctr.s);
        // loop {
        //     let seqs = ctr.count();
        //     // let seqs = ctr.count_min();
        //     if seqs == 0 {
        //         break;
        //     }
        //     println!(
        //         "Finished part {}. with {} sequences",
        //         ctr.chunk.load(Ordering::Acquire),
        //         seqs
        //     );
        // }
        println!("n_parts: {}", ctr.chunk.load(Ordering::Acquire));
    }

    #[test]
    fn merge_test() {
        let ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts".to_owned(),
            15,
        );
        println!("n_parts: {}", ctr.n_parts);
        ctr.merge()
    }
}
