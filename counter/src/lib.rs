use dashmap::DashMap as DMap;
use indicatif::ProgressBar;
use kmer::{kmer::KmerGenerator, kmer_minimisers::KmerMinimiserGenerator, minimiser, Kmer};
use ktio::seq::{SeqFormat, Sequences};
use scc::HashMap as SccMap;
use std::{
    cmp::min,
    collections::HashMap,
    fs,
    hash::Hash,
    io::{BufReader, BufWriter, Read, Write},
    sync::{
        atomic::{AtomicU64, Ordering},
        mpsc::{channel, sync_channel},
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
    table_part: Arc<AtomicU64>,
}

impl CountComputer {
    pub fn new(in_path: String, out_path: String, ksize: usize) -> Self {
        let format = SeqFormat::get(&in_path).unwrap();
        let reader = ktio::seq::get_reader(&in_path).unwrap();
        let records = Sequences::new(format, reader).unwrap();

        Self {
            in_path,
            out_path,
            ksize,
            threads: rayon::current_num_threads(),
            records: Arc::new(Mutex::new(records)),
            table_part: Arc::new(AtomicU64::new(0)),
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
        let counts_table: SccMap<Kmer, u32> = SccMap::with_capacity(5 * (1 << 20));
        let counts_table_arc = Arc::new(counts_table);

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
                        if total_kmers_so_far_clone.load(Ordering::Relaxed) > 2 * (1 << 30) {
                            break;
                        }
                        let record = { records_arc_clone.lock().unwrap().next() };
                        if let Some(record) = record {
                            total_records_clone.fetch_add(1, Ordering::Relaxed);
                            for (fmer, rmer) in KmerGenerator::new(&record.seq, self.ksize) {
                                let min_mer = min(fmer, rmer);
                                // counts_table_arc_clone
                                //     .entry(min_mer)
                                //     .and_modify(|v| *v += 1)
                                //     .or_insert(1);
                            }
                            total_kmers_so_far_clone
                                .fetch_add(record.seq.len() as u64, Ordering::Relaxed);
                            if total_records_clone.load(Ordering::Acquire) % 10000 == 0 {
                                pbar_clone.set_message(format!(
                                    "Processed no. of sequences from part: {}: {}",
                                    self.table_part.load(Ordering::Relaxed),
                                    record.n + 1
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

        pbar.set_message(format!(
            "Processed no. of sequences from part: {}: {}",
            self.table_part.load(Ordering::Relaxed),
            total_records.load(Ordering::Acquire)
        ));
        pbar.finish();

        let outf = fs::File::create(format!(
            "{}.part_{}",
            self.out_path,
            self.table_part.load(Ordering::Acquire)
        ))
        .unwrap();
        let mut buff = BufWriter::new(outf);
        self.table_part.fetch_add(1, Ordering::SeqCst);

        counts_table_arc.scan(|k, v| {
            buff.write_all(format!("{k}\t{v:?}\n").as_bytes()).unwrap();
        });

        counts_table_arc.scan(|k, v| {
            buff.write_all(format!("{k}\t{v:?}\n").as_bytes()).unwrap();
        });

        total_records.load(Ordering::Acquire)
    }

    pub fn count_min(&self) -> u64 {
        let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();
        let total_records = Arc::new(AtomicU64::new(0));
        // an estimate of worse case kmer count
        let total_kmers_so_far = Arc::new(AtomicU64::new(0));
        let pbar = ProgressBar::new_spinner();
        let counts_table: SccMap<Kmer, SccMap<Kmer, u32>> = SccMap::with_capacity(5 * (1 << 20));
        let counts_table_arc = Arc::new(counts_table);
        let minimiser_tables: Vec<SccMap<Kmer, u32>> = vec![SccMap::new(); 1048576];
        let minimiser_tables_arc = Arc::new(minimiser_tables);

        pool.scope(|scope| {
            for _ in 0..self.threads {
                let records_arc_clone = Arc::clone(&self.records);
                let total_records_clone = Arc::clone(&total_records);
                let pbar_clone = pbar.clone();
                let counts_table_arc_clone = Arc::clone(&counts_table_arc);
                let total_kmers_so_far_clone = Arc::clone(&total_kmers_so_far);
                let minimiser_tables_arc_clone = Arc::clone(&minimiser_tables_arc);

                scope.spawn(move |_| {
                    loop {
                        // when limit reached exit without further reads
                        if total_kmers_so_far_clone.load(Ordering::Relaxed) > 2 * (1 << 30) {
                            break;
                        }
                        let record = { records_arc_clone.lock().unwrap().next() };
                        if let Some(record) = record {
                            total_records_clone.fetch_add(1, Ordering::Relaxed);
                            for (mmer, _, _, kmers) in
                                KmerMinimiserGenerator::new(&record.seq, self.ksize, 10)
                            {
                                unsafe {
                                    let map =
                                        minimiser_tables_arc_clone.get_unchecked(mmer as usize);
                                    for k in kmers {
                                        map.entry(k).and_modify(|v| *v += 1).or_insert(1);
                                    }
                                }
                            }
                            total_kmers_so_far_clone
                                .fetch_add(record.seq.len() as u64, Ordering::Relaxed);
                            if total_records_clone.load(Ordering::Acquire) % 10000 == 0 {
                                pbar_clone.set_message(format!(
                                    "Processed no. of sequences from part: {}: {}",
                                    self.table_part.load(Ordering::Relaxed),
                                    record.n + 1
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

        pbar.set_message(format!(
            "Processed no. of sequences from part: {}: {}",
            self.table_part.load(Ordering::Relaxed),
            total_records.load(Ordering::Acquire)
        ));
        pbar.finish();

        let outf = fs::File::create(format!(
            "{}.part_{}",
            self.out_path,
            self.table_part.load(Ordering::Acquire)
        ))
        .unwrap();
        let mut buff = BufWriter::new(outf);
        self.table_part.fetch_add(1, Ordering::SeqCst);

        for map in minimiser_tables_arc.iter() {
            map.scan(|k, v| {
                buff.write_all(format!("{k}\t{v:?}\n").as_bytes()).unwrap();
            });
        }
        total_records.load(Ordering::Acquire)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // const PATH_FQ: &str = "../test_data/reads.fq";
    const PATH_FQ: &str = "/home/anuvini/Downloads/reads.fasta";

    #[test]
    fn count_test() {
        let ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts".to_owned(),
            15,
        );
        loop {
            let seqs = ctr.count();
            // let seqs = ctr.count_min();
            if seqs == 0 {
                break;
            }
            println!(
                "Finished part {}. with {} sequences",
                ctr.table_part.load(Ordering::Acquire),
                seqs
            );
        }
        // let result = add(2, 2);
        // assert_eq!(result, 4);
    }
}
