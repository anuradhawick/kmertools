use indicatif::{ProgressBar, ProgressStyle};
use kmer::{kmer::KmerGenerator, numeric_to_kmer, Kmer};
use ktio::{
    fops::delete_file_if_exists,
    seq::{get_reader, SeqFormat, Sequences},
};
use rayon::prelude::*;
use scc::HashMap as SccMap;
use std::{
    cmp::{max, min},
    fs,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc, Mutex,
    },
};

// only to make code more readable
type SeqArc = Arc<Mutex<Sequences<BufReader<Box<dyn Read + Sync + Send>>>>>;

pub struct CountComputer {
    in_path: String,
    out_dir: String,
    ksize: usize,
    threads: usize,
    records: SeqArc,
    chunks: u64,
    n_parts: u64,
    memory_ceil_gb: f64,
    seq_count: u64,
    debug: bool,
    acgt: bool,
}

impl CountComputer {
    pub fn new(in_path: String, out_dir: String, ksize: usize) -> Self {
        let format = SeqFormat::get(&in_path).unwrap();
        let reader = ktio::seq::get_reader(&in_path).unwrap();
        let records = Sequences::new(format, reader).unwrap();

        Self {
            in_path,
            out_dir,
            ksize,
            threads: rayon::current_num_threads(),
            records: Arc::new(Mutex::new(records)),
            chunks: 0,
            n_parts: 0,
            seq_count: 0,
            memory_ceil_gb: 6_f64,
            debug: false,
            acgt: false,
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn set_max_memory(&mut self, memory_ceil_gb: f64) {
        self.memory_ceil_gb = memory_ceil_gb;
    }

    pub fn set_acgt_output(&mut self, acgt: bool) {
        self.acgt = acgt;
    }

    pub fn count(&mut self) {
        self.init();
        let pbar = ProgressBar::new(self.seq_count);
        pbar.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}",
            )
            .unwrap()
            .progress_chars("#>-"),
        );
        loop {
            // TODO have to fix below line being called even the next chunk does not exist
            pbar.set_message(format!("Processing chunk: {}", self.chunks + 1));
            let records = self.count_chunk(&pbar);
            if records > 0 {
                self.chunks += 1;
            } else {
                break;
            }
        }
        pbar.finish();
    }

    fn count_chunk(&self, pbar: &ProgressBar) -> u64 {
        let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();
        let total_records = Arc::new(AtomicU64::new(0));
        // an estimate of worse case kmer count
        let total_kmers_so_far = Arc::new(AtomicU64::new(0));
        let counts_table: Vec<SccMap<Kmer, u32>> = vec![SccMap::new(); self.n_parts as usize];
        let counts_table_arc = Arc::new(counts_table);
        // make pbar for all bases struct wide

        pool.scope(|scope| {
            for _ in 0..self.threads {
                let records_arc_clone = Arc::clone(&self.records);
                let total_records_clone = Arc::clone(&total_records);
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
                            pbar.inc(1);
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
                        } else {
                            // end of iteration
                            break;
                        }
                    }
                });
            }
        });

        let recs = total_records.load(Ordering::Acquire);

        if recs == 0 {
            return 0;
        }

        pool.scope(|_| {
            counts_table_arc
                .par_iter()
                .enumerate()
                .for_each(|(part, map)| {
                    let outf = fs::File::create(format!(
                        "{}/temp_kmers.part_{}_chunk_{}",
                        self.out_dir, part, self.chunks
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

    pub fn merge(&self, delete: bool) {
        let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();
        let outf = fs::File::create(format!("{}/kmers.counts", self.out_dir)).unwrap();
        let mut buff = BufWriter::new(outf);
        let pbar = ProgressBar::new(self.n_parts * self.chunks);
        pbar.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}",
            )
            .unwrap()
            .progress_chars("#>-"),
        );

        for part in 0..self.n_parts {
            let completed = Arc::new(AtomicU64::new(0));
            let map: SccMap<Kmer, u32> = SccMap::new();
            let map_arc = Arc::new(map);
            pbar.set_message(format!("Merging partition: {}", part + 1));

            pool.scope(|scope| {
                for chunk in 0..self.chunks {
                    let map_arc_clone = Arc::clone(&map_arc);
                    let pbar_clone = pbar.clone();
                    let completed_clone = Arc::clone(&completed);

                    scope.spawn(move |_| {
                        let path =
                            format!("{}/temp_kmers.part_{}_chunk_{}", self.out_dir, part, chunk);
                        let file = fs::File::open(&path).unwrap();
                        let buff = BufReader::new(file);
                        for line in buff.lines().map_while(Result::ok) {
                            let mut parts = line.trim().split('\t');
                            let kmer: Kmer = parts.next().unwrap().parse().unwrap();
                            let count: u32 = parts.next().unwrap().parse().unwrap();
                            *map_arc_clone.entry(kmer).or_insert(0) += count;
                        }
                        if delete {
                            delete_file_if_exists(&path).expect("file must be removable");
                        }
                        completed_clone.fetch_add(1, Ordering::Acquire);
                        pbar_clone.inc(1);
                    });
                }
            });

            map_arc.scan(|k, v| {
                if self.acgt {
                    buff.write_all(
                        format!("{}\t{:?}\n", numeric_to_kmer(*k, self.ksize), v).as_bytes(),
                    )
                    .unwrap();
                } else {
                    buff.write_all(format!("{}\t{:?}\n", k, v).as_bytes())
                        .unwrap();
                }
            });
        }

        pbar.finish();
    }

    pub fn init(&mut self) {
        let reader = get_reader(&self.in_path).unwrap();
        let format = SeqFormat::get(&self.in_path).unwrap();
        let stats = Sequences::seq_stats(format, reader);
        let data_size_gb = stats.total_length as f64 / (1 << 30) as f64;
        // assuming 8 bytes per kmer
        // at least this should be the num threads for fastest possible merging
        let n_parts = max(
            if self.debug { 1 } else { self.threads as u64 },
            (8_f64 * data_size_gb / (2_f64 * self.memory_ceil_gb)).ceil() as u64,
        );
        self.n_parts = n_parts;
        self.seq_count = stats.seq_count as u64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ktio::fops::{create_directory, load_lines_sorted};

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn count_test() {
        create_directory("../test_data/computed_counts").expect("Directory must be creatable");
        let mut ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts".to_owned(),
            15,
        );
        ctr.debug = true;
        ctr.count();
        assert_eq!(ctr.n_parts, 1);
        assert_eq!(ctr.chunks, 1);
        let exp = load_lines_sorted("../test_data/expected_counts.part_0_chunk_0");
        let res = load_lines_sorted("../test_data/computed_counts/temp_kmers.part_0_chunk_0");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }

    #[test]
    fn merge_test() {
        let mut ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts_test".to_owned(),
            15,
        );
        ctr.chunks = 2;
        ctr.n_parts = 2;
        ctr.merge(false);
        let exp = load_lines_sorted("../test_data/expected_counts_test.counts");
        let res = load_lines_sorted("../test_data/computed_counts_test/kmers.counts");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }

    #[test]
    fn merge_acgt_test() {
        let mut ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts_acgt_test".to_owned(),
            15,
        );
        ctr.chunks = 2;
        ctr.n_parts = 2;
        ctr.set_acgt_output(true);
        ctr.merge(false);
        let exp = load_lines_sorted("../test_data/expected_counts_acgt_test.counts");
        let res = load_lines_sorted("../test_data/computed_counts_acgt_test/kmers.counts");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }
}
