use indicatif::ProgressBar;
use kmer::{minimiser::MinimiserGenerator, numeric_to_kmer};
use ktio::seq::*;
use scc::HashMap as SccMap;
use std::{
    fs,
    io::{BufReader, BufWriter, Read, Write},
    sync::{atomic::AtomicU64, Arc, Mutex},
};

pub fn bin_sequences(wsize: usize, msize: usize, in_path: &str, out_path: &str, threads: usize) {
    let mut threads = threads;
    if threads == 0 {
        threads = rayon::current_num_threads();
    }
    let format = SeqFormat::get(in_path).unwrap();
    let reader = ktio::seq::get_reader(in_path).unwrap();
    let records: Sequences<BufReader<Box<dyn Read + Sync + Send>>> =
        Sequences::new(format, reader).unwrap();
    let pbar = ProgressBar::new_spinner();
    let result: SccMap<String, Vec<(String, usize, usize)>> = SccMap::new();
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let records_arc = Arc::new(Mutex::new(records));
    let result_arc = Arc::new(result);
    let total_records = Arc::new(AtomicU64::new(0));

    pool.scope(|scope| {
        for _ in 0..threads {
            let records_arc_clone = Arc::clone(&records_arc);
            let result_arc_clone = Arc::clone(&result_arc);
            let total_records_clone = Arc::clone(&total_records);
            let pbar_clone = pbar.clone();

            scope.spawn(move |_| {
                loop {
                    let record = {
                        total_records_clone.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        records_arc_clone.lock().unwrap().next()
                    };
                    if let Some(record) = record {
                        let mgen = if wsize == 0 {
                            MinimiserGenerator::new(&record.seq, record.seq.len(), msize)
                        } else {
                            MinimiserGenerator::new(&record.seq, wsize, msize)
                        };
                        for (k, s, e) in mgen {
                            result_arc_clone
                                .entry_sync(numeric_to_kmer(k, msize))
                                .and_modify(|v| v.push((record.id.clone(), s, e)))
                                .or_insert(vec![(record.id.clone(), s, e)]);
                        }
                        // buff.write_all("\n".as_bytes()).unwrap();
                        if (record.n + 1) % 10000 == 0 {
                            pbar_clone.set_message(format!(
                                "Processed no. of sequences: {}",
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
        "Processed no. of sequences: {}",
        total_records.load(std::sync::atomic::Ordering::Acquire)
    ));
    pbar.finish();

    let outf = fs::File::create(out_path).unwrap();
    let mut buff = BufWriter::new(outf);

    result_arc.iter_sync(|k, v| {
        buff.write_all(format!("{k}\t{v:?}\n").as_bytes()).unwrap();
        true
    });
}

pub fn seq_to_min(wsize: usize, msize: usize, in_path: &str, out_path: &str, threads: usize) {
    let mut threads = threads;
    if threads == 0 {
        threads = rayon::current_num_threads();
    }
    let format = SeqFormat::get(in_path).unwrap();
    let reader = ktio::seq::get_reader(in_path).unwrap();
    let records: Sequences<BufReader<Box<dyn Read + Sync + Send>>> =
        Sequences::new(format, reader).unwrap();
    let pbar = ProgressBar::new_spinner();
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let records_arc = Arc::new(Mutex::new(records));
    let total_records = Arc::new(AtomicU64::new(0));
    let outf = fs::File::create(out_path).unwrap();
    let buff = Arc::new(Mutex::new(BufWriter::new(outf)));

    pool.scope(|scope| {
        for _ in 0..threads {
            let records_arc_clone = Arc::clone(&records_arc);
            let total_records_clone = Arc::clone(&total_records);
            let pbar_clone = pbar.clone();
            let buff_clone = Arc::clone(&buff);

            scope.spawn(move |_| {
                loop {
                    let record = {
                        total_records_clone.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        records_arc_clone.lock().unwrap().next()
                    };
                    if let Some(record) = record {
                        let mgen = if wsize == 0 {
                            MinimiserGenerator::new(&record.seq, record.seq.len(), msize)
                        } else {
                            MinimiserGenerator::new(&record.seq, wsize, msize)
                        };
                        let mut mins = Vec::new();
                        mins.push(record.id);

                        for (k, s, e) in mgen {
                            mins.push(format!("{}:{}-{}", numeric_to_kmer(k, msize), s, e));
                        }
                        mins.push("\n".to_string());
                        {
                            buff_clone
                                .lock()
                                .unwrap()
                                .write_all(mins.join("\t").as_bytes())
                                .unwrap();
                        }
                        if (record.n + 1) % 10000 == 0 {
                            pbar_clone.set_message(format!(
                                "Processed no. of sequences: {}",
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
        "Processed no. of sequences: {}",
        total_records.load(std::sync::atomic::Ordering::Acquire)
    ));
    pbar.finish();
}

#[cfg(test)]
mod tests {
    use super::*;
    use ktio::fops::load_lines_sorted;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn bin_sequences_test() {
        bin_sequences(0, 10, PATH_FQ, "../test_data/computed_minimisers", 32);
        let exp = load_lines_sorted("../test_data/expected_minimisers");
        let res = load_lines_sorted("../test_data/computed_minimisers");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }

    #[test]
    fn seq_to_min_test() {
        seq_to_min(31, 7, PATH_FQ, "../test_data/computed_seq_minimisers", 32);
        let exp = load_lines_sorted("../test_data/expected_seq_minimisers");
        let res = load_lines_sorted("../test_data/computed_seq_minimisers");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }
}
