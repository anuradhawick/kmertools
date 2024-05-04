use ktio::seq::{get_reader, SeqFormat, Sequences};

pub struct CountComputer {
    in_path: String,
    out_path: String,
    ksize: usize,
    threads: usize,
}

impl CountComputer {
    pub fn new(in_path: String, out_path: String, ksize: usize) -> Self {
        Self {
            in_path,
            out_path,
            ksize,
            threads: rayon::current_num_threads(),
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn count(&self) {
        let format = SeqFormat::get(&self.in_path).unwrap();
        let reader = get_reader(&self.in_path).unwrap();
        let stats = Sequences::seq_stats(format, reader);
        let k_mer_count_estimate =
            stats.total_length as i64 + stats.seq_count as i64 * (1 - self.ksize as i64);
        let n_iterations = (((k_mer_count_estimate as f64)
            * 2_f64.powf((2_f64 * (self.ksize as f64)).log2().ceil())
            + 32_f64)
            / (4_f64 * (1 << 30) as f64))
            .ceil();

        println!(
            "Count = {} Length = {}",
            stats.seq_count, stats.total_length
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn count_test() {
        let ctr = CountComputer::new(PATH_FQ.to_owned(), "../test_data/counts".to_owned(), 15);
        ctr.count();
        // let result = add(2, 2);
        // assert_eq!(result, 4);
    }
}
