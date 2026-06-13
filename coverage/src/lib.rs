use counter::CountComputer;
use kmer::{kmer::KmerGenerator, KmerU1024, KmerU2048, KmerU256, KmerU512, KmerWord};
use ktio::seq::{SeqFormat, Sequences};
use rayon::prelude::*;
use std::{
    cmp::min,
    collections::HashMap,
    fmt::Debug,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    str::FromStr,
};

const NUMBER_SIZE: usize = 8;

trait CoverageKmer: KmerWord + FromStr
where
    <Self as FromStr>::Err: Debug,
{
}

impl<K> CoverageKmer for K
where
    K: KmerWord + FromStr,
    <K as FromStr>::Err: Debug,
{
}

pub struct CovComputer {
    in_path: String,
    in_path_kmer: String,
    out_dir: String,
    ksize: usize,
    threads: usize,
    norm: bool,
    delim: String,
    bin_size: usize,
    bin_count: usize,
    memory_ceil_gb: f64,
}

impl CovComputer {
    pub fn new(
        in_path: String,
        out_dir: String,
        ksize: usize,
        bin_size: usize,
        bin_count: usize,
    ) -> Self {
        Self {
            in_path: in_path.clone(),
            in_path_kmer: in_path,
            out_dir,
            ksize,
            threads: rayon::current_num_threads(),
            norm: true,
            delim: " ".to_owned(),
            bin_size,
            bin_count,
            memory_ceil_gb: 6_f64,
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

    pub fn set_kmer_path(&mut self, path: String) {
        self.in_path_kmer = path;
    }

    pub fn set_max_memory(&mut self, memory_ceil_gb: f64) {
        self.memory_ceil_gb = memory_ceil_gb;
    }

    pub fn build_table(&self) -> Result<(), String> {
        let mut ctr =
            CountComputer::new(self.in_path_kmer.clone(), self.out_dir.clone(), self.ksize);
        ctr.set_threads(self.threads);
        ctr.set_max_memory(self.memory_ceil_gb);
        ctr.count();
        ctr.merge(true);
        Ok(())
    }

    pub fn compute_coverages(&self) {
        match self.ksize {
            1..=32 => self.compute_coverages_with::<u64>(),
            33..=64 => self.compute_coverages_with::<u128>(),
            65..=128 => self.compute_coverages_with::<KmerU256>(),
            129..=256 => self.compute_coverages_with::<KmerU512>(),
            257..=512 => self.compute_coverages_with::<KmerU1024>(),
            513..=1024 => self.compute_coverages_with::<KmerU2048>(),
            _ => panic!("k-mer size must be between 1 and 1024 bases"),
        }
    }

    fn compute_coverages_with<K: CoverageKmer>(&self)
    where
        <K as FromStr>::Err: Debug,
    {
        let kmer_path = format!("{}/kmers.counts", self.out_dir);
        let vec_path = format!("{}/kmers.vectors", self.out_dir);
        let file = fs::File::open(kmer_path).unwrap();
        let buff = BufReader::new(file);
        let mut counts = HashMap::new();

        for line in buff.lines().map_while(Result::ok) {
            let mut parts = line.trim().split('\t');
            let kmer: K = parts.next().unwrap().parse().unwrap();
            let count: u32 = parts.next().unwrap().parse().unwrap();
            counts.insert(kmer, count);
        }

        let reader = ktio::seq::get_reader(&self.in_path).unwrap();
        let format = SeqFormat::get(&self.in_path).unwrap();
        let records = Sequences::new(format, reader).unwrap();
        let file = File::create(vec_path).unwrap();
        let mut out_buffer = BufWriter::new(file);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();

        pool.install(|| {
            rayon::scope(|_| {
                let mut buffer = Vec::with_capacity(1000);
                let mut total = 0_usize;

                for record in records {
                    total += record.seq.len();
                    buffer.push(record);

                    if total as u64 >= self.memory_ceil_gb as u64 * (1 << 30) {
                        let result = buffer
                            .par_iter()
                            .map(|seq| {
                                let kvec = self.vectorise_one(&seq.seq, &counts);
                                // optimise this with pre-sized string
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
                        buffer.clear();
                        total = 0;
                    }
                }

                if total > 0 {
                    // optimise this with pre-sized string
                    let result = buffer
                        .par_iter()
                        .map(|seq| {
                            let kvec = self.vectorise_one(&seq.seq, &counts);
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
                    buffer.clear();
                }
            });
        });
    }

    fn vectorise_one<K: KmerWord>(&self, seq: &[u8], counts: &HashMap<K, u32>) -> Vec<f64> {
        let mut vec = vec![0_f64; self.bin_count];
        let mut total = 0_f64;

        for (fmer, rmer) in KmerGenerator::<K>::new(seq, self.ksize) {
            let min_mer = min(fmer, rmer);
            let count = *counts.get(&min_mer).unwrap_or(&0);
            let kmer_bin = (count as f64 / self.bin_size as f64).floor() as usize;
            let vec_bin = min(kmer_bin, self.bin_count - 1);
            unsafe {
                // we already know the size of the vector and
                *vec.get_unchecked_mut(vec_bin) += 1_f64;
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
    use ktio::fops::create_directory;
    use std::fs;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn kmer_count_vecs_test() {
        create_directory("../test_data/computed_coverages").expect("Directory must be creatable");
        let mut cov = CovComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_coverages".to_owned(),
            4,
            2,
            3,
        );
        cov.build_table().unwrap();
        cov.memory_ceil_gb = 0.1;
        cov.compute_coverages();

        assert_eq!(
            fs::read("../test_data/expected_counts.vectors").unwrap(),
            fs::read("../test_data/computed_coverages/kmers.vectors").unwrap()
        );

        cov.memory_ceil_gb = 1.0;
        cov.compute_coverages();

        assert_eq!(
            fs::read("../test_data/expected_counts.vectors").unwrap(),
            fs::read("../test_data/computed_coverages/kmers.vectors").unwrap()
        );
    }

    #[test]
    fn kmer_count_vecs_unnorm_test() {
        create_directory("../test_data/computed_coverage_unnorm")
            .expect("Directory must be creatable");
        let mut cov = CovComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_coverage_unnorm".to_owned(),
            4,
            2,
            3,
        );
        cov.set_norm(false);
        cov.build_table().unwrap();
        cov.compute_coverages();

        assert_eq!(
            fs::read("../test_data/expected_counts_unnorm.vectors").unwrap(),
            fs::read("../test_data/computed_coverage_unnorm/kmers.vectors").unwrap()
        );
    }

    #[test]
    fn coverage_supports_three_limb_kmers() {
        create_directory("../test_data/computed_coverage_unnorm_wide")
            .expect("Directory must be creatable");
        let mut cov = CovComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_coverage_unnorm_wide".to_owned(),
            50,
            2,
            3,
        );
        cov.set_norm(false);
        cov.build_table().unwrap();
        cov.compute_coverages();

        assert_eq!(
            fs::read("../test_data/expected_counts_unnorm_wide.vectors").unwrap(),
            fs::read("../test_data/computed_coverage_unnorm_wide/kmers.vectors").unwrap()
        );
    }
}
