use ktio::seq::{SeqFormat, Sequence, Sequences};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufWriter, Write},
};

const GB_4: usize = 4 * (1 << 30);
type Point = (f64, f64);

// Code and test adopted from https://github.com/skatila/pycgr (as of 2:55 am Friday, 7 June 2024 Coordinated Universal Time (UTC))
// Git tag 922ebd4bec482ae6522453c14d3f9dc5d1c99995
// Under full compliance of GPL-3.0 license (https://github.com/skatila/pycgr/blob/922ebd4bec482ae6522453c14d3f9dc5d1c99995/LICENSE)
pub struct CgrComputer {
    in_path: String,
    out_path: String,
    threads: usize,
    memory: usize,
    cgr_center: Point,
    cgr_map: HashMap<u8, Point>,
}

impl CgrComputer {
    pub fn new(in_path: String, out_path: String, vecsize: usize) -> Self {
        let (cgr_center, cgr_map) = CgrComputer::cgr_maps(vecsize as f64);
        Self {
            in_path,
            out_path,
            threads: rayon::current_num_threads(),
            memory: GB_4,
            cgr_center,
            cgr_map,
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn vectorise(&self) -> Result<(), String> {
        let mut reader = ktio::seq::get_reader(&self.in_path).unwrap();
        let buffer = reader
            .fill_buf()
            .map_err(|_| String::from("Invalid stream"))?;
        let format = if buffer[0] == b'>' {
            SeqFormat::Fasta
        } else {
            SeqFormat::Fastq
        };
        let records = Sequences::new(format, reader).unwrap();
        let file = File::create(&self.out_path)
            .map_err(|_| format!("Unable to write to file: {}", self.out_path))?;
        let mut out_buffer = BufWriter::new(file);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();

        pool.install(|| {
            rayon::scope(|_| {
                let mut buffer = Vec::with_capacity(1000);
                let mut total = 0_usize;

                // Define a closure to handle buffer processing
                let mut process_buffer = |buffer: &Vec<Sequence>| {
                    let result = buffer
                        .par_iter()
                        .map(|seq| {
                            let kvec = self.vectorise_one(&seq.seq).unwrap();
                            let kvec_str: Vec<String> = kvec
                                .iter()
                                .map(|val| format!("({},{})", val.0, val.1))
                                .collect();
                            format!("{}\n", kvec_str.join(" "))
                        })
                        .collect::<Vec<String>>()
                        .join("");
                    out_buffer.write_all(result.as_bytes()).unwrap();
                };

                for record in records {
                    total += record.seq.len();
                    buffer.push(record);

                    if total >= self.memory {
                        process_buffer(&buffer);
                        buffer.clear();
                        total = 0;
                    }
                }

                if !buffer.is_empty() {
                    process_buffer(&buffer);
                }
            });
        });
        Ok(())
    }

    fn vectorise_one(&self, seq: &[u8]) -> Result<Vec<Point>, String> {
        let mut cgr = Vec::with_capacity(seq.len());
        let mut cgr_marker = self.cgr_center;

        for s in seq.iter() {
            if let Some(&cgr_corner) = self.cgr_map.get(s) {
                cgr_marker = (
                    (cgr_corner.0 + cgr_marker.0) / 2.0,
                    (cgr_corner.1 + cgr_marker.1) / 2.0,
                );
                cgr.push(cgr_marker);
            } else {
                return Err("Bad nucleotide, unable to proceed".to_string());
            }
        }

        Ok(cgr)
    }

    fn cgr_maps(vecsize: f64) -> (Point, HashMap<u8, Point>) {
        let cgr_a: Point = (0.0, 0.0);
        let cgr_t: Point = (vecsize, 0.0);
        let cgr_g: Point = (vecsize, vecsize);
        let cgr_c: Point = (0.0, vecsize);
        let cgr_center: Point = (vecsize / 2.0, vecsize / 2.0);

        let cgr_dict: HashMap<u8, Point> = [
            (b'A', cgr_a), // Adenine
            (b'T', cgr_t), // Thymine
            (b'G', cgr_g), // Guanine
            (b'C', cgr_c), // Cytosine
            (b'U', cgr_t), // Uracil (demethylated form of thymine)
            (b'a', cgr_a), // Adenine
            (b't', cgr_t), // Thymine
            (b'g', cgr_g), // Guanine
            (b'c', cgr_c), // Cytosine
            (b'u', cgr_t), // Uracil/Thymine
        ]
        .iter()
        .cloned()
        .collect();

        (cgr_center, cgr_dict)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn cgr_vec_test() {
        let cgr = CgrComputer::new(PATH_FQ.to_owned(), "".to_owned(), 1);
        let vec = cgr
            .vectorise_one("atgatgaaatagagagactttat".as_bytes())
            .unwrap();
        let res = vec![
            (0.25, 0.25),
            (0.625, 0.125),
            (0.8125, 0.5625),
            (0.40625, 0.28125),
            (0.703125, 0.140625),
            (0.8515625, 0.5703125),
            (0.42578125, 0.28515625),
            (0.212890625, 0.142578125),
            (0.1064453125, 0.0712890625),
            (0.55322265625, 0.03564453125),
            (0.276611328125, 0.017822265625),
            (0.6383056640625, 0.5089111328125),
            (0.31915283203125, 0.25445556640625),
            (0.659576416015625, 0.627227783203125),
            (0.3297882080078125, 0.3136138916015625),
            (0.6648941040039062, 0.6568069458007812),
            (0.3324470520019531, 0.3284034729003906),
            (0.16622352600097656, 0.6642017364501953),
            (0.5831117630004883, 0.33210086822509766),
            (0.7915558815002441, 0.16605043411254883),
            (0.8957779407501221, 0.08302521705627441),
            (0.44788897037506104, 0.04151260852813721),
            (0.7239444851875305, 0.020756304264068604),
        ];

        assert_eq!(vec, res);
    }

    #[test]
    fn cgr_complete_unnorm_test() {
        let mut cgr = CgrComputer::new(PATH_FQ.to_owned(), "../test_data/reads.cgr".to_owned(), 1);
        cgr.set_threads(4);
        cgr.vectorise().unwrap();

        assert_eq!(
            fs::read("../test_data/expected_reads.cgr").unwrap(),
            fs::read("../test_data/reads.cgr").unwrap()
        )
    }
}
