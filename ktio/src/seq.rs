use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

use crate::io::FastxReader;

#[derive(PartialEq, Debug)]
pub struct Sequence {
    pub n: usize, // 0 based index of the sequence
    pub id: String,
    pub desc: String,
    pub seq: String,
    pub qual: String,
}

impl Sequence {
    /// Create a new sequence record
    pub fn new(id: String, desc: String, seq: String, qual: String) -> Self {
        Self {
            n: 0,
            id,
            desc,
            seq,
            qual,
        }
    }
}

pub struct SeqStats {
    pub seq_count: usize,
    pub total_length: usize,
}

#[derive(Debug, Clone, Copy)]
pub enum SeqFormat {
    Fasta,
    Fastq,
}

impl SeqFormat {
    pub fn get(path: &str) -> Option<SeqFormat> {
        let mut path = path;
        if path.ends_with(".gz") {
            path = path.trim_end_matches(".gz");
        }
        if path.ends_with(".fq") || path.ends_with(".fastq") {
            return Some(SeqFormat::Fastq);
        } else if path.ends_with(".fasta") || path.ends_with(".fa") || path.ends_with(".fna") {
            return Some(SeqFormat::Fasta);
        }
        None
    }
}

pub struct Sequences<R: BufRead> {
    pub current_record: usize,
    pub records: FastxReader<R>,
}

impl<R: BufRead> Sequences<R> {
    pub fn new(format: SeqFormat, reader: R) -> Self {
        Self {
            current_record: 0,
            records: FastxReader::new(format, reader),
        }
    }

    pub fn seq_stats(format: SeqFormat, reader: R) -> SeqStats {
        let mut total_length = 0_usize;
        let mut seq_count = 0_usize;
        let mut records = FastxReader::new(format, reader);

        while let Ok(record) = records.read_next()
            && let Some(record) = record
        {
            total_length += record.seq.len();
            seq_count += 1;
        }

        SeqStats {
            seq_count,
            total_length,
        }
    }

    pub fn read_next(&mut self) -> io::Result<Option<Sequence>> {
        let mut res = self.records.read_next();
        if let Ok(Some(seq)) = res.as_mut() {
            seq.n = self.current_record;
        }
        self.current_record += 1;
        res
    }
}

impl<R: BufRead> Iterator for Sequences<R> {
    type Item = Sequence;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next().expect("Error reading sequence")
    }
}

pub fn get_reader(path: &str) -> Result<BufReader<Box<dyn Read + Sync + Send>>, String> {
    if path == "-" {
        let stdin = io::stdin();
        return Ok(BufReader::new(Box::new(stdin)));
    }
    let file = File::open(path).map_err(|e| format!("Unable to open {}: {}", path, e))?;
    let input: Box<dyn Read + Sync + Send> = if path.ends_with(".gz") {
        // Gzip members include BGZF blocks; decode all members, not just the first.
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::with_capacity(64 * 1024, input))
}

#[cfg(test)]
mod tests {
    use super::*;
    const PATH_FQ: &str = "../test_data/reads.fq";
    const PATH_FA: &str = "../test_data/reads.fa";
    const PATH_FQ_GZ: &str = "../test_data/reads.fq.gz";

    #[test]
    fn seq_stats_test() {
        // fastq
        let reader = get_reader(PATH_FQ).unwrap();
        let stats = Sequences::seq_stats(SeqFormat::Fastq, reader);
        assert_eq!(stats.seq_count, 2);
        assert_eq!(stats.total_length, 144);
        // fasta
        let reader = get_reader(PATH_FA).unwrap();
        let stats = Sequences::seq_stats(SeqFormat::Fasta, reader);
        assert_eq!(stats.seq_count, 2);
        assert_eq!(stats.total_length, 144);
        // gz
        let reader = get_reader(PATH_FQ_GZ).unwrap();
        let stats = Sequences::seq_stats(SeqFormat::Fastq, reader);
        assert_eq!(stats.seq_count, 2);
        assert_eq!(stats.total_length, 144);
    }

    #[test]
    fn load_fq_file_test() {
        let reader = get_reader(PATH_FQ).unwrap();
        let mut seqs = Sequences::new(SeqFormat::Fastq, reader);
        let record_1 = seqs.next().unwrap();
        assert_eq!("Read_1", record_1.id);
        assert_eq!(
            "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Read_2", record_2.id);
        assert_eq!(
            "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT",
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn load_fa_file_test() {
        let reader = get_reader(PATH_FA).unwrap();
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader);
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!(
            "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Record_2", record_2.id);
        assert_eq!(
            "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT",
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn load_fa_stdin_test() {
        let input = ">Record_1\nACGTACGTACGT";
        let reader = BufReader::new(input.as_bytes());
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader);
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!("ACGTACGTACGT", record_1.seq);
        let finish = seqs.next();
        assert!(finish.is_none());
    }
}
