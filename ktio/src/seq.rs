use bio::io::fasta::{Reader as FastaReader, Records as FastaRecords};
use bio::io::fastq::{Reader as FastqReader, Records as FastqRecords};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

// Record set entries of type R, which implement BufRead trait (stdin/file)
pub enum RecordSet<R: BufRead> {
    Fasta(FastaRecords<BufReader<R>>),
    Fastq(FastqRecords<BufReader<R>>),
}

pub struct Sequence {
    pub n: usize,
    pub id: String,
    pub seq: Vec<u8>,
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
    pub records: RecordSet<R>,
}

impl<R: BufRead> Sequences<R> {
    pub fn new(format: SeqFormat, reader: R) -> Result<Self, String> {
        match format {
            SeqFormat::Fastq => {
                let fastq_reader = FastqReader::new(reader);
                Ok(Sequences {
                    current_record: 0,
                    records: RecordSet::Fastq(fastq_reader.records()),
                })
            }
            SeqFormat::Fasta => {
                let fasta_reader = FastaReader::new(reader);
                Ok(Sequences {
                    current_record: 0,
                    records: RecordSet::Fasta(fasta_reader.records()),
                })
            }
        }
    }

    pub fn seq_stats(format: SeqFormat, reader: R) -> SeqStats {
        let mut total_length = 0_usize;
        let mut seq_count = 0_usize;

        match format {
            SeqFormat::Fastq => {
                let fastq_reader = FastqReader::new(reader);
                for record in fastq_reader.records() {
                    total_length += record.unwrap().seq().len();
                    seq_count += 1;
                }
            }
            SeqFormat::Fasta => {
                let fasta_reader = FastaReader::new(reader);
                for record in fasta_reader.records() {
                    total_length += record.unwrap().seq().len();
                    seq_count += 1;
                }
            }
        }

        SeqStats {
            seq_count,
            total_length,
        }
    }
}

impl<R: BufRead> Iterator for Sequences<R> {
    type Item = Sequence;

    fn next(&mut self) -> Option<Self::Item> {
        // records do not have a common trait to get id and seq, we can create one
        // but this looks simpler for the time being
        match self.records {
            RecordSet::Fastq(ref mut records) => {
                let next_record = records.next();
                if let Some(record) = next_record {
                    let record = record.unwrap();
                    self.current_record += 1;
                    return Some(Sequence {
                        n: self.current_record - 1,
                        id: record.id().to_string(),
                        seq: record.seq().to_vec(),
                    });
                }
                None
            }
            RecordSet::Fasta(ref mut records) => {
                let next_record = records.next();
                if let Some(record) = next_record {
                    let record = record.unwrap();
                    self.current_record += 1;
                    return Some(Sequence {
                        n: self.current_record - 1,
                        id: record.id().to_string(),
                        seq: record.seq().to_vec(),
                    });
                }
                None
            }
        }
    }

    fn count(self) -> usize
    where
        Self: Sized,
    {
        unimplemented!("Count cannot be performed without always having a rewindable input stream, stdin is not!");
    }
}

pub fn get_reader(path: &str) -> Result<BufReader<Box<dyn Read + Sync + Send>>, String> {
    if path == "-" {
        let stdin = io::stdin();
        Ok(BufReader::new(Box::new(stdin)))
    } else {
        let is_zip = path.ends_with(".gz");
        let file = File::open(path).map_err(|_| format!("Unable to open: {}", path))?;
        if is_zip {
            let decoder = flate2::read::GzDecoder::new(file);
            Ok(BufReader::new(Box::new(decoder)))
        } else {
            Ok(BufReader::new(Box::new(file)))
        }
    }
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
        let mut seqs = Sequences::new(SeqFormat::Fastq, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Read_1", record_1.id);
        assert_eq!(
            b"GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA".to_vec(),
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Read_2", record_2.id);
        assert_eq!(
            b"GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT".to_vec(),
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn load_fa_file_test() {
        let reader = get_reader(PATH_FA).unwrap();
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!(
            b"GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA".to_vec(),
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Record_2", record_2.id);
        assert_eq!(
            b"GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT".to_vec(),
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn load_fa_stdin_test() {
        let input = ">Record_1\nACGTACGTACGT";
        let reader = BufReader::new(input.as_bytes());
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!(b"ACGTACGTACGT".to_vec(), record_1.seq);
        let finish = seqs.next();
        assert!(finish.is_none());
    }
}
