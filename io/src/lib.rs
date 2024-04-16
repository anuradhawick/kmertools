use bio::io::fasta::{Reader as FastaReader, Records as FastaRecords};
use bio::io::fastq::{Reader as FastqReader, Records as FastqRecords};
use std::io::{BufRead, BufReader};

// Record set entries of type R, which implement BufRead trait (stdin/file)
pub enum RecordSet<R: BufRead> {
    Fasta(FastaRecords<BufReader<R>>),
    Fastq(FastqRecords<BufReader<R>>),
}

pub struct Sequences<R: BufRead> {
    pub current_record: usize,
    pub records: RecordSet<R>,
}

pub struct Sequence {
    pub n: usize,
    pub id: String,
    pub seq: String,
}

pub enum SeqFormat {
    Fasta,
    Fastq,
}

impl SeqFormat {
    pub fn get(path: &str) -> Option<SeqFormat> {
        if path.ends_with("fq") || path.ends_with("fastq") {
            return Some(SeqFormat::Fastq);
        } else if path.ends_with("fasta") || path.ends_with("fa") || path.ends_with("fna") {
            return Some(SeqFormat::Fasta);
        }
        None
    }
}

impl<R: BufRead> Sequences<R> {
    pub fn new(format: SeqFormat, reader: R) -> Result<Self, String> {
        match format {
            SeqFormat::Fastq => {
                // let file = fs::File::open(path).map_err(|_| "Error".to_string());
                let fastq_reader = FastqReader::new(reader);
                // let s = fastq_reader.records();
                Ok(Sequences {
                    current_record: 0,
                    records: RecordSet::Fastq(fastq_reader.records()),
                })
            }
            SeqFormat::Fasta => {
                // let file = fs::File::open(path).map_err(|_| "Error".to_string());
                let fasta_reader = FastaReader::new(reader);
                Ok(Sequences {
                    current_record: 0,
                    records: RecordSet::Fasta(fasta_reader.records()),
                })
            }
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
                        seq: String::from_utf8(record.seq().to_owned()).unwrap(),
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
                        seq: String::from_utf8(record.seq().to_owned()).unwrap(),
                    });
                }
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn load_fq_file_test() {
        let file = fs::File::open("../test_data/reads.fq").unwrap();
        let reader = BufReader::new(file);
        let mut seqs = Sequences::new(SeqFormat::Fastq, reader).unwrap();
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
        let file = fs::File::open("../test_data/reads.fa").unwrap();
        let reader = BufReader::new(file);
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader).unwrap();
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
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!("ACGTACGTACGT", record_1.seq);
        let finish = seqs.next();
        assert!(finish.is_none());
    }
}
