use bio::io::fasta::{Reader as FastaReader, Records as FastaRecords};
use bio::io::fastq::{Reader as FastqReader, Records as FastqRecords};
use std::fmt::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader};

enum RecordSet<R: BufRead> {
    Fasta(FastaRecords<R>),
    Fastq(FastqRecords<R>),
}

struct Seq<'a, R: BufRead> {
    path: &'a str,
    pub records: RecordSet<R>,
}

impl<'a> Seq<'a, BufReader<File>> {
    fn new(path: &'a str) -> Result<Self, String> {
        if path.ends_with("fq") || path.ends_with("fastq") {
            let file = fs::File::open(path).map_err(|_| "Error".to_string());
            let fastq_reader = FastqReader::new(file?);
            // let s = fastq_reader.records();
            Ok(Seq {
                path,
                records: RecordSet::Fastq(fastq_reader.records()),
            })
        } else if path.ends_with("fa") || path.ends_with("fasta") {
            let file = fs::File::open(path).map_err(|_| "Error".to_string());
            let fasta_reader = FastaReader::new(file?);
            Ok(Seq {
                path,
                records: RecordSet::Fasta(fasta_reader.records()),
            })
        } else {
            Err("Unsupported file format".to_string())
        }
    }
}

// impl<'a> Iterator for Seq<'a, R: BufRead> {
//     type Item = ;
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_fq() {
        let seqs = Seq::new("../test_data/reads.fq").unwrap();
        let records = seqs.records;

        match records {
            RecordSet::Fastq(records) => {
                for record in records {
                    println!("{:?}", record);
                }
            }
            _ => {}
        }

        // for record in seqs {}
    }
}
