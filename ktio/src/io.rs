//! Buffered FASTA/FASTQ parsing inspired by klib's kseq.h:
//! https://github.com/attractivechaos/klib/blob/master/kseq.h
//!
//! Read into owned record buffers, concatenate sequence lines, and consume quality by
//! sequence length. This is an independent Rust implementation with strict
//! format validation rather than kseq's permissive header scanning.
use crate::seq::{SeqFormat, Sequence};
use std::io::{self, BufRead};

pub struct FastxReader<R> {
    reader: R,
    format: SeqFormat,
    line: Vec<u8>,
    seq_capacity: usize,
    qual_capacity: usize,
    pending_header: bool,
    finished: bool,
}

impl<R: BufRead> FastxReader<R> {
    pub fn new(format: SeqFormat, reader: R) -> Self {
        Self {
            reader,
            format,
            line: Vec::with_capacity(256),
            seq_capacity: 256,
            qual_capacity: 256,
            pending_header: false,
            finished: false,
        }
    }

    // BufRead::read_until scans buffered chunks for LF and copies slices in
    // bulk. No per-byte I/O, UTF-8 decoding, or per-line allocation is needed.
    fn read_line(&mut self) -> io::Result<bool> {
        self.line.clear();
        Self::append_line(&mut self.reader, &mut self.line)
    }

    // Append directly to the record's eventual String allocation, then remove
    // only this line's terminator. Empty lines must not trim earlier data.
    fn append_line(reader: &mut R, bytes: &mut Vec<u8>) -> io::Result<bool> {
        let start = bytes.len();
        if reader.read_until(b'\n', bytes)? == 0 {
            return Ok(false);
        }
        if bytes.last() == Some(&b'\n') {
            bytes.pop();
            if bytes.len() > start && bytes.last() == Some(&b'\r') {
                bytes.pop();
            }
        }
        Ok(true)
    }

    // Valid UTF-8 transfers ownership without copying. Preserve the existing
    // replacement-character behavior for malformed UTF-8 input.
    fn into_string(bytes: Vec<u8>) -> String {
        String::from_utf8(bytes)
            .unwrap_or_else(|error| String::from_utf8_lossy(error.as_bytes()).into_owned())
    }

    fn invalid(message: &str) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidData, message)
    }

    /// Read one record. EOF and errors are terminal; errors never produce a
    /// partially parsed record. Sequence/quality buffers become the returned
    /// strings; the previous record's size guides the next allocation.
    pub fn read_next(&mut self) -> io::Result<Option<Sequence>> {
        if self.finished {
            return Ok(None);
        }
        let result = self.read_record();
        if !matches!(result, Ok(Some(_))) {
            self.finished = true;
        }
        result
    }

    fn read_record(&mut self) -> io::Result<Option<Sequence>> {
        if self.pending_header {
            self.pending_header = false;
        } else {
            loop {
                if !self.read_line()? {
                    return Ok(None);
                }
                if !self.line.is_empty() {
                    break;
                }
            }
        }
        let marker = match self.format {
            SeqFormat::Fasta => b'>',
            SeqFormat::Fastq => b'@',
        };
        if self.line.first() != Some(&marker) {
            return Err(Self::invalid("Expected a FASTA/FASTQ header"));
        }
        let header = &self.line[1..];
        let split = header
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(header.len());
        if split == 0 {
            return Err(Self::invalid("Sequence identifier is empty"));
        }
        let id = String::from_utf8_lossy(&header[..split]).into_owned();
        let desc = if split < header.len() {
            String::from_utf8_lossy(&header[split + 1..]).into_owned()
        } else {
            String::new()
        };
        let mut seq = Vec::with_capacity(self.seq_capacity);
        let mut qual = Vec::new();
        match self.format {
            SeqFormat::Fasta => loop {
                let start = seq.len();
                if !Self::append_line(&mut self.reader, &mut seq)? {
                    break;
                }
                if seq.get(start) == Some(&b'>') {
                    self.line.clear();
                    self.line.extend_from_slice(&seq[start..]);
                    seq.truncate(start);
                    self.pending_header = true;
                    break;
                }
            },
            SeqFormat::Fastq => {
                loop {
                    let start = seq.len();
                    if !Self::append_line(&mut self.reader, &mut seq)? {
                        return Err(Self::invalid("FASTQ record is missing its '+' separator"));
                    }
                    if seq.get(start) == Some(&b'+') {
                        seq.truncate(start);
                        break;
                    }
                    if seq.get(start) == Some(&b'@') {
                        return Err(Self::invalid("FASTQ record is missing its '+' separator"));
                    }
                }
                // Read at least one quality line, including for empty records.
                qual.reserve(self.qual_capacity.max(seq.len().saturating_add(2)));
                loop {
                    if !Self::append_line(&mut self.reader, &mut qual)? {
                        return Err(Self::invalid("FASTQ quality string is truncated"));
                    }
                    if qual.len() >= seq.len() {
                        break;
                    }
                }
                if qual.len() != seq.len() {
                    return Err(Self::invalid("FASTQ sequence and quality lengths differ"));
                }
            }
        }
        // Leave room for a lookahead header or '+' separator without doubling
        // the next record's allocation when its sequence has the same length.
        self.seq_capacity = seq.len().saturating_add(256);
        self.qual_capacity = qual.len().saturating_add(2);
        Ok(Some(Sequence {
            n: 0,
            id,
            desc,
            seq: Self::into_string(seq),
            qual: Self::into_string(qual),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn parse(format: SeqFormat, input: &[u8], capacity: usize) -> io::Result<Vec<Sequence>> {
        let mut reader = FastxReader::new(format, BufReader::with_capacity(capacity, input));
        let mut records = Vec::new();
        while let Some(record) = reader.read_next()? {
            records.push(record);
        }
        assert!(reader.read_next()?.is_none());
        Ok(records)
    }

    #[test]
    fn wrapped_records_and_buffer_boundaries() {
        for capacity in [1, 2, 7, 64 * 1024] {
            let fasta = parse(
                SeqFormat::Fasta,
                b"\n>one\tdescription with spaces\r\nAC\r\n\r\nGT\r\n>empty\r\n>last\nN",
                capacity,
            )
            .unwrap();
            assert_eq!(fasta.len(), 3);
            assert_eq!(
                fasta[0],
                Sequence::new(
                    "one".into(),
                    "description with spaces".into(),
                    "ACGT".into(),
                    "".into()
                )
            );
            assert!(fasta[1].seq.is_empty());
            assert_eq!(fasta[2].seq, "N");
            let fastq = parse(
                SeqFormat::Fastq,
                b"@one desc\r\nAC\r\nGT\r\n+one\r\n@+\r\n!I\r\n@two\nN\n+\n#",
                capacity,
            )
            .unwrap();
            assert_eq!(fastq.len(), 2);
            assert_eq!(fastq[0].seq, "ACGT");
            assert_eq!(fastq[0].qual, "@+!I");
            assert_eq!(fastq[1].qual, "#");
        }
    }

    #[test]
    fn long_lines_and_owned_records() {
        let bases = "ACGT".repeat(100_000);
        let input = format!(">long\n{}\n>short\nN\n", bases);
        let records = parse(SeqFormat::Fasta, input.as_bytes(), 17).unwrap();
        assert_eq!(records[0].seq, bases);
        assert_eq!(records[1].seq, "N");
    }

    #[test]
    fn direct_buffers_preserve_utf8_replacement_and_bare_carriage_returns() {
        for capacity in [1, 7, 64 * 1024] {
            let records =
                parse(SeqFormat::Fasta, b">first\nAC\xff\n>second\nN\r", capacity).unwrap();
            assert_eq!(records[0].seq, "AC\u{fffd}");
            assert_eq!(records[1].seq, "N\r");

            let records = parse(
                SeqFormat::Fastq,
                b"@first\nA\xff\n+first\n!\xfe\n@second\nT\n+\n@",
                capacity,
            )
            .unwrap();
            assert_eq!(records[0].seq, "A\u{fffd}");
            assert_eq!(records[0].qual, "!\u{fffd}");
            assert_eq!(records[1].qual, "@");
        }
    }

    #[test]
    fn empty_input_and_empty_fastq_record() {
        assert!(parse(SeqFormat::Fasta, b"", 1).unwrap().is_empty());
        let records = parse(SeqFormat::Fastq, b"@empty\n\n+\n\n@next\nA\n+\n!\n", 1).unwrap();
        assert!(records[0].seq.is_empty());
        assert!(records[0].qual.is_empty());
        assert_eq!(records[1].seq, "A");
    }

    #[test]
    fn malformed_records_fail_once() {
        for input in [
            b"bad\n".as_slice(),
            b"@\nA\n+\n!\n",
            b"@x\nAC\n",
            b"@x\nAC\n@next\n",
            b"@x\nAC\n+\n!",
            b"@x\nA\n+\n!!\n",
            b"@x\nA\n+",
        ] {
            let mut reader = FastxReader::new(SeqFormat::Fastq, Cursor::new(input));
            assert_eq!(
                reader.read_next().unwrap_err().kind(),
                io::ErrorKind::InvalidData
            );
            assert!(reader.read_next().unwrap().is_none());
        }
        assert!(parse(SeqFormat::Fasta, b">\nAC\n", 1).is_err());
    }
}
