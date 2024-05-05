use super::Kmer;
use std::cmp::min;
use std::collections::VecDeque;
use std::iter::Iterator;

// https://github.com/lh3/minimap2/blob/0cc3cdca27f050fb80a19c90d25ecc6ab0b0907b/sketch.c#L9C1-L26C3
const SEQ_NT4_TABLE: [u8; 256] = [
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];
const REV_MASK: u64 = 3;

pub struct MinimiserGenerator<'a> {
    seq: &'a [u8],
    pos: usize,
    wsize: usize,
    msize: usize,
    m_mask: u64,
    m_window_start: usize,
    m_window_end: usize,
    m_val_f: u64,
    m_val_r: u64,
    m_val_l: usize,
    m_active: u64,
    m_shift: u64,
    buff: VecDeque<u64>,
    buff_pos: usize,
}

impl<'a> MinimiserGenerator<'a> {
    pub fn new(seq: &'a [u8], wsize: usize, msize: usize) -> Self {
        MinimiserGenerator {
            seq,
            wsize,
            msize,
            pos: 0,
            buff_pos: 0,
            m_active: u64::MAX,
            m_mask: (1_u64 << (2 * msize)) - 1,
            m_val_f: 0,
            m_val_r: 0,
            m_val_l: 0,
            m_window_end: 0,
            m_window_start: 0,
            buff: VecDeque::with_capacity(wsize - msize + 1),
            m_shift: 2 * (msize - 1) as u64,
        }
    }
}

// technique adopted from https://github.com/lh3/minimap2/blob/0cc3cdca27f050fb80a19c90d25ecc6ab0b0907b/sketch.c#L77
impl<'a> Iterator for MinimiserGenerator<'a> {
    type Item = (Kmer, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let mut min_m_val: u64;
        let mut prev_m_val: u64;
        let mut prev_w_start: usize;
        let mut prev_w_end: usize;

        loop {
            if self.pos == self.seq.len() {
                return None;
            }
            let pos_char = self.seq[self.pos];
            let pos_f_val = SEQ_NT4_TABLE[pos_char as usize] as u64;
            let pos_r_val = pos_f_val ^ REV_MASK;

            if pos_f_val < 4 {
                // non ambiguous
                // minimiser
                self.m_val_f = ((self.m_val_f << 2) | pos_f_val) & self.m_mask;
                self.m_val_r = (self.m_val_r >> 2) | (pos_r_val << self.m_shift);
                self.m_val_l += 1;
            } else {
                // ambiguous
                // check if we have passed a good complete window
                let should_return = self.buff.len() == self.wsize - self.msize + 1;
                prev_m_val = self.m_active;
                prev_w_start = self.m_window_start;
                prev_w_end = self.pos;
                self.buff_pos = 0;
                self.m_active = u64::MAX;
                self.m_val_f = 0;
                self.m_val_r = 0;
                self.m_val_l = 0;
                self.m_window_end = 0;
                self.m_window_start = self.pos + 1;

                self.buff.clear();
                self.pos += 1;
                if should_return {
                    return Some((prev_m_val, prev_w_start, prev_w_end));
                }
                continue;
            }

            if self.m_val_l < self.msize {
                self.pos += 1;
                continue;
            }

            self.m_val_l -= 1;
            // self.w_val_l -= 1;
            min_m_val = min(self.m_val_f, self.m_val_r);

            // minimiser buffer is full
            if self.buff.len() == self.wsize - self.msize + 1 {
                // pop and push to avoid resizing
                self.buff.pop_front();
                self.buff.push_back(min_m_val);

                // we have removed the minimum
                if self.buff_pos == 0 {
                    let mut new_min = u64::MAX;
                    for j in 0..self.buff.len() {
                        if *self.buff.get(j).unwrap() < new_min {
                            self.buff_pos = j;
                            new_min = *self.buff.get(j).unwrap();
                        }
                    }
                    // minimiser changed
                    if new_min != self.m_active {
                        self.m_window_end = self.pos;
                        prev_m_val = self.m_active;
                        prev_w_start = self.m_window_start;
                        prev_w_end = self.m_window_end;
                        self.m_active = new_min;
                        self.m_window_start = self.pos - self.wsize + 1;
                        self.pos += 1;
                        return Some((prev_m_val, prev_w_start, prev_w_end));
                    }
                } else if min_m_val < self.m_active {
                    // break the window
                    self.m_window_end = self.pos;
                    prev_m_val = self.m_active;
                    prev_w_start = self.m_window_start;
                    prev_w_end = self.m_window_end;
                    self.m_active = min_m_val;
                    self.buff_pos = self.buff.len() - 1;
                    self.m_window_start = self.pos - self.wsize + 1;
                    self.pos += 1;
                    return Some((prev_m_val, prev_w_start, prev_w_end));
                } else {
                    self.buff_pos -= 1;
                }
            } else {
                // add new minimizer to buffer
                self.buff.push_back(min_m_val);
            }

            // first time we are experiencing all minimizers
            if self.m_active == u64::MAX && self.buff.len() == self.wsize - self.msize + 1 {
                for j in 0..self.buff.len() {
                    if *self.buff.get(j).unwrap() < self.m_active {
                        self.buff_pos = j;
                        self.m_active = *self.buff.get(j).unwrap();
                    }
                }
            }

            if self.pos == self.seq.len() - 1 {
                self.pos += 1;
                return Some((self.m_active, self.m_window_start, self.seq.len()));
            }

            self.pos += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numeric_to_kmer;

    #[test]
    fn minimisers_generated_test() {
        // Acquired from https://homolog.us/blogs/bioinfo/2017/10/25/intro-minimizer/
        let mut mg = MinimiserGenerator::new(b"ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA", 31, 7);
        let (kmer, start, end) = mg.next().unwrap();
        assert_eq!(numeric_to_kmer(kmer, 7), "ACGATAT");
        assert_eq!(
            &mg.seq[start..end],
            "ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCG".as_bytes()
        );
        println!(
            "Window start:{:5} end:{:5} - {:60} : {:10}",
            start,
            end,
            String::from_utf8(mg.seq[start..end].to_vec()).unwrap(),
            numeric_to_kmer(kmer, 7)
        );
        let (kmer, start, end) = mg.next().unwrap();
        assert_eq!(numeric_to_kmer(kmer, 7), "ACGCCTA");
        assert_eq!(
            &mg.seq[start..end],
            "TATCGTAGGCGTCGATGGAGAGCTAGATCGATCGAT".as_bytes()
        );
        println!(
            "Window start:{:5} end:{:5} - {:60} : {:10}",
            start,
            end,
            String::from_utf8(mg.seq[start..end].to_vec()).unwrap(),
            numeric_to_kmer(kmer, 7)
        );
        let (kmer, start, end) = mg.next().unwrap();
        assert_eq!(numeric_to_kmer(kmer, 7), "AGAGCTA");
        assert_eq!(
            &mg.seq[start..end],
            "AGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCC".as_bytes()
        );
        println!(
            "Window start:{:5} end:{:5} - {:60} : {:10}",
            start,
            end,
            String::from_utf8(mg.seq[start..end].to_vec()).unwrap(),
            numeric_to_kmer(kmer, 7)
        );
        let (kmer, start, end) = mg.next().unwrap();
        assert_eq!(numeric_to_kmer(kmer, 7), "AAATCCC");
        assert_eq!(
            &mg.seq[start..end],
            "ATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAG".as_bytes()
        );
        println!(
            "Window start:{:5} end:{:5} - {:60} : {:10}",
            start,
            end,
            String::from_utf8(mg.seq[start..end].to_vec()).unwrap(),
            numeric_to_kmer(kmer, 7)
        );
        let (kmer, start, end) = mg.next().unwrap();
        assert_eq!(numeric_to_kmer(kmer, 7), "AATCCCG");
        assert_eq!(
            &mg.seq[start..end],
            "AATCCCGATCGATTCCGAGCGCGATCAAAGC".as_bytes()
        );
        println!(
            "Window start:{:5} end:{:5} - {:60} : {:10}",
            start,
            end,
            String::from_utf8(mg.seq[start..end].to_vec()).unwrap(),
            numeric_to_kmer(kmer, 7)
        );
        let (kmer, start, end) = mg.next().unwrap();
        assert_eq!(numeric_to_kmer(kmer, 7), "AATCGAT");
        assert_eq!(
            &mg.seq[start..end],
            "ATCCCGATCGATTCCGAGCGCGATCAAAGCG".as_bytes()
        );
        println!(
            "Window start:{:5} end:{:5} - {:60} : {:10}",
            start,
            end,
            String::from_utf8(mg.seq[start..end].to_vec()).unwrap(),
            numeric_to_kmer(kmer, 7)
        );
        let (kmer, start, end) = mg.next().unwrap();
        assert_eq!(numeric_to_kmer(kmer, 7), "AAAGCGC");
        assert_eq!(
            &mg.seq[start..end],
            "TCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA".as_bytes()
        );
        println!(
            "Window start:{:5} end:{:5} - {:60} : {:10}",
            start,
            end,
            String::from_utf8(mg.seq[start..end].to_vec()).unwrap(),
            numeric_to_kmer(kmer, 7)
        );
        let res = mg.next();
        assert_eq!(res, None);
    }

    #[test]
    fn minimisers_generated_with_error_test() {
        // Acquired from https://homolog.us/blogs/bioinfo/2017/10/25/intro-minimizer/
        let mg = MinimiserGenerator::new(b"ATGCGATATCGNTAGGCGTCGATGGA", 8, 5);
        let seq = mg.seq;
        let expected = [
            ("ATGCGATA", "ATCGC"),
            ("TGCGATATCG", "ATATC"),
            ("TAGGCGTCGA", "ACGCC"),
            ("GCGTCGATGGA", "ATCGA"),
        ];

        for (i, (kmer, start, end)) in mg.enumerate() {
            assert_eq!(&seq[start..end], expected[i].0.as_bytes());
            assert_eq!(numeric_to_kmer(kmer, 5), expected[i].1);
            println!(
                "Window start:{:5} end:{:5} - {:60} : {:10}",
                start,
                end,
                String::from_utf8(seq[start..end].to_vec()).unwrap(),
                numeric_to_kmer(kmer, 5)
            );
        }
    }
}
