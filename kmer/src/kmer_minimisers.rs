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

pub struct KmerMinimiserGenerator<'a> {
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
    k_mask: u64,
    k_val_f: u64,
    k_val_r: u64,
    k_val_l: usize,
    k_shift: u64,
    buff: VecDeque<u64>,
    buff_pos: usize,
}

impl<'a> KmerMinimiserGenerator<'a> {
    pub fn new(seq: &'a [u8], wsize: usize, msize: usize) -> Self {
        KmerMinimiserGenerator {
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
            k_mask: (1_u64 << (2 * wsize)) - 1,
            k_val_f: 0,
            k_val_r: 0,
            k_val_l: 0,
            m_window_end: 0,
            m_window_start: 0,
            buff: VecDeque::with_capacity(wsize - msize + 1),
            m_shift: 2 * (msize - 1) as u64,
            k_shift: 2 * (wsize - 1) as u64,
        }
    }
}

// technique adopted from https://github.com/lh3/minimap2/blob/0cc3cdca27f050fb80a19c90d25ecc6ab0b0907b/sketch.c#L77
impl<'a> Iterator for KmerMinimiserGenerator<'a> {
    type Item = (Kmer, usize, usize, Vec<Kmer>);

    fn next(&mut self) -> Option<Self::Item> {
        let mut min_m_val: u64;
        let mut k_buff = Vec::new();
        let mut prev_k_buff = Vec::new();
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
                // kmer
                self.k_val_f = ((self.k_val_f << 2) | pos_f_val) & self.k_mask;
                self.k_val_r = (self.k_val_r >> 2) | (pos_r_val << self.k_shift);
                self.k_val_l += 1;
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
                if should_return {
                    prev_k_buff.clone_from(&k_buff)
                }
                self.buff_pos = 0;
                self.m_active = u64::MAX;
                self.m_val_f = 0;
                self.m_val_r = 0;
                self.m_val_l = 0;
                self.k_val_f = 0;
                self.k_val_r = 0;
                self.k_val_l = 0;
                self.m_window_end = 0;
                self.m_window_start = self.pos + 1;
                self.buff.clear();
                k_buff.clear();
                self.pos += 1;

                if should_return {
                    return Some((prev_m_val, prev_w_start, prev_w_end, prev_k_buff));
                }
                continue;
            }

            if self.m_val_l < self.msize {
                self.pos += 1;
                continue;
            }

            self.m_val_l -= 1;

            min_m_val = min(self.m_val_f, self.m_val_r);

            if self.k_val_l == self.wsize {
                k_buff.push(min(self.k_val_f, self.k_val_r));
                self.k_val_l -= 1;
            }

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
                        return Some((prev_m_val, prev_w_start, prev_w_end, k_buff));
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
                    return Some((prev_m_val, prev_w_start, prev_w_end, k_buff));
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
                // TODO return strand (implement only when needed)
                return Some((self.m_active, self.m_window_start, self.seq.len(), k_buff));
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
        let mg = KmerMinimiserGenerator::new(b"ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA", 31, 7);
        let expected = [
            ("ATGCGATATCGTAGGCGTCGATGGAGAGCTA", "ACGATAT"),
            ("CTAGCTCTCCATCGACGCCTACGATATCGCA", "ACGATAT"),
            ("GCGATATCGTAGGCGTCGATGGAGAGCTAGA", "ACGATAT"),
            ("ATCTAGCTCTCCATCGACGCCTACGATATCG", "ACGATAT"),
            ("GATATCGTAGGCGTCGATGGAGAGCTAGATC", "ACGATAT"),
            ("ATATCGTAGGCGTCGATGGAGAGCTAGATCG", "ACGATAT"),
            ("TATCGTAGGCGTCGATGGAGAGCTAGATCGA", "ACGATAT"),
            ("ATCGATCTAGCTCTCCATCGACGCCTACGAT", "ACGCCTA"),
            ("GATCGATCTAGCTCTCCATCGACGCCTACGA", "ACGCCTA"),
            ("CGATCGATCTAGCTCTCCATCGACGCCTACG", "ACGCCTA"),
            ("GTAGGCGTCGATGGAGAGCTAGATCGATCGA", "ACGCCTA"),
            ("ATCGATCGATCTAGCTCTCCATCGACGCCTA", "ACGCCTA"),
            ("AGGCGTCGATGGAGAGCTAGATCGATCGATC", "ACGCCTA"),
            ("AGATCGATCGATCTAGCTCTCCATCGACGCC", "AGAGCTA"),
            ("GCGTCGATGGAGAGCTAGATCGATCGATCTA", "AGAGCTA"),
            ("CGTCGATGGAGAGCTAGATCGATCGATCTAA", "AGAGCTA"),
            ("GTCGATGGAGAGCTAGATCGATCGATCTAAA", "AGAGCTA"),
            ("ATTTAGATCGATCGATCTAGCTCTCCATCGA", "AGAGCTA"),
            ("CGATGGAGAGCTAGATCGATCGATCTAAATC", "AGAGCTA"),
            ("GATGGAGAGCTAGATCGATCGATCTAAATCC", "AGAGCTA"),
            ("ATGGAGAGCTAGATCGATCGATCTAAATCCC", "AGAGCTA"),
            ("CGGGATTTAGATCGATCGATCTAGCTCTCCA", "AAATCCC"),
            ("GGAGAGCTAGATCGATCGATCTAAATCCCGA", "AAATCCC"),
            ("ATCGGGATTTAGATCGATCGATCTAGCTCTC", "AAATCCC"),
            ("AGAGCTAGATCGATCGATCTAAATCCCGATC", "AAATCCC"),
            ("CGATCGGGATTTAGATCGATCGATCTAGCTC", "AAATCCC"),
            ("AGCTAGATCGATCGATCTAAATCCCGATCGA", "AAATCCC"),
            ("ATCGATCGGGATTTAGATCGATCGATCTAGC", "AAATCCC"),
            ("AATCGATCGGGATTTAGATCGATCGATCTAG", "AAATCCC"),
            ("GAATCGATCGGGATTTAGATCGATCGATCTA", "AAATCCC"),
            ("AGATCGATCGATCTAAATCCCGATCGATTCC", "AAATCCC"),
            ("CGGAATCGATCGGGATTTAGATCGATCGATC", "AAATCCC"),
            ("ATCGATCGATCTAAATCCCGATCGATTCCGA", "AAATCCC"),
            ("CTCGGAATCGATCGGGATTTAGATCGATCGA", "AAATCCC"),
            ("CGATCGATCTAAATCCCGATCGATTCCGAGC", "AAATCCC"),
            ("CGCTCGGAATCGATCGGGATTTAGATCGATC", "AAATCCC"),
            ("ATCGATCTAAATCCCGATCGATTCCGAGCGC", "AAATCCC"),
            ("CGCGCTCGGAATCGATCGGGATTTAGATCGA", "AAATCCC"),
            ("CGATCTAAATCCCGATCGATTCCGAGCGCGA", "AAATCCC"),
            ("ATCGCGCTCGGAATCGATCGGGATTTAGATC", "AAATCCC"),
            ("ATCTAAATCCCGATCGATTCCGAGCGCGATC", "AAATCCC"),
            ("TCTAAATCCCGATCGATTCCGAGCGCGATCA", "AAATCCC"),
            ("CTAAATCCCGATCGATTCCGAGCGCGATCAA", "AAATCCC"),
            ("TAAATCCCGATCGATTCCGAGCGCGATCAAA", "AAATCCC"),
            ("AAATCCCGATCGATTCCGAGCGCGATCAAAG", "AAATCCC"),
            ("AATCCCGATCGATTCCGAGCGCGATCAAAGC", "AAATCCC"),
            ("ATCCCGATCGATTCCGAGCGCGATCAAAGCG", "AATCCCG"),
            ("GCGCTTTGATCGCGCTCGGAATCGATCGGGA", "AATCGAT"),
            ("CCCGATCGATTCCGAGCGCGATCAAAGCGCG", "AAAGCGC"),
            ("CCGATCGATTCCGAGCGCGATCAAAGCGCGA", "AAAGCGC"),
            ("ATCGCGCTTTGATCGCGCTCGGAATCGATCG", "AAAGCGC"),
            ("GATCGATTCCGAGCGCGATCAAAGCGCGATA", "AAAGCGC"),
            ("ATCGATTCCGAGCGCGATCAAAGCGCGATAG", "AAAGCGC"),
            ("CCTATCGCGCTTTGATCGCGCTCGGAATCGA", "AAAGCGC"),
            ("CGATTCCGAGCGCGATCAAAGCGCGATAGGC", "AAAGCGC"),
            ("AGCCTATCGCGCTTTGATCGCGCTCGGAATC", "AAAGCGC"),
            ("ATTCCGAGCGCGATCAAAGCGCGATAGGCTA", "AAAGCGC"),
            ("CTAGCCTATCGCGCTTTGATCGCGCTCGGAA", "AAAGCGC"),
            ("GCTAGCCTATCGCGCTTTGATCGCGCTCGGA", "AAAGCGC"),
            ("AGCTAGCCTATCGCGCTTTGATCGCGCTCGG", "AAAGCGC"),
            ("CGAGCGCGATCAAAGCGCGATAGGCTAGCTA", "AAAGCGC"),
            ("GAGCGCGATCAAAGCGCGATAGGCTAGCTAA", "AAAGCGC"),
            ("AGCGCGATCAAAGCGCGATAGGCTAGCTAAA", "AAAGCGC"),
            ("CTTTAGCTAGCCTATCGCGCTTTGATCGCGC", "AAAGCGC"),
            ("CGCGATCAAAGCGCGATAGGCTAGCTAAAGC", "AAAGCGC"),
            ("AGCTTTAGCTAGCCTATCGCGCTTTGATCGC", "AAAGCGC"),
            ("CGATCAAAGCGCGATAGGCTAGCTAAAGCTA", "AAAGCGC"),
            ("CTAGCTTTAGCTAGCCTATCGCGCTTTGATC", "AAAGCGC"),
            ("ATCAAAGCGCGATAGGCTAGCTAAAGCTAGC", "AAAGCGC"),
            ("TCAAAGCGCGATAGGCTAGCTAAAGCTAGCA", "AAAGCGC"),
        ];
        let mut i = 0;
        for (m, _, _, ks) in mg {
            for k in ks {
                println!("{}, {}", numeric_to_kmer(k, 31), numeric_to_kmer(m, 7));
                assert_eq!(numeric_to_kmer(k, 31), expected[i].0);
                assert_eq!(numeric_to_kmer(m, 7), expected[i].1);
                i += 1;
            }
        }
    }

    #[test]
    fn minimisers_generated_with_error_test() {
        // Acquired from https://homolog.us/blogs/bioinfo/2017/10/25/intro-minimizer/
        let mg = KmerMinimiserGenerator::new(b"ATGCGATATCGNTAGGCGTCGATGGA", 8, 5);
        let expected = [
            ("ATGCGATA", "ATCGC"),
            ("ATATCGCA", "ATCGC"),
            ("GATATCGC", "ATATC"),
            ("CGATATCG", "ATATC"),
            ("GACGCCTA", "ACGCC"),
            ("AGGCGTCG", "ACGCC"),
            ("GGCGTCGA", "ACGCC"),
            ("ATCGACGC", "ACGCC"),
            ("CATCGACG", "ATCGA"),
            ("CCATCGAC", "ATCGA"),
            ("TCCATCGA", "ATCGA"),
        ];
        let mut i = 0;
        for (m, _, _, ks) in mg {
            for k in ks {
                println!(
                    "Kmer: {}, Minimiser: {}",
                    numeric_to_kmer(k, 8),
                    numeric_to_kmer(m, 5)
                );
                assert_eq!(numeric_to_kmer(k, 8), expected[i].0);
                assert_eq!(numeric_to_kmer(m, 5), expected[i].1);
                i += 1;
            }
        }
    }
}
