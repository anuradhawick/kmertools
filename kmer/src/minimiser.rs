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
    ksize: usize,
    msize: usize,
    k_mask: u64,
    k_val_f: u64,
    k_val_r: u64,
    k_val_l: usize,
    m_mask: u64,
    m_window_start: usize,
    m_window_end: usize,
    m_val_f: u64,
    m_val_r: u64,
    m_val_l: usize,
    m_active: u64,
    m_shift: u64,
    k_shift: u64,
    buff: VecDeque<u64>,
    buff_pos: usize,
}

impl<'a> MinimiserGenerator<'a> {
    pub fn new(seq: &'a [u8], ksize: usize, msize: usize) -> Self {
        MinimiserGenerator {
            seq,
            ksize,
            msize,
            pos: 0,
            buff_pos: 0,
            k_mask: (1_u64 << (2 * ksize)) - 1,
            k_val_f: 0,
            k_val_r: 0,
            k_val_l: 0,
            m_active: u64::MAX,
            m_mask: (1_u64 << (2 * msize)) - 1,
            m_val_f: 0,
            m_val_r: 0,
            m_val_l: 0,
            m_window_end: 0,
            m_window_start: 0,
            buff: VecDeque::with_capacity(ksize - msize + 1),
            k_shift: 2 * (ksize - 1) as u64,
            m_shift: 2 * (msize - 1) as u64,
        }
    }
}

// technique adopted from https://github.com/lh3/minimap2/blob/0cc3cdca27f050fb80a19c90d25ecc6ab0b0907b/sketch.c#L77
impl<'a> Iterator for MinimiserGenerator<'a> {
    type Item = (Kmer, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let mut min_m_val = 0;
        // valid base
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
                self.buff_pos = 0;
                self.k_mask = (1_u64 << (2 * self.ksize)) - 1;
                self.k_val_f = 0;
                self.k_val_r = 0;
                self.k_val_l = 0;
                self.m_active = u64::MAX;
                self.m_mask = (1_u64 << (2 * self.msize)) - 1;
                self.m_val_f = 0;
                self.m_val_r = 0;
                self.m_val_l = 0;
                self.m_window_end = 0;
                self.m_window_start = 0;
                self.buff.clear();
                self.pos += 1;
                continue;
            }

            if self.m_val_l < self.msize {
                self.pos += 1;
                continue;
            }

            self.m_val_l -= 1;
            self.k_val_l -= 1;
            min_m_val = min(self.m_val_f, self.m_val_r);
            // min_m_val = self.m_val_f;

            // minimiser buffer is full
            if self.buff.len() == self.ksize - self.msize + 1 {
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
                        println!(
                            "Window start:{:5} end:{:5} - {:60} : {:10}",
                            self.m_window_start,
                            self.m_window_end,
                            String::from_utf8(
                                self.seq[self.m_window_start..self.m_window_end].to_vec()
                            )
                            .unwrap(),
                            super::numeric_to_kmer(self.m_active, self.msize)
                        );
                        self.m_active = new_min;
                        self.m_window_start = self.pos - self.ksize + 1;
                        self.pos += 1;
                        return Some((min_m_val, self.m_window_start, self.m_window_end));
                    }
                } else if min_m_val < self.m_active {
                    // break the window
                    self.m_window_end = self.pos;
                    println!(
                        "Window start:{:5} end:{:5} - {:60} : {:10}",
                        self.m_window_start,
                        self.m_window_end,
                        String::from_utf8(
                            self.seq[self.m_window_start..self.m_window_end].to_vec()
                        )
                        .unwrap(),
                        super::numeric_to_kmer(self.m_active, self.msize)
                    );
                    self.m_active = min_m_val;
                    self.buff_pos = self.buff.len() - 1;
                    self.m_window_start = self.pos - self.ksize + 1;
                    self.pos += 1;
                    return Some((min_m_val, self.m_window_start, self.m_window_end));
                } else {
                    self.buff_pos -= 1;
                }
            } else {
                // add new minimizer to buffer
                self.buff.push_back(min_m_val);
            }

            // first time we are experiencing all minimizers
            if self.m_active == u64::MAX && self.buff.len() == self.ksize - self.msize + 1 {
                for j in 0..self.buff.len() {
                    if *self.buff.get(j).unwrap() < self.m_active {
                        self.buff_pos = j;
                        self.m_active = *self.buff.get(j).unwrap();
                    }
                }
            }

            if self.pos == self.seq.len() - 1 {
                println!(
                    "Window start:{:5} end:{:5} - {:60} : {:10}",
                    self.m_window_start,
                    self.seq.len(),
                    String::from_utf8(self.seq[self.m_window_start..self.seq.len()].to_vec())
                        .unwrap(),
                    super::numeric_to_kmer(self.m_active, self.msize)
                );
                self.pos += 1;
                return Some((min_m_val, self.m_window_start, self.m_window_end));
            }

            self.pos += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimisers_generated_test() {
        let mut mg = MinimiserGenerator::new(b"ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA", 31, 7);
        // let kmer1 = mg.next();
        // let kmer2 = mg.next();
        // let kmer3 = mg.next();
        // let kmer4 = mg.next();
        // let kmer4 = mg.next();
        // let kmer4 = mg.next();
        for (i, (k, s, e)) in mg.enumerate() {}
        // // AC 00 01, GT 10 11
        // assert_eq!(kmer1, Some((1, 11)));
        // // CG 01 10, CG 01 10
        // assert_eq!(kmer2, Some((6, 6)));
        // // GT 10 11, AC 00 01
        // assert_eq!(kmer3, Some((11, 1)));
        // // None
        // assert_eq!(kmer4, None);
    }
}
