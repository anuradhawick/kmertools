use super::Kmer;
use std::collections::{HashMap, HashSet};
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

pub struct KmerGenerator<'a> {
    seq: &'a [u8],
    fval: u64,
    rval: u64,
    len: usize,
    pos: usize,
    ksize: usize,
    mask: u64,
    shift: u64,
}

impl<'a> KmerGenerator<'a> {
    pub fn new(seq: &'a [u8], ksize: usize) -> Self {
        KmerGenerator {
            seq,
            fval: 0,
            rval: 0,
            len: 0,
            pos: 0,
            ksize,
            mask: (1_u64 << (2 * ksize)) - 1,
            shift: 2 * (ksize - 1) as u64,
        }
    }

    pub fn rev_comp(kmer: Kmer, ksize: usize) -> Kmer {
        let mut rkmer = 0;
        let mut kmer = kmer;
        for _ in 0..ksize as u64 {
            rkmer <<= 2;
            rkmer |= (kmer & 3) ^ REV_MASK;
            kmer >>= 2;
        }
        rkmer
    }

    pub fn kmer_pos_maps(ksize: usize) -> (Vec<usize>, HashMap<usize, u64>, usize) {
        // this function returns a vector of size 4 ^ k that maps to (4 ^ k) / 2 or (4 ^ k) / 2 + 4 ^ (k / 2) / 2
        // behaves as a non-min but perfect hash (MPHF) function that maps kmers to an index of vector we want
        let mut min_mer_set = HashSet::new();
        let mut pos_min_mer_map = HashMap::new();
        let mut min_mer_pos_map = vec![0_usize; 4_u64.pow(ksize as u32) as usize];
        for kmer in 0..(4_u64.pow(ksize as u32)) {
            let min_mer = u64::min(kmer, KmerGenerator::rev_comp(kmer, ksize));
            min_mer_set.insert(min_mer);
        }
        let count = min_mer_set.len();
        let mut min_mer_vec: Vec<u64> = Vec::from_iter(min_mer_set);
        min_mer_vec.sort();

        for (pos, &kmer) in min_mer_vec.iter().enumerate() {
            min_mer_pos_map[kmer as usize] = pos;
            pos_min_mer_map.insert(pos, kmer);
        }
        (min_mer_pos_map, pos_min_mer_map, count)
    }
}

// technique adopted from https://github.com/lh3/minimap2/blob/0cc3cdca27f050fb80a19c90d25ecc6ab0b0907b/sketch.c#L77
impl<'a> Iterator for KmerGenerator<'a> {
    type Item = (Kmer, Kmer);

    fn next(&mut self) -> Option<(Kmer, Kmer)> {
        // valid base
        loop {
            if self.pos == self.seq.len() {
                return None;
            }
            let pos_char = self.seq[self.pos];
            let pos_f_val = SEQ_NT4_TABLE[pos_char as usize] as u64;
            let pos_r_val = pos_f_val ^ REV_MASK;
            self.pos += 1;

            if pos_f_val < 4 {
                // non ambiguous
                self.fval = ((self.fval << 2) | pos_f_val) & self.mask;
                self.rval = (self.rval >> 2) | (pos_r_val << self.shift);
                self.len += 1;
            } else {
                // ambiguous
                self.len = 0;
            }

            if self.len == self.ksize {
                self.len -= 1;
                return Some((self.fval, self.rval));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmers_generated_test() {
        let mut kg = KmerGenerator::new(b"ACGT", 2);
        let kmer1 = kg.next();
        let kmer2 = kg.next();
        let kmer3 = kg.next();
        let kmer4 = kg.next();
        // AC 00 01, GT 10 11
        assert_eq!(kmer1, Some((1, 11)));
        // CG 01 10, CG 01 10
        assert_eq!(kmer2, Some((6, 6)));
        // GT 10 11, AC 00 01
        assert_eq!(kmer3, Some((11, 1)));
        // None
        assert_eq!(kmer4, None);
    }

    #[test]
    fn kmers_generated_ambiguous_test() {
        let mut kg = KmerGenerator::new(b"ACNGTT", 2);
        let kmer1 = kg.next();
        let kmer2 = kg.next();
        let kmer3 = kg.next();
        let kmer4 = kg.next();
        // AC 00 01, GT 10 11
        assert_eq!(kmer1, Some((1, 11)));
        // GT 10 11, AC 00 01
        assert_eq!(kmer2, Some((11, 1)));
        // TT 11 11, AA 00 00
        assert_eq!(kmer3, Some((15, 0)));
        // None
        assert_eq!(kmer4, None);
    }

    #[test]
    fn rev_comp_test() {
        // ACGT 00 01 10 11 -> ACGT 00 01 10 11
        assert_eq!(KmerGenerator::rev_comp(0b00011011, 4), 0b00011011);
        // ATCGGT 00 11 01 10 10 11 -> ACCGAT 00 01 01 10 00 11
        assert_eq!(KmerGenerator::rev_comp(0b001101101011, 6), 0b000101100011);
    }

    #[test]
    fn pos_map_test() {
        let (min_mer_pos_map, pos_min_mer_map, min_mer_count) = KmerGenerator::kmer_pos_maps(4);
        let mut pos_index_count = 0;

        assert_eq!(min_mer_count, 136);
        assert_eq!(pos_min_mer_map.len(), 136);

        for pos in min_mer_pos_map.clone() {
            if pos > 0 {
                pos_index_count += 1;
            }
            assert!(pos < 136);
        }
        assert!(pos_index_count == 135);
        // AAAA -> 0
        assert_eq!(min_mer_pos_map[0], 0);
        // TTTT -> 0
        assert_eq!(min_mer_pos_map[0b11111111], 0);
        // AAAT -> 11
        assert_eq!(min_mer_pos_map[0b11], 0b11);
    }
}
