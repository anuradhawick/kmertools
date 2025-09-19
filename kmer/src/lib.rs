pub mod kmer;
pub mod kmer_minimisers;
pub mod minimiser;
pub type Kmer = u64;

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

pub fn numeric_to_kmer(kmer: u64, k: usize) -> String {
    let mut s = String::new();
    let mut kmer = kmer;
    for _ in 0..k {
        let c = match kmer & 0b11 {
            0b00 => 'A',
            0b01 => 'C',
            0b10 => 'G',
            0b11 => 'T',
            _ => panic!("Impossible!"),
        };
        s.push(c);
        kmer >>= 2;
    }
    s.chars().rev().collect()
}

pub fn kmer_to_numeric(kmer: &str) -> (u64, u64) {
    let mut fval = 0;
    let mut rval = 0;
    let shift = 2 * (kmer.len() - 1) as u64;
    let mask = (1_u64 << (2 * kmer.len())) - 1;

    for c in kmer.chars() {
        let pos_f_val = SEQ_NT4_TABLE[c as usize] as u64;
        let pos_r_val = pos_f_val ^ REV_MASK;
        fval = ((fval << 2) | pos_f_val) & mask;
        rval = (rval >> 2) | (pos_r_val << shift);
    }

    (fval, rval)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_to_kmer_test() {
        let kmer_1 = numeric_to_kmer(0b0001101111, 5);
        let kmer_2 = numeric_to_kmer(0b0000011011, 5);

        assert_eq!(kmer_1, "ACGTT");
        assert_eq!(kmer_2, "AACGT");
    }

    #[test]
    fn kmer_to_numeric_test() {
        let (fval, rval) = kmer_to_numeric("ACGTT"); // rev AACGT

        assert_eq!(fval, 0b0001101111);
        assert_eq!(rval, 0b0000011011);
    }
}
