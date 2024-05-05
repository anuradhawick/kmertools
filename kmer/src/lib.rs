pub mod kmer;
pub mod kmer_minimisers;
pub mod minimiser;
pub type Kmer = u64;

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
