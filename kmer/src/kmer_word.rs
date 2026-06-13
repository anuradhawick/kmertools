use ruint::Uint;
use std::{
    fmt::Debug,
    hash::Hash,
    ops::{BitAnd, BitOr, BitXor, Rem, Shl, Shr},
};

/// Integer operations required by a packed, two-bit k-mer.
pub trait KmerWord:
    Copy
    + Debug
    + Eq
    + Hash
    + Ord
    + Send
    + Sync
    + Shl<usize, Output = Self>
    + Shr<usize, Output = Self>
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + Rem<Output = Self>
    + 'static
{
    const BITS: usize;
    const ZERO: Self;
    const MAX: Self;

    fn from_u64(value: u64) -> Self;
    fn to_u64(self) -> u64;
}

// Native integers use Rust's built-in constants and lossless `as` conversion
// from u64. Keeping this separate preserves the direct u64/u128 fast path.
macro_rules! impl_primitive_kmer_word {
    ($word:ty) => {
        impl KmerWord for $word {
            const BITS: usize = <$word>::BITS as usize;
            const ZERO: Self = 0;
            const MAX: Self = <$word>::MAX;

            #[inline]
            fn from_u64(value: u64) -> Self {
                value as Self
            }

            #[inline]
            fn to_u64(self) -> u64 {
                self as u64
            }
        }
    };
}

impl_primitive_kmer_word!(u64);
impl_primitive_kmer_word!(u128);

pub type KmerU128 = Uint<128, 2>;
pub type KmerU256 = Uint<256, 4>;
pub type KmerU512 = Uint<512, 8>;
pub type KmerU1024 = Uint<1024, 16>;
pub type KmerU2048 = Uint<2048, 32>;

// `ruint::Uint` does not implement `From<u64>`; conversion goes through its
// own API. These implementations cover only the preset widths we dispatch to.
macro_rules! impl_uint_kmer_word {
    ($word:ty, $bits:expr) => {
        impl KmerWord for $word {
            const BITS: usize = $bits;
            const ZERO: Self = Self::ZERO;
            const MAX: Self = Self::MAX;

            #[inline]
            fn from_u64(value: u64) -> Self {
                Self::wrapping_from(value)
            }

            #[inline]
            fn to_u64(self) -> u64 {
                self.wrapping_to()
            }
        }
    };
}

impl_uint_kmer_word!(KmerU128, 128);
impl_uint_kmer_word!(KmerU256, 256);
impl_uint_kmer_word!(KmerU512, 512);
impl_uint_kmer_word!(KmerU1024, 1024);
impl_uint_kmer_word!(KmerU2048, 2048);

#[cfg(test)]
mod tests {
    use super::*;

    fn encode<K: KmerWord>(bases: &[u8]) -> K {
        let used_bits = bases.len() * 2;
        assert!(used_bits <= K::BITS);
        let mask = K::MAX >> (K::BITS - used_bits);

        bases.iter().fold(K::ZERO, |kmer, &base| {
            ((kmer << 2) | K::from_u64(base as u64)) & mask
        })
    }

    #[test]
    fn raw_integer_operations_match_u64() {
        let bases = [0, 1, 2, 3, 3, 2, 1, 0];
        let native = encode::<u64>(&bases);
        let wide = encode::<KmerU256>(&bases);

        assert_eq!(wide.as_limbs()[0], native);
        assert!(wide.as_limbs()[1..].iter().all(|&limb| limb == 0));
    }

    #[test]
    fn raw_reverse_shift_and_or_match_u64() {
        #[allow(clippy::identity_op)]
        // Here I am validating actual operations, nothing to lint away.
        let native = (0b0001_u64 >> 2) | (0b11_u64 << 6);
        let wide: KmerU256 =
            (KmerU256::from_u64(0b0001) >> 2_usize) | (KmerU256::from_u64(0b11) << 6_usize);

        assert_eq!(wide.as_limbs()[0], native);
    }
}
