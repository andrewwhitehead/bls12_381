//! This module provides an implementation of the BLS12-381 scalar field $\mathbb{F}_q$
//! where `q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`

use core::fmt;
use core::ops;

use crypto_bigint::{Encoding, Zero, U256};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::util::{uint_is_nonzero_vartime, Montgomery};

pub(crate) type Inner = U256;

/// Represents an element of the scalar field $\mathbb{F}_q$ of the BLS12-381 elliptic
/// curve construction.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Scalar` values are always in
// Montgomery form; i.e., Scalar(a) = aR mod q, with R = 2^256.
#[derive(Clone, Copy, Eq)]
#[repr(transparent)]
pub struct Scalar(pub(crate) Inner);

impl fmt::Debug for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl From<u64> for Scalar {
    #[inline]
    fn from(val: u64) -> Scalar {
        Scalar::from_canonical(Inner::from_u64(val))
    }
}

impl From<Inner> for Scalar {
    #[inline]
    fn from(value: Inner) -> Scalar {
        Scalar::from_canonical(value)
    }
}

impl From<&Inner> for Scalar {
    #[inline]
    fn from(value: &Inner) -> Scalar {
        Scalar::from_canonical(*value)
    }
}

impl ConstantTimeEq for Scalar {
    #[inline]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl PartialEq for Scalar {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Scalar {
    #[inline]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Scalar(Inner::conditional_select(&a.0, &b.0, choice))
    }
}

/// The field definition with associated helper functions.
pub(crate) const FIELD: Montgomery<{ Inner::LIMBS }> = Montgomery::new(Inner::from_be_hex(
    "73eda753299d7d48\
     3339d80809a1d805\
     53bda402fffe5bfe\
     ffffffff00000001",
));

/// GENERATOR = 7 (multiplicative generator of r-1 order, that is also quadratic nonresidue)
const GENERATOR: Inner = Inner::from_be_hex(
    "351332208fc5a8c4\
     ff9c57876f8457b0\
     17e363d300189c0f\
     0000000efffffff1",
);

impl<'a> ops::Neg for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        self.neg()
    }
}

impl ops::Neg for Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        -&self
    }
}

impl<'a, 'b> ops::Sub<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn sub(self, rhs: &'b Scalar) -> Scalar {
        self.sub(rhs)
    }
}

impl<'a, 'b> ops::Add<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn add(self, rhs: &'b Scalar) -> Scalar {
        self.add(rhs)
    }
}

impl<'a, 'b> ops::Mul<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn mul(self, rhs: &'b Scalar) -> Scalar {
        self.mul(rhs)
    }
}

impl_binops_additive!(Scalar, Scalar);
impl_binops_multiplicative!(Scalar, Scalar);

// 2^S * t = MODULUS - 1 with t odd
const S: u32 = 32;

/// GENERATOR^t where t * 2^s + 1 = q
/// with t odd. In other words, this
/// is a 2^s root of unity.
///
/// `GENERATOR = 7 mod q` is a generator
/// of the q - 1 order multiplicative
/// subgroup.
const ROOT_OF_UNITY: Inner = Inner::from_be_hex(
    "5bf3adda19e9b27b\
     0af53ae352a31e64\
     5b1b4c801819d7ec\
     b9b58d8c5f0e466a",
);

impl Default for Scalar {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(feature = "zeroize")]
impl zeroize::DefaultIsZeroes for Scalar {}

impl Scalar {
    /// The maximum size of a scalar value in bits.
    pub const BIT_SIZE: usize = 255;

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Scalar {
        Scalar(Inner::ZERO)
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Scalar {
        Scalar(FIELD.one())
    }

    /// Returns true iff this element is zero.
    #[inline]
    pub fn is_zero(&self) -> Choice {
        self.0.is_zero()
    }

    /// Returns true iff this element is zero.
    #[inline]
    pub fn is_zero_vartime(&self) -> bool {
        !uint_is_nonzero_vartime(&self.0)
    }

    /// Doubles this field element.
    #[inline]
    pub const fn double(&self) -> Scalar {
        Scalar(FIELD.double(&self.0))
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 32]) -> CtOption<Scalar> {
        FIELD
            .try_from_canonical(&Inner::from_le_bytes(*bytes))
            .map(|s| Scalar(s))
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes_vartime(bytes: &[u8; 32]) -> Option<Scalar> {
        FIELD
            .try_from_canonical_vartime(&Inner::from_le_bytes(*bytes))
            .map(|s| Scalar(s))
    }

    /// Convert a big-endian hex representation of a scalar into a `Scalar`
    /// without checking that it is canonical.
    pub const fn from_hex_unchecked(s: &str) -> Scalar {
        Scalar(FIELD.from_canonical(&Inner::from_be_hex(s)))
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 32] {
        self.to_canonical().to_le_bytes()
    }

    /// Converts a 512-bit little endian integer into a `Scalar`
    /// by reducing by the modulus.
    pub fn from_bytes_wide(bytes: &[u8; 64]) -> Scalar {
        let lo = Inner::from_le_bytes(bytes[0..32].try_into().unwrap());
        let hi = Inner::from_le_bytes(bytes[32..64].try_into().unwrap());
        Scalar(FIELD.from_canonical_wide(&lo, &hi))
    }

    /// Converts from a canonical scalar represented by a U256.
    #[inline]
    pub(crate) const fn from_canonical(val: Inner) -> Self {
        Scalar(FIELD.from_canonical(&val))
    }

    /// Turn into canonical form by computing
    /// (a.R) / R = a
    #[inline]
    pub const fn to_canonical(&self) -> Inner {
        FIELD.to_canonical(&self.0)
    }

    /// Converts from an integer represented in little endian
    /// into its (congruent) `Scalar` representation.
    pub const fn from_raw(val: [u64; 4]) -> Self {
        Scalar::from_canonical(inner_from_raw(val))
    }

    #[inline]
    #[allow(unused)]
    pub(crate) const fn to_raw(&self) -> [u64; 4] {
        inner_to_raw(self.0)
    }

    #[inline]
    #[allow(unused)]
    pub(crate) const fn from_raw_unchecked(val: [u64; 4]) -> Self {
        Scalar(inner_from_raw(val))
    }

    /// Squares this element.
    #[inline]
    pub const fn square(&self) -> Scalar {
        Scalar(FIELD.square(&self.0))
    }

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // https://eprint.iacr.org/2012/685.pdf (page 12, algorithm 5)

        // exp = (t - 1) // 2 = ((q - 1) >> S) - 1) // 2
        //     = 6104339283789297388802252303364915521546564123189034618274734669823
        const EXP: Inner = Inner::from_be_hex(
            "0000000039f6d3a9\
             94cebea4199cec04\
             04d0ec02a9ded201\
             7fff2dff7fffffff",
        );

        let w = FIELD.pow_vartime(&self.0, &EXP);
        let mut v = S;
        let mut x = FIELD.mul_inline(&self.0, &w);
        let mut b = FIELD.mul_inline(&x, &w);

        // Initialize z as the 2^S root of unity.
        let mut z = ROOT_OF_UNITY;

        for max_v in (1..=S).rev() {
            let mut k = 1;
            let mut tmp = FIELD.square_inline(&b);
            let mut j_less_than_v: Choice = 1.into();

            for j in 2..max_v {
                let tmp_is_one = tmp.ct_eq(&FIELD.one());
                let sel = Inner::conditional_select(&tmp, &z, tmp_is_one);
                let squared = FIELD.square_inline(&sel);
                tmp = Inner::conditional_select(&squared, &tmp, tmp_is_one);
                let new_z = Inner::conditional_select(&z, &squared, tmp_is_one);
                j_less_than_v &= !j.ct_eq(&v);
                k = u32::conditional_select(&j, &k, tmp_is_one);
                z = Inner::conditional_select(&z, &new_z, j_less_than_v);
            }

            let result = FIELD.mul_inline(&x, &z);
            x = Inner::conditional_select(&result, &x, b.ct_eq(&FIELD.one()));
            z = FIELD.square_inline(&z);
            b = FIELD.mul_inline(&b, &z);
            v = k;
        }

        let xsq = FIELD.square_inline(&x);
        CtOption::new(
            Scalar(x),
            xsq.ct_eq(&self.0), // Only return Some if it's the square root.
        )
    }

    /// Exponentiates `self` by `by`, where `by` is a
    /// little-endian order integer exponent.
    pub fn pow(&self, by: &[u64; 4]) -> Self {
        let by = inner_from_raw(*by);
        Scalar(FIELD.pow(&self.0, &by))
    }

    /// Exponentiates `self` by `by`, where `by` is a
    /// little-endian order integer exponent.
    ///
    /// **This operation is variable time with respect
    /// to the exponent.** If the exponent is fixed,
    /// this operation is effectively constant time.
    pub const fn pow_vartime(&self, by: &[u64; 4]) -> Self {
        let by = inner_from_raw(*by);
        Scalar(FIELD.pow_vartime(&self.0, &by))
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn invert(&self) -> CtOption<Self> {
        #[inline(always)]
        fn mul(n: &Inner, other: &Inner) -> Inner {
            FIELD.mul_inline(&n, &other)
        }
        #[inline(always)]
        fn pow2k(n: &Inner, num_times: usize) -> Inner {
            let mut res = *n;
            for _ in 0..num_times {
                res = FIELD.square_inline(&res);
            }
            res
        }
        // found using https://github.com/kwantam/addchain
        let mut t0 = pow2k(&self.0, 1);
        let mut t1 = mul(&t0, &self.0);
        let mut t16 = pow2k(&t0, 1);
        let mut t6 = pow2k(&t16, 1);
        let mut t5 = mul(&t6, &t0);
        t0 = mul(&t6, &t16);
        let mut t12 = mul(&t5, &t16);
        let mut t2 = pow2k(&t6, 1);
        let mut t7 = mul(&t5, &t6);
        let mut t15 = mul(&t0, &t5);
        let mut t17 = pow2k(&t12, 1);
        t1 = mul(&t1, &t17);
        let mut t3 = mul(&t7, &t2);
        let t8 = mul(&t1, &t17);
        let t4 = mul(&t8, &t2);
        let t9 = mul(&t8, &t7);
        t7 = mul(&t4, &t5);
        let t11 = mul(&t4, &t17);
        t5 = mul(&t9, &t17);
        let t14 = mul(&t7, &t15);
        let t13 = mul(&t11, &t12);
        t12 = mul(&t11, &t17);
        t15 = mul(&t15, &t12);
        t16 = mul(&t16, &t15);
        t3 = mul(&t3, &t16);
        t17 = mul(&t17, &t3);
        t0 = mul(&t0, &t17);
        t6 = mul(&t6, &t0);
        t2 = mul(&t2, &t6);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t17);
        t0 = pow2k(&t0, 9);
        t0 = mul(&t0, &t16);
        t0 = pow2k(&t0, 9);
        t0 = mul(&t0, &t15);
        t0 = pow2k(&t0, 9);
        t0 = mul(&t0, &t15);
        t0 = pow2k(&t0, 7);
        t0 = mul(&t0, &t14);
        t0 = pow2k(&t0, 7);
        t0 = mul(&t0, &t13);
        t0 = pow2k(&t0, 10);
        t0 = mul(&t0, &t12);
        t0 = pow2k(&t0, 9);
        t0 = mul(&t0, &t11);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t8);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &self.0);
        t0 = pow2k(&t0, 14);
        t0 = mul(&t0, &t9);
        t0 = pow2k(&t0, 10);
        t0 = mul(&t0, &t8);
        t0 = pow2k(&t0, 15);
        t0 = mul(&t0, &t7);
        t0 = pow2k(&t0, 10);
        t0 = mul(&t0, &t6);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t5);
        t0 = pow2k(&t0, 16);
        t0 = mul(&t0, &t3);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t2);
        t0 = pow2k(&t0, 7);
        t0 = mul(&t0, &t4);
        t0 = pow2k(&t0, 9);
        t0 = mul(&t0, &t2);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t3);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t2);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t2);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t2);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t3);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t2);
        t0 = pow2k(&t0, 8);
        t0 = mul(&t0, &t2);
        t0 = pow2k(&t0, 5);
        t0 = mul(&t0, &t1);
        t0 = pow2k(&t0, 5);
        t0 = mul(&t0, &t1);

        CtOption::new(Scalar(t0), !self.0.is_zero())
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub const fn invert_vartime(&self) -> Option<Self> {
        match FIELD.invert_vartime(&self.0) {
            Some(s) => Some(Scalar(s)),
            None => None,
        }
    }

    /// Efficiently invert an array of scalars.
    pub fn batch_invert(batch: &[Self], output: &mut [Self]) {
        use crypto_bigint::{Limb, Word};
        // use high bit to store skip flag
        const SKIP_WORD: Word = 1 << (Limb::BIT_SIZE - 1);
        const SKIP: Inner = {
            let mut s = [Limb::ZERO; Inner::LIMBS];
            s[Inner::LIMBS - 1].0 = SKIP_WORD;
            Inner::new(s)
        };
        let mut acc = Scalar::one();
        for (a, b) in batch.iter().zip(output.iter_mut()) {
            let skip = a.is_zero();
            b.0 = Inner::conditional_select(&acc.0, &SKIP, skip);
            acc *= Scalar::conditional_select(&a, &Scalar::one(), skip);
        }
        acc = acc.invert().unwrap();
        for (a, b) in batch.iter().zip(output.iter_mut()).rev() {
            let hi_word = b.0.limbs()[Inner::LIMBS - 1].0;
            let skip = Choice::from((hi_word >> (Word::BITS - 1)) as u8);
            b.0.limbs_mut()[Inner::LIMBS - 1].0 = hi_word & !SKIP_WORD;
            *b *= acc;
            acc *= Scalar::conditional_select(&a, &Scalar::one(), skip);
        }
    }

    /// Efficiently invert an array of scalars in-place.
    #[cfg(feature = "alloc")]
    #[inline]
    pub fn batch_invert_in_place(batch: &mut [Self]) {
        let inputs = alloc::vec::Vec::from_iter(batch.iter().copied());
        Self::batch_invert(&inputs, batch)
    }

    /// Multiplies `rhs` by `self`, returning the result.
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        Scalar(FIELD.mul(&self.0, &rhs.0))
    }

    /// Subtracts `rhs` from `self`, returning the result.
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        Scalar(FIELD.sub(&self.0, &rhs.0))
    }

    /// Adds `rhs` to `self`, returning the result.
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        Scalar(FIELD.add(&self.0, &rhs.0))
    }

    /// Negates `self`.
    #[inline]
    pub const fn neg(&self) -> Self {
        Scalar(FIELD.neg(&self.0))
    }

    /// Get the bit value at index `index`.
    #[inline]
    pub const fn bit(&self, index: usize) -> u8 {
        self.0.bit_vartime(index) as u8
    }

    /// Generate a random scalar.
    pub fn random(mut rng: impl RngCore) -> Self {
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
    }
}

impl From<Scalar> for [u8; 32] {
    fn from(value: Scalar) -> [u8; 32] {
        value.to_bytes()
    }
}

impl<'a> From<&Scalar> for [u8; 32] {
    fn from(value: &Scalar) -> [u8; 32] {
        value.to_bytes()
    }
}

impl ff::Field for Scalar {
    fn random(rng: impl RngCore) -> Self {
        Self::random(rng)
    }

    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::one()
    }

    fn is_zero(&self) -> Choice {
        self.is_zero()
    }

    fn is_zero_vartime(&self) -> bool {
        self.is_zero_vartime()
    }

    #[must_use]
    fn square(&self) -> Self {
        self.square()
    }

    #[must_use]
    fn double(&self) -> Self {
        self.double()
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }

    fn sqrt(&self) -> CtOption<Self> {
        self.sqrt()
    }
}

impl ff::PrimeField for Scalar {
    type Repr = [u8; 32];

    fn from_repr(r: Self::Repr) -> CtOption<Self> {
        Self::from_bytes(&r)
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_bytes()
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_bytes()[0] & 1)
    }

    const NUM_BITS: u32 = FIELD.modulus.bits_vartime() as u32;
    const CAPACITY: u32 = Self::NUM_BITS - 1;

    fn multiplicative_generator() -> Self {
        Scalar(GENERATOR)
    }

    const S: u32 = S;

    fn root_of_unity() -> Self {
        Scalar(ROOT_OF_UNITY)
    }
}

#[cfg(all(feature = "bits", target_pointer_width = "32"))]
type ReprBits = [u32; 8];

#[cfg(all(feature = "bits", target_pointer_width = "64"))]
type ReprBits = [u64; 4];

#[cfg(feature = "bits")]
impl ff::PrimeFieldBits for Scalar {
    type ReprBits = ReprBits;

    fn to_le_bits(&self) -> ff::FieldBits<Self::ReprBits> {
        ff::FieldBits::new(self.to_canonical().to_words())
    }

    fn char_le_bits() -> ff::FieldBits<Self::ReprBits> {
        ff::FieldBits::new(FIELD.modulus.to_words())
    }
}

impl<T> core::iter::Sum<T> for Scalar
where
    T: core::borrow::Borrow<Scalar>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
    }
}

#[inline(always)]
#[cfg(target_pointer_width = "32")]
const fn inner_from_raw(arr: [u64; 4]) -> Inner {
    const MASK: u64 = u32::MAX as u64;
    Inner::from_words([
        (arr[0] & MASK) as u32,
        (arr[0] >> 32) as u32,
        (arr[1] & MASK) as u32,
        (arr[1] >> 32) as u32,
        (arr[2] & MASK) as u32,
        (arr[2] >> 32) as u32,
        (arr[3] & MASK) as u32,
        (arr[3] >> 32) as u32,
    ])
}

#[inline(always)]
#[cfg(target_pointer_width = "64")]
const fn inner_from_raw(arr: [u64; 4]) -> Inner {
    Inner::from_words(arr)
}

#[inline]
#[cfg(target_pointer_width = "32")]
const fn inner_to_raw(uint: Inner) -> [u64; 4] {
    let words = uint.as_words();
    [
        (words[0] as u64) | ((words[1] as u64) << 32),
        (words[2] as u64) | ((words[3] as u64) << 32),
        (words[4] as u64) | ((words[5] as u64) << 32),
        (words[6] as u64) | ((words[7] as u64) << 32),
    ]
}

#[inline]
#[cfg(target_pointer_width = "64")]
const fn inner_to_raw(uint: Inner) -> [u64; 4] {
    uint.to_words()
}

#[test]
fn test_inv() {
    // Compute -(q^{-1} mod 2^{Limb::BIT_SIZE}) mod 2^{Limb::BIT_SIZE} by
    // exponentiating by totient(2**{Limb::BIT_SIZE}) - 1
    use crypto_bigint::Limb;

    let mut inv = Limb::ONE;
    for _ in 0..(Limb::BIT_SIZE - 1) {
        inv = inv.wrapping_mul(inv);
        inv = inv.wrapping_mul(FIELD.modulus.limbs()[0]);
    }
    inv = Limb::ZERO.wrapping_sub(inv);
    assert_eq!(inv, FIELD.inv);
}

#[test]
fn test_debug() {
    assert_eq!(
        format!("{:?}", Scalar::zero()),
        "0x0000000000000000000000000000000000000000000000000000000000000000"
    );
    assert_eq!(
        format!("{:?}", Scalar::one()),
        "0x0000000000000000000000000000000000000000000000000000000000000001"
    );
    assert_eq!(
        format!("{:?}", Scalar(FIELD.r2)),
        "0x1824b159acc5056f998c4fefecbc4ff55884b7fa0003480200000001fffffffe"
    );
}

#[test]
fn test_equality() {
    assert_eq!(Scalar::zero(), Scalar::zero());
    assert_eq!(Scalar::one(), Scalar::one());
    assert_eq!(FIELD.r2, FIELD.r2);

    assert!(Scalar::zero() != Scalar::one());
    assert!(Scalar::one() != Scalar(FIELD.r2));
}

#[test]
fn test_to_bytes() {
    assert_eq!(Scalar::zero().to_bytes(), [0u8; 32]);

    assert_eq!(
        Scalar::one().to_bytes(),
        [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0
        ]
    );

    assert_eq!(
        Scalar(FIELD.r2).to_bytes(),
        [
            254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88, 245, 79, 188, 236, 239,
            79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24
        ]
    );

    assert_eq!(
        (-&Scalar::one()).to_bytes(),
        [
            0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ]
    );
}

#[test]
fn test_from_bytes() {
    assert_eq!(Scalar::from_bytes(&[0u8; 32]).unwrap(), Scalar::zero());

    assert_eq!(
        Scalar::from_bytes(&[
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0
        ])
        .unwrap(),
        Scalar::one()
    );

    assert_eq!(
        Scalar::from_bytes(&[
            254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88, 245, 79, 188, 236, 239,
            79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24
        ])
        .unwrap(),
        Scalar(FIELD.r2)
    );

    // -1 should work
    assert!(bool::from(
        Scalar::from_bytes(&[
            0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_some()
    ));

    // modulus is invalid
    assert!(bool::from(
        Scalar::from_bytes(&[
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_none()
    ));

    // Anything larger than the modulus is invalid
    assert!(bool::from(
        Scalar::from_bytes(&[
            2, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_none()
    ));
    assert!(bool::from(
        Scalar::from_bytes(&[
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 58, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_none()
    ));
    assert!(bool::from(
        Scalar::from_bytes(&[
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 116
        ])
        .is_none()
    ));
}

#[test]
fn test_from_bytes_wide_r() {
    assert_eq!(
        Scalar(FIELD.one()),
        Scalar::from_bytes_wide(&[
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0
        ])
    );
}

#[test]
fn test_from_bytes_wide_r2() {
    assert_eq!(
        Scalar(FIELD.r2),
        Scalar::from_bytes_wide(&[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0
        ])
    );
}

#[test]
fn test_from_bytes_wide_negative_one() {
    assert_eq!(
        -&Scalar::one(),
        Scalar::from_bytes_wide(&[
            0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ])
    );
}

#[test]
fn test_from_bytes_wide_maximum() {
    assert_eq!(
        Scalar(FIELD.r3) - Scalar::one(),
        Scalar::from_bytes_wide(&[0xff; 64])
    );
}

#[test]
fn test_zero() {
    assert_eq!(Scalar::zero(), -&Scalar::zero());
    assert_eq!(Scalar::zero(), Scalar::zero() + Scalar::zero());
    assert_eq!(Scalar::zero(), Scalar::zero() - Scalar::zero());
    assert_eq!(Scalar::zero(), Scalar::zero() * Scalar::zero());
}

#[cfg(test)]
const LARGEST: Scalar = Scalar::from_raw_unchecked([
    0xffff_ffff_0000_0000,
    0x53bd_a402_fffe_5bfe,
    0x3339_d808_09a1_d805,
    0x73ed_a753_299d_7d48,
]);

#[test]
fn test_addition() {
    let mut tmp = LARGEST;
    tmp += &LARGEST;

    assert_eq!(
        tmp,
        Scalar::from_raw_unchecked([
            0xffff_fffe_ffff_ffff,
            0x53bd_a402_fffe_5bfe,
            0x3339_d808_09a1_d805,
            0x73ed_a753_299d_7d48,
        ])
    );

    let mut tmp = LARGEST;
    tmp += &Scalar::from_raw_unchecked([1, 0, 0, 0]);

    assert_eq!(tmp, Scalar::zero());
}

#[test]
fn test_negation() {
    let tmp = -&LARGEST;

    assert_eq!(tmp, Scalar::from_raw_unchecked([1, 0, 0, 0]));

    let tmp = -&Scalar::zero();
    assert_eq!(tmp, Scalar::zero());
    let tmp = -&Scalar::from_raw_unchecked([1, 0, 0, 0]);
    assert_eq!(tmp, LARGEST);
}

#[test]
fn test_subtraction() {
    let mut tmp = LARGEST;
    tmp -= &LARGEST;

    assert_eq!(tmp, Scalar::zero());

    let mut tmp = Scalar::zero();
    tmp -= &LARGEST;

    let mut tmp2 = Scalar(FIELD.modulus);
    tmp2 -= &LARGEST;

    assert_eq!(tmp, tmp2);
}

#[test]
fn test_multiplication() {
    let mut cur = LARGEST;

    for _ in 0..100 {
        let mut tmp = cur;
        tmp *= &cur;

        let mut tmp2 = Scalar::zero();
        for b in cur
            .to_bytes()
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
        {
            let tmp3 = tmp2;
            tmp2 += tmp3;

            if b {
                tmp2 += cur;
            }
        }

        assert_eq!(tmp, tmp2);

        cur += LARGEST;
    }
}

#[test]
fn test_squaring() {
    let mut cur = LARGEST;

    for _ in 0..100 {
        let mut tmp = cur;
        tmp = tmp.square();

        let mut tmp2 = Scalar::zero();
        for b in cur
            .to_bytes()
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
        {
            let tmp3 = tmp2;
            tmp2 += tmp3;

            if b {
                tmp2 += &cur;
            }
        }

        assert_eq!(tmp, tmp2);

        cur += LARGEST;
    }
}

#[test]
fn test_inversion() {
    assert!(bool::from(Scalar::zero().invert().is_none()));
    assert_eq!(Scalar::one().invert().unwrap(), Scalar::one());
    assert_eq!((-&Scalar::one()).invert().unwrap(), -&Scalar::one());

    let mut tmp = Scalar(FIELD.r2);

    for _ in 0..100 {
        let mut tmp2 = tmp.invert().unwrap();
        tmp2 *= tmp;

        assert_eq!(tmp2, Scalar::one());

        tmp += Scalar(FIELD.r2);
    }

    assert!(Scalar::zero().invert_vartime().is_none());
    assert_eq!(Scalar::one().invert_vartime(), Some(Scalar::one()));
    assert_eq!((-&Scalar::one()).invert_vartime(), Some(-&Scalar::one()));

    let mut tmp = Scalar(FIELD.r2);

    for _ in 0..100 {
        let mut tmp2 = tmp.invert_vartime().unwrap();
        tmp2 *= tmp;

        assert_eq!(tmp2, Scalar::one());

        tmp += Scalar(FIELD.r2);
    }
}

#[test]
fn test_batch_inversion() {
    let a = Scalar::from_raw_unchecked([1, 2, 3, 4]);
    let ai = a.invert().unwrap();
    let b = Scalar::from_raw_unchecked([5, 6, 7, 8]);
    let bi = b.invert().unwrap();
    let tests = &[(
        vec![a, b, Scalar::one(), Scalar::zero()],
        vec![ai, bi, Scalar::one(), Scalar::zero()],
    )];

    for (batch, expect) in tests {
        let mut out = vec![Scalar::zero(); batch.len()];
        Scalar::batch_invert(&batch, &mut out);
        assert_eq!(&out, expect);

        let mut out = batch.clone();
        Scalar::batch_invert_in_place(&mut out);
        assert_eq!(&out, expect);
    }
}

#[test]
fn test_invert_is_pow() {
    let q_minus_2 = [
        0xffff_fffe_ffff_ffff,
        0x53bd_a402_fffe_5bfe,
        0x3339_d808_09a1_d805,
        0x73ed_a753_299d_7d48,
    ];

    let mut r1 = Scalar::one();
    let mut r2 = Scalar::one();
    let mut r3 = Scalar::one();

    for i in 0..100 {
        r1 = r1.invert().unwrap();
        r2 = r2.pow_vartime(&q_minus_2);
        r3 = r3.pow(&q_minus_2);

        assert_eq!(r1, r2, "failed on {}", i);
        assert_eq!(r2, r3);
        // Add R so we check something different next time around
        r1 += Scalar::one();
        r2 = r1;
        r3 = r1;
    }
}

#[test]
fn test_sqrt() {
    {
        assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
    }

    let mut square = Scalar::from_raw_unchecked([
        0x46cd_85a5_f273_077e,
        0x1d30_c47d_d68f_c735,
        0x77f6_56f6_0bec_a0eb,
        0x494a_a01b_df32_468d,
    ]);

    let mut none_count = 0;

    for _ in 0..100 {
        let square_root = square.sqrt();
        if bool::from(square_root.is_none()) {
            none_count += 1;
        } else {
            assert_eq!(square_root.unwrap() * square_root.unwrap(), square);
        }
        square -= Scalar::one();
    }

    assert_eq!(49, none_count);
}

#[test]
fn test_from_raw() {
    assert_eq!(
        Scalar::from_raw([
            0x0001_ffff_fffd,
            0x5884_b7fa_0003_4802,
            0x998c_4fef_ecbc_4ff5,
            0x1824_b159_acc5_056f,
        ]),
        Scalar::from_raw([u64::MAX; 4])
    );

    assert_eq!(
        Scalar::from_raw(inner_to_raw(FIELD.modulus)),
        Scalar::zero()
    );

    assert_eq!(Scalar::from_raw([1, 0, 0, 0]), Scalar::one());
}

#[test]
fn test_from_canonical() {
    assert_eq!(
        Scalar::from_raw([
            0x0001_ffff_fffd,
            0x5884_b7fa_0003_4802,
            0x998c_4fef_ecbc_4ff5,
            0x1824_b159_acc5_056f,
        ]),
        Scalar::from_canonical(Inner::MAX)
    );

    assert_eq!(Scalar::from_canonical(FIELD.modulus), Scalar::zero());

    assert_eq!(Scalar::from_canonical(Inner::ONE), Scalar::one());
}

#[test]
fn test_double() {
    let a = Scalar::from_raw([
        0x1fff_3231_233f_fffd,
        0x4884_b7fa_0003_4802,
        0x998c_4fef_ecbc_4ff3,
        0x1824_b159_acc5_0562,
    ]);

    assert_eq!(a.double(), a + a);
}

#[cfg(feature = "zeroize")]
#[test]
fn test_zeroize() {
    use zeroize::Zeroize;

    let mut a = Scalar::from_raw([
        0x1fff_3231_233f_fffd,
        0x4884_b7fa_0003_4802,
        0x998c_4fef_ecbc_4ff3,
        0x1824_b159_acc5_0562,
    ]);
    a.zeroize();
    assert!(bool::from(a.is_zero()));
}
