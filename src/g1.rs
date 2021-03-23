//! This module provides an implementation of the $\mathbb{G}_1$ group of BLS12-381.

use core::borrow::Borrow;
use core::fmt;
use core::iter::Sum;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use group::{
    prime::{PrimeCurve, PrimeCurveAffine, PrimeGroup},
    Curve, Group, GroupEncoding, UncompressedEncoding,
};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "alloc")]
use group::WnafGroup;

use crate::fp::Fp;
use crate::hash_to_curve::IsogenyMap;
use crate::Scalar;

/// This is an element of $\mathbb{G}_1$ represented in the affine coordinate space.
/// It is ideal to keep elements in this representation to reduce memory usage and
/// improve performance through the use of mixed curve model arithmetic.
///
/// Values of `G1Affine` are guaranteed to be in the $q$-order subgroup unless an
/// "unchecked" API was misused.
#[cfg_attr(docsrs, doc(cfg(feature = "groups")))]
#[derive(Copy, Clone, Debug)]
pub struct G1Affine {
    pub(crate) x: Fp,
    pub(crate) y: Fp,
    infinity: Choice,
}

impl Default for G1Affine {
    fn default() -> G1Affine {
        G1Affine::identity()
    }
}

impl fmt::Display for G1Affine {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<'a> From<&'a G1Projective> for G1Affine {
    fn from(p: &'a G1Projective) -> G1Affine {
        let zinv = p.z.invert().unwrap_or(Fp::zero());
        let x = p.x * zinv;
        let y = p.y * zinv;

        let tmp = G1Affine {
            x,
            y,
            infinity: Choice::from(0u8),
        };

        G1Affine::conditional_select(&tmp, &G1Affine::identity(), zinv.is_zero())
    }
}

impl From<G1Projective> for G1Affine {
    fn from(p: G1Projective) -> G1Affine {
        G1Affine::from(&p)
    }
}

impl ConstantTimeEq for G1Affine {
    fn ct_eq(&self, other: &Self) -> Choice {
        // The only cases in which two points are equal are
        // 1. infinity is set on both
        // 2. infinity is not set on both, and their coordinates are equal

        (self.infinity & other.infinity)
            | ((!self.infinity)
                & (!other.infinity)
                & self.x.ct_eq(&other.x)
                & self.y.ct_eq(&other.y))
    }
}

impl ConditionallySelectable for G1Affine {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        G1Affine {
            x: Fp::conditional_select(&a.x, &b.x, choice),
            y: Fp::conditional_select(&a.y, &b.y, choice),
            infinity: Choice::conditional_select(&a.infinity, &b.infinity, choice),
        }
    }
}

impl Eq for G1Affine {}
impl PartialEq for G1Affine {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a G1Affine {
    type Output = G1Affine;

    #[inline]
    fn neg(self) -> G1Affine {
        G1Affine {
            x: self.x,
            y: Fp::conditional_select(&-self.y, &Fp::one(), self.infinity),
            infinity: self.infinity,
        }
    }
}

impl Neg for G1Affine {
    type Output = G1Affine;

    #[inline]
    fn neg(self) -> G1Affine {
        -&self
    }
}

impl<'a, 'b> Add<&'b G1Projective> for &'a G1Affine {
    type Output = G1Projective;

    #[inline]
    fn add(self, rhs: &'b G1Projective) -> G1Projective {
        rhs.add_mixed(self)
    }
}

impl<'a, 'b> Add<&'b G1Affine> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn add(self, rhs: &'b G1Affine) -> G1Projective {
        self.add_mixed(rhs)
    }
}

impl<'a, 'b> Sub<&'b G1Projective> for &'a G1Affine {
    type Output = G1Projective;

    #[inline]
    fn sub(self, rhs: &'b G1Projective) -> G1Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Sub<&'b G1Affine> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn sub(self, rhs: &'b G1Affine) -> G1Projective {
        self + (-rhs)
    }
}

impl<T> Sum<T> for G1Projective
where
    T: Borrow<G1Projective>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::identity(), |acc, item| acc + item.borrow())
    }
}

impl_binops_additive!(G1Projective, G1Affine);
impl_binops_additive_specify_output!(G1Affine, G1Projective, G1Projective);

const B: Fp = Fp::from_raw_unchecked([
    0xaa27_0000_000c_fff3,
    0x53cc_0032_fc34_000a,
    0x478f_e97a_6b0a_807f,
    0xb1d3_7ebe_e6ba_24d7,
    0x8ec9_733b_bf78_ab2f,
    0x09d6_4551_3d83_de7e,
]);

const G1_GENERATOR_X: Fp = Fp::from_raw_unchecked([
    0x5cb3_8790_fd53_0c16,
    0x7817_fc67_9976_fff5,
    0x154f_95c7_143b_a1c1,
    0xf0ae_6acd_f3d0_e747,
    0xedce_6ecc_21db_f440,
    0x1201_7741_9e0b_fb75,
]);

const G1_GENERATOR_Y: Fp = Fp::from_raw_unchecked([
    0xbaac_93d5_0ce7_2271,
    0x8c22_631a_7918_fd8e,
    0xdd59_5f13_5707_25ce,
    0x51ac_5829_5040_5194,
    0x0e1c_8c3f_ad00_59c0,
    0x0bbc_3efc_5008_a26a,
]);

impl G1Affine {
    /// Returns the identity of the group: the point at infinity.
    pub fn identity() -> G1Affine {
        G1Affine {
            x: Fp::zero(),
            y: Fp::one(),
            infinity: Choice::from(1u8),
        }
    }

    /// Returns a fixed generator of the group. See [`notes::design`](notes/design/index.html#fixed-generators)
    /// for how this generator is chosen.
    pub fn generator() -> G1Affine {
        G1Affine {
            x: G1_GENERATOR_X,
            y: G1_GENERATOR_Y,
            infinity: Choice::from(0u8),
        }
    }

    /// Serializes this element into compressed form. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn to_compressed(&self) -> [u8; 48] {
        // Strictly speaking, self.x is zero already when self.infinity is true, but
        // to guard against implementation mistakes we do not assume this.
        let mut res = Fp::conditional_select(&self.x, &Fp::zero(), self.infinity).to_bytes();

        // This point is in compressed form, so we set the most significant bit.
        res[0] |= 1u8 << 7;

        // Is this point at infinity? If so, set the second-most significant bit.
        res[0] |= u8::conditional_select(&0u8, &(1u8 << 6), self.infinity);

        // Is the y-coordinate the lexicographically largest of the two associated with the
        // x-coordinate? If so, set the third-most significant bit so long as this is not
        // the point at infinity.
        res[0] |= u8::conditional_select(
            &0u8,
            &(1u8 << 5),
            (!self.infinity) & self.y.lexicographically_largest(),
        );

        res
    }

    /// Serializes this element into uncompressed form. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn to_uncompressed(&self) -> [u8; 96] {
        let mut res = [0; 96];

        res[0..48].copy_from_slice(
            &Fp::conditional_select(&self.x, &Fp::zero(), self.infinity).to_bytes()[..],
        );
        res[48..96].copy_from_slice(
            &Fp::conditional_select(&self.y, &Fp::zero(), self.infinity).to_bytes()[..],
        );

        // Is this point at infinity? If so, set the second-most significant bit.
        res[0] |= u8::conditional_select(&0u8, &(1u8 << 6), self.infinity);

        res
    }

    /// Attempts to deserialize an uncompressed element. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn from_uncompressed(bytes: &[u8; 96]) -> CtOption<Self> {
        Self::from_uncompressed_unchecked(bytes)
            .and_then(|p| CtOption::new(p, p.is_on_curve() & p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 96]) -> CtOption<Self> {
        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = Choice::from((bytes[0] >> 7) & 1);
        let infinity_flag_set = Choice::from((bytes[0] >> 6) & 1);
        let sort_flag_set = Choice::from((bytes[0] >> 5) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[0] &= 0b0001_1111;

            Fp::from_bytes(&tmp)
        };

        // Attempt to obtain the y-coordinate
        let y = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[48..96]);

            Fp::from_bytes(&tmp)
        };

        x.and_then(|x| {
            y.and_then(|y| {
                // Create a point representing this value
                let p = G1Affine::conditional_select(
                    &G1Affine {
                        x,
                        y,
                        infinity: infinity_flag_set,
                    },
                    &G1Affine::identity(),
                    infinity_flag_set,
                );

                CtOption::new(
                    p,
                    // If the infinity flag is set, the x and y coordinates should have been zero.
                    ((!infinity_flag_set) | (infinity_flag_set & x.is_zero() & y.is_zero())) &
                    // The compression flag should not have been set, as this is an uncompressed element
                    (!compression_flag_set) &
                    // The sort flag should not have been set, as this is an uncompressed element
                    (!sort_flag_set),
                )
            })
        })
    }

    /// Attempts to deserialize a compressed element. See [`notes::serialization`](crate::notes::serialization)
    /// for details about how group elements are serialized.
    pub fn from_compressed(bytes: &[u8; 48]) -> CtOption<Self> {
        // We already know the point is on the curve because this is established
        // by the y-coordinate recovery procedure in from_compressed_unchecked().

        Self::from_compressed_unchecked(bytes).and_then(|p| CtOption::new(p, p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(bytes: &[u8; 48]) -> CtOption<Self> {
        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = Choice::from((bytes[0] >> 7) & 1);
        let infinity_flag_set = Choice::from((bytes[0] >> 6) & 1);
        let sort_flag_set = Choice::from((bytes[0] >> 5) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[0] &= 0b0001_1111;

            Fp::from_bytes(&tmp)
        };

        x.and_then(|x| {
            // If the infinity flag is set, return the value assuming
            // the x-coordinate is zero and the sort bit is not set.
            //
            // Otherwise, return a recovered point (assuming the correct
            // y-coordinate can be found) so long as the infinity flag
            // was not set.
            CtOption::new(
                G1Affine::identity(),
                infinity_flag_set & // Infinity flag should be set
                compression_flag_set & // Compression flag should be set
                (!sort_flag_set) & // Sort flag should not be set
                x.is_zero(), // The x-coordinate should be zero
            )
            .or_else(|| {
                // Recover a y-coordinate given x by y = sqrt(x^3 + 4)
                ((x.square() * x) + B).sqrt().and_then(|y| {
                    // Switch to the correct y-coordinate if necessary.
                    let y = Fp::conditional_select(
                        &y,
                        &-y,
                        y.lexicographically_largest() ^ sort_flag_set,
                    );

                    CtOption::new(
                        G1Affine {
                            x,
                            y,
                            infinity: infinity_flag_set,
                        },
                        (!infinity_flag_set) & // Infinity flag should not be set
                        compression_flag_set, // Compression flag should be set
                    )
                })
            })
        })
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.infinity
    }

    /// Returns true if this point is free of an $h$-torsion component, and so it
    /// exists within the $q$-order subgroup $\mathbb{G}_1$. This should always return true
    /// unless an "unchecked" API was used.
    pub fn is_torsion_free(&self) -> Choice {
        const FQ_MODULUS_BYTES: [u8; 32] = [
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115,
        ];

        // Clear the r-torsion from the point and check if it is the identity
        G1Projective::from(*self)
            .multiply(&FQ_MODULUS_BYTES)
            .is_identity()
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // y^2 - x^3 ?= 4
        (self.y.square() - (self.x.square() * self.x)).ct_eq(&B) | self.infinity
    }
}

/// This is an element of $\mathbb{G}_1$ represented in the projective coordinate space.
#[cfg_attr(docsrs, doc(cfg(feature = "groups")))]
#[derive(Copy, Clone, Debug)]
pub struct G1Projective {
    x: Fp,
    y: Fp,
    z: Fp,
}

impl Default for G1Projective {
    fn default() -> G1Projective {
        G1Projective::identity()
    }
}

impl fmt::Display for G1Projective {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<'a> From<&'a G1Affine> for G1Projective {
    fn from(p: &'a G1Affine) -> G1Projective {
        G1Projective {
            x: p.x,
            y: p.y,
            z: Fp::conditional_select(&Fp::one(), &Fp::zero(), p.infinity),
        }
    }
}

impl From<G1Affine> for G1Projective {
    fn from(p: G1Affine) -> G1Projective {
        G1Projective::from(&p)
    }
}

impl ConstantTimeEq for G1Projective {
    fn ct_eq(&self, other: &Self) -> Choice {
        // Is (xz, yz, z) equal to (x'z', y'z', z') when converted to affine?

        let x1 = self.x * other.z;
        let x2 = other.x * self.z;

        let y1 = self.y * other.z;
        let y2 = other.y * self.z;

        let self_is_zero = self.z.is_zero();
        let other_is_zero = other.z.is_zero();

        (self_is_zero & other_is_zero) // Both point at infinity
            | ((!self_is_zero) & (!other_is_zero) & x1.ct_eq(&x2) & y1.ct_eq(&y2))
        // Neither point at infinity, coordinates are the same
    }
}

impl ConditionallySelectable for G1Projective {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        G1Projective {
            x: Fp::conditional_select(&a.x, &b.x, choice),
            y: Fp::conditional_select(&a.y, &b.y, choice),
            z: Fp::conditional_select(&a.z, &b.z, choice),
        }
    }
}

impl Eq for G1Projective {}
impl PartialEq for G1Projective {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn neg(self) -> G1Projective {
        G1Projective {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
    }
}

impl Neg for G1Projective {
    type Output = G1Projective;

    #[inline]
    fn neg(self) -> G1Projective {
        -&self
    }
}

impl<'a, 'b> Add<&'b G1Projective> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn add(self, rhs: &'b G1Projective) -> G1Projective {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b G1Projective> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn sub(self, rhs: &'b G1Projective) -> G1Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G1Projective {
    type Output = G1Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        self.multiply(&other.to_bytes())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G1Affine {
    type Output = G1Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        G1Projective::from(self).multiply(&other.to_bytes())
    }
}

impl_binops_additive!(G1Projective, G1Projective);
impl_binops_multiplicative!(G1Projective, Scalar);
impl_binops_multiplicative_mixed!(G1Affine, Scalar, G1Projective);

#[inline(always)]
fn mul_by_3b(a: Fp) -> Fp {
    let a = a + a; // 2
    let a = a + a; // 4
    let a = a + a + a; // 12
    a
}

impl G1Projective {
    /// Returns the identity of the group: the point at infinity.
    pub fn identity() -> G1Projective {
        G1Projective {
            x: Fp::zero(),
            y: Fp::one(),
            z: Fp::zero(),
        }
    }

    /// Returns a fixed generator of the group. See [`notes::design`](notes/design/index.html#fixed-generators)
    /// for how this generator is chosen.
    pub fn generator() -> G1Projective {
        G1Projective {
            x: G1_GENERATOR_X,
            y: G1_GENERATOR_Y,
            z: Fp::one(),
        }
    }

    /// Computes the doubling of this point.
    pub fn double(&self) -> G1Projective {
        // Algorithm 9, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.y.square();
        let z3 = t0 + t0;
        let z3 = z3 + z3;
        let z3 = z3 + z3;
        let t1 = self.y * self.z;
        let t2 = self.z.square();
        let t2 = mul_by_3b(t2);
        let x3 = t2 * z3;
        let y3 = t0 + t2;
        let z3 = t1 * z3;
        let t1 = t2 + t2;
        let t2 = t1 + t2;
        let t0 = t0 - t2;
        let y3 = t0 * y3;
        let y3 = x3 + y3;
        let t1 = self.x * self.y;
        let x3 = t0 * t1;
        let x3 = x3 + x3;

        let tmp = G1Projective {
            x: x3,
            y: y3,
            z: z3,
        };

        G1Projective::conditional_select(&tmp, &G1Projective::identity(), self.is_identity())
    }

    /// Adds this point to another point.
    pub fn add(&self, rhs: &G1Projective) -> G1Projective {
        // Algorithm 7, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t2 = self.z * rhs.z;
        let t3 = self.x + self.y;
        let t4 = rhs.x + rhs.y;
        let t3 = t3 * t4;
        let t4 = t0 + t1;
        let t3 = t3 - t4;
        let t4 = self.y + self.z;
        let x3 = rhs.y + rhs.z;
        let t4 = t4 * x3;
        let x3 = t1 + t2;
        let t4 = t4 - x3;
        let x3 = self.x + self.z;
        let y3 = rhs.x + rhs.z;
        let x3 = x3 * y3;
        let y3 = t0 + t2;
        let y3 = x3 - y3;
        let x3 = t0 + t0;
        let t0 = x3 + t0;
        let t2 = mul_by_3b(t2);
        let z3 = t1 + t2;
        let t1 = t1 - t2;
        let y3 = mul_by_3b(y3);
        let x3 = t4 * y3;
        let t2 = t3 * t1;
        let x3 = t2 - x3;
        let y3 = y3 * t0;
        let t1 = t1 * z3;
        let y3 = t1 + y3;
        let t0 = t0 * t3;
        let z3 = z3 * t4;
        let z3 = z3 + t0;

        G1Projective {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &G1Affine) -> G1Projective {
        // Algorithm 8, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t3 = rhs.x + rhs.y;
        let t4 = self.x + self.y;
        let t3 = t3 * t4;
        let t4 = t0 + t1;
        let t3 = t3 - t4;
        let t4 = rhs.y * self.z;
        let t4 = t4 + self.y;
        let y3 = rhs.x * self.z;
        let y3 = y3 + self.x;
        let x3 = t0 + t0;
        let t0 = x3 + t0;
        let t2 = mul_by_3b(self.z);
        let z3 = t1 + t2;
        let t1 = t1 - t2;
        let y3 = mul_by_3b(y3);
        let x3 = t4 * y3;
        let t2 = t3 * t1;
        let x3 = t2 - x3;
        let y3 = y3 * t0;
        let t1 = t1 * z3;
        let y3 = t1 + y3;
        let t0 = t0 * t3;
        let z3 = z3 * t4;
        let z3 = z3 + t0;

        let tmp = G1Projective {
            x: x3,
            y: y3,
            z: z3,
        };

        G1Projective::conditional_select(&tmp, &self, rhs.is_identity())
    }

    fn multiply(&self, by: &[u8; 32]) -> G1Projective {
        let mut acc = G1Projective::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the leading bit because it's always unset for Fq
        // elements.
        for bit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = G1Projective::conditional_select(&acc, &(acc + self), bit);
        }

        acc
    }

    /// Multiply `self` by `crate::BLS_X`, using double and add.
    fn mul_by_x(&self) -> G1Projective {
        let mut xself = G1Projective::identity();
        // NOTE: in BLS12-381 we can just skip the first bit.
        let mut x = crate::BLS_X >> 1;
        let mut tmp = *self;
        while x != 0 {
            tmp = tmp.double();

            if x % 2 == 1 {
                xself += tmp;
            }
            x >>= 1;
        }
        // finally, flip the sign
        if crate::BLS_X_IS_NEGATIVE {
            xself = -xself;
        }
        xself
    }

    /// Multiplies by $(1 - z)$, where $z$ is the parameter of BLS12-381, which
    /// [suffices to clear](https://ia.cr/2019/403) the cofactor and map
    /// elliptic curve points to elements of $\mathbb{G}\_1$.
    pub fn clear_cofactor(&self) -> G1Projective {
        self - self.mul_by_x()
    }

    /// Converts a batch of `G1Projective` elements into `G1Affine` elements. This
    /// function will panic if `p.len() != q.len()`.
    pub fn batch_normalize(p: &[Self], q: &mut [G1Affine]) {
        assert_eq!(p.len(), q.len());

        let mut acc = Fp::one();
        for (p, q) in p.iter().zip(q.iter_mut()) {
            // We use the `x` field of `G1Affine` to store the product
            // of previous z-coordinates seen.
            q.x = acc;

            // We will end up skipping all identities in p
            acc = Fp::conditional_select(&(acc * p.z), &acc, p.is_identity());
        }

        // This is the inverse, as all z-coordinates are nonzero and the ones
        // that are not are skipped.
        acc = acc.invert().unwrap();

        for (p, q) in p.iter().rev().zip(q.iter_mut().rev()) {
            let skip = p.is_identity();

            // Compute tmp = 1/z
            let tmp = q.x * acc;

            // Cancel out z-coordinate in denominator of `acc`
            acc = Fp::conditional_select(&(acc * p.z), &acc, skip);

            // Set the coordinates to the correct value
            q.x = p.x * tmp;
            q.y = p.y * tmp;
            q.infinity = Choice::from(0u8);

            *q = G1Affine::conditional_select(&q, &G1Affine::identity(), skip);
        }
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.z.is_zero()
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // Y^2 Z = X^3 + b Z^3

        (self.y.square() * self.z).ct_eq(&(self.x.square() * self.x + self.z.square() * self.z * B))
            | self.z.is_zero()
    }

    /// Isogeny evaluation function.
    pub(crate) fn eval_iso(&mut self, coeffs: [&[Fp]; 4]) {
        let mut tmp = [Fp::zero(); 16];
        let mut mapvals = [Fp::zero(); 4];

        // unpack input point
        let x = &mut self.x;
        let y = &mut self.y;
        let z = &mut self.z;

        // precompute powers of z
        let zpows = {
            let mut zpows = [Fp::zero(); 15];
            zpows[0] = z.square(); // z^2
            zpows[1] = zpows[0].square(); // z^4
            {
                let (z_squared, rest) = zpows.split_at_mut(1);
                for idx in 1..coeffs[2].len() - 2 {
                    if idx % 2 == 0 {
                        rest[idx] = rest[idx / 2 - 1].square();
                    } else {
                        rest[idx] = rest[idx - 1];
                        rest[idx].mul_assign(&z_squared[0]);
                    }
                }
            }
            zpows
        };

        for idx in 0..4 {
            let clen = coeffs[idx].len() - 1;
            // multiply coeffs by powers of Z
            for jdx in 0..clen {
                tmp[jdx] = coeffs[idx][clen - 1 - jdx];
                tmp[jdx].mul_assign(&zpows[jdx]);
            }
            // compute map value by Horner's rule
            mapvals[idx] = coeffs[idx][clen];
            for tmpval in &tmp[..clen] {
                mapvals[idx].mul_assign(&*x);
                mapvals[idx].add_assign(tmpval);
            }
        }

        // x denominator is order 1 less than x numerator, so we need an extra factor of Z^2
        mapvals[1].mul_assign(&zpows[0]);

        // multiply result of Y map by the y-coord, y / z^3
        mapvals[2].mul_assign(&*y);
        mapvals[3].mul_assign(&*z);
        mapvals[3].mul_assign(&zpows[0]);

        // compute Jacobian coordinates of resulting point
        *z = mapvals[1] * mapvals[3]; // Zout = xden * yden

        *x = mapvals[0];
        x.mul_assign(&mapvals[3]); // xnum * yden
        x.mul_assign(&*z); // xnum * xden * yden^2

        *y = z.square(); // xden^2 * yden^2
        y.mul_assign(&mapvals[2]); // ynum * xden^2 * yden^2
        y.mul_assign(&mapvals[1]); // ynum * xden^3 * yden^2
    }
}

/// Coefficients of the 11-isogeny x map's numerator
const XNUM: [Fp; 12] = [
    Fp::from_raw_unchecked([
        0x4d18b6f3af00131c,
        0x19fa219793fee28c,
        0x3f2885f1467f19ae,
        0x23dcea34f2ffb304,
        0xd15b58d2ffc00054,
        0x0913be200a20bef4,
    ]),
    Fp::from_raw_unchecked([
        0x898985385cdbbd8b,
        0x3c79e43cc7d966aa,
        0x1597e193f4cd233a,
        0x8637ef1e4d6623ad,
        0x11b22deed20d827b,
        0x07097bc5998784ad,
    ]),
    Fp::from_raw_unchecked([
        0xa542583a480b664b,
        0xfc7169c026e568c6,
        0x5ba2ef314ed8b5a6,
        0x5b5491c05102f0e7,
        0xdf6e99707d2a0079,
        0x0784151ed7605524,
    ]),
    Fp::from_raw_unchecked([
        0x494e212870f72741,
        0xab9be52fbda43021,
        0x26f5577994e34c3d,
        0x049dfee82aefbd60,
        0x65dadd7828505289,
        0x0e93d431ea011aeb,
    ]),
    Fp::from_raw_unchecked([
        0x90ee774bd6a74d45,
        0x7ada1c8a41bfb185,
        0x0f1a8953b325f464,
        0x104c24211be4805c,
        0x169139d319ea7a8f,
        0x09f20ead8e532bf6,
    ]),
    Fp::from_raw_unchecked([
        0x6ddd93e2f43626b7,
        0xa5482c9aa1ccd7bd,
        0x143245631883f4bd,
        0x2e0a94ccf77ec0db,
        0xb0282d480e56489f,
        0x18f4bfcbb4368929,
    ]),
    Fp::from_raw_unchecked([
        0x23c5f0c953402dfd,
        0x7a43ff6958ce4fe9,
        0x2c390d3d2da5df63,
        0xd0df5c98e1f9d70f,
        0xffd89869a572b297,
        0x1277ffc72f25e8fe,
    ]),
    Fp::from_raw_unchecked([
        0x79f4f0490f06a8a6,
        0x85f894a88030fd81,
        0x12da3054b18b6410,
        0xe2a57f6505880d65,
        0xbba074f260e400f1,
        0x08b76279f621d028,
    ]),
    Fp::from_raw_unchecked([
        0xe67245ba78d5b00b,
        0x8456ba9a1f186475,
        0x7888bff6e6b33bb4,
        0xe21585b9a30f86cb,
        0x05a69cdcef55feee,
        0x09e699dd9adfa5ac,
    ]),
    Fp::from_raw_unchecked([
        0x0de5c357bff57107,
        0x0a0db4ae6b1a10b2,
        0xe256bb67b3b3cd8d,
        0x8ad456574e9db24f,
        0x0443915f50fd4179,
        0x098c4bf7de8b6375,
    ]),
    Fp::from_raw_unchecked([
        0xe6b0617e7dd929c7,
        0xfe6e37d442537375,
        0x1dafdeda137a489e,
        0xe4efd1ad3f767ceb,
        0x4a51d8667f0fe1cf,
        0x054fdf4bbf1d821c,
    ]),
    Fp::from_raw_unchecked([
        0x72db2a50658d767b,
        0x8abf91faa257b3d5,
        0xe969d6833764ab47,
        0x464170142a1009eb,
        0xb14f01aadb30be2f,
        0x18ae6a856f40715d,
    ]),
];

/// Coefficients of the 11-isogeny x map's denominator
const XDEN: [Fp; 11] = [
    Fp::from_raw_unchecked([
        0xb962a077fdb0f945,
        0xa6a9740fefda13a0,
        0xc14d568c3ed6c544,
        0xb43fc37b908b133e,
        0x9c0b3ac929599016,
        0x0165aa6c93ad115f,
    ]),
    Fp::from_raw_unchecked([
        0x23279a3ba506c1d9,
        0x92cfca0a9465176a,
        0x3b294ab13755f0ff,
        0x116dda1c5070ae93,
        0xed4530924cec2045,
        0x083383d6ed81f1ce,
    ]),
    Fp::from_raw_unchecked([
        0x9885c2a6449fecfc,
        0x4a2b54ccd37733f0,
        0x17da9ffd8738c142,
        0xa0fba72732b3fafd,
        0xff364f36e54b6812,
        0x0f29c13c660523e2,
    ]),
    Fp::from_raw_unchecked([
        0xe349cc118278f041,
        0xd487228f2f3204fb,
        0xc9d325849ade5150,
        0x43a92bd69c15c2df,
        0x1c2c7844bc417be4,
        0x12025184f407440c,
    ]),
    Fp::from_raw_unchecked([
        0x587f65ae6acb057b,
        0x1444ef325140201f,
        0xfbf995e71270da49,
        0xccda066072436a42,
        0x7408904f0f186bb2,
        0x13b93c63edf6c015,
    ]),
    Fp::from_raw_unchecked([
        0xfb918622cd141920,
        0x4a4c64423ecaddb4,
        0x0beb232927f7fb26,
        0x30f94df6f83a3dc2,
        0xaeedd424d780f388,
        0x06cc402dd594bbeb,
    ]),
    Fp::from_raw_unchecked([
        0xd41f761151b23f8f,
        0x32a92465435719b3,
        0x64f436e888c62cb9,
        0xdf70a9a1f757c6e4,
        0x6933a38d5b594c81,
        0x0c6f7f7237b46606,
    ]),
    Fp::from_raw_unchecked([
        0x693c08747876c8f7,
        0x22c9850bf9cf80f0,
        0x8e9071dab950c124,
        0x89bc62d61c7baf23,
        0xbc6be2d8dad57c23,
        0x17916987aa14a122,
    ]),
    Fp::from_raw_unchecked([
        0x1be3ff439c1316fd,
        0x9965243a7571dfa7,
        0xc7f7f62962f5cd81,
        0x32c6aa9af394361c,
        0xbbc2ee18e1c227f4,
        0x0c102cbac531bb34,
    ]),
    Fp::from_raw_unchecked([
        0x997614c97bacbf07,
        0x61f86372b99192c0,
        0x5b8c95fc14353fc3,
        0xca2b066c2a87492f,
        0x16178f5bbf698711,
        0x12a6dcd7f0f4e0e8,
    ]),
    Fp::from_raw_unchecked([
        0x760900000002fffd,
        0xebf4000bc40c0002,
        0x5f48985753c758ba,
        0x77ce585370525745,
        0x5c071a97a256ec6d,
        0x15f65ec3fa80e493,
    ]),
];

/// Coefficients of the 11-isogeny y map's numerator
const YNUM: [Fp; 16] = [
    Fp::from_raw_unchecked([
        0x2b567ff3e2837267,
        0x1d4d9e57b958a767,
        0xce028fea04bd7373,
        0xcc31a30a0b6cd3df,
        0x7d7b18a682692693,
        0x0d300744d42a0310,
    ]),
    Fp::from_raw_unchecked([
        0x99c2555fa542493f,
        0xfe7f53cc4874f878,
        0x5df0608b8f97608a,
        0x14e03832052b49c8,
        0x706326a6957dd5a4,
        0x0a8dadd9c2414555,
    ]),
    Fp::from_raw_unchecked([
        0x13d942922a5cf63a,
        0x357e33e36e261e7d,
        0xcf05a27c8456088d,
        0x0000bd1de7ba50f0,
        0x83d0c7532f8c1fde,
        0x13f70bf38bbf2905,
    ]),
    Fp::from_raw_unchecked([
        0x5c57fd95bfafbdbb,
        0x28a359a65e541707,
        0x3983ceb4f6360b6d,
        0xafe19ff6f97e6d53,
        0xb3468f4550192bf7,
        0x0bb6cde49d8ba257,
    ]),
    Fp::from_raw_unchecked([
        0x590b62c7ff8a513f,
        0x314b4ce372cacefd,
        0x6bef32ce94b8a800,
        0x6ddf84a095713d5f,
        0x64eace4cb0982191,
        0x0386213c651b888d,
    ]),
    Fp::from_raw_unchecked([
        0xa5310a31111bbcdd,
        0xa14ac0f5da148982,
        0xf9ad9cc95423d2e9,
        0xaa6ec095283ee4a7,
        0xcf5b1f022e1c9107,
        0x01fddf5aed881793,
    ]),
    Fp::from_raw_unchecked([
        0x65a572b0d7a7d950,
        0xe25c2d8183473a19,
        0xc2fcebe7cb877dbd,
        0x05b2d36c769a89b0,
        0xba12961be86e9efb,
        0x07eb1b29c1dfde1f,
    ]),
    Fp::from_raw_unchecked([
        0x93e09572f7c4cd24,
        0x364e929076795091,
        0x8569467e68af51b5,
        0xa47da89439f5340f,
        0xf4fa918082e44d64,
        0x0ad52ba3e6695a79,
    ]),
    Fp::from_raw_unchecked([
        0x911429844e0d5f54,
        0xd03f51a3516bb233,
        0x3d587e5640536e66,
        0xfa86d2a3a9a73482,
        0xa90ed5adf1ed5537,
        0x149c9c326a5e7393,
    ]),
    Fp::from_raw_unchecked([
        0x462bbeb03c12921a,
        0xdc9af5fa0a274a17,
        0x9a558ebde836ebed,
        0x649ef8f11a4fae46,
        0x8100e1652b3cdc62,
        0x1862bd62c291dacb,
    ]),
    Fp::from_raw_unchecked([
        0x05c9b8ca89f12c26,
        0x0194160fa9b9ac4f,
        0x6a643d5a6879fa2c,
        0x14665bdd8846e19d,
        0xbb1d0d53af3ff6bf,
        0x12c7e1c3b28962e5,
    ]),
    Fp::from_raw_unchecked([
        0xb55ebf900b8a3e17,
        0xfedc77ec1a9201c4,
        0x1f07db10ea1a4df4,
        0x0dfbd15dc41a594d,
        0x389547f2334a5391,
        0x02419f98165871a4,
    ]),
    Fp::from_raw_unchecked([
        0xb416af000745fc20,
        0x8e563e9d1ea6d0f5,
        0x7c763e17763a0652,
        0x01458ef0159ebbef,
        0x8346fe421f96bb13,
        0x0d2d7b829ce324d2,
    ]),
    Fp::from_raw_unchecked([
        0x93096bb538d64615,
        0x6f2a2619951d823a,
        0x8f66b3ea59514fa4,
        0xf563e63704f7092f,
        0x724b136c4cf2d9fa,
        0x046959cfcfd0bf49,
    ]),
    Fp::from_raw_unchecked([
        0xea748d4b6e405346,
        0x91e9079c2c02d58f,
        0x41064965946d9b59,
        0xa06731f1d2bbe1ee,
        0x07f897e267a33f1b,
        0x1017290919210e5f,
    ]),
    Fp::from_raw_unchecked([
        0x872aa6c17d985097,
        0xeecc53161264562a,
        0x07afe37afff55002,
        0x54759078e5be6838,
        0xc4b92d15db8acca8,
        0x106d87d1b51d13b9,
    ]),
];

/// Coefficients of the 11-isogeny y map's denominator
const YDEN: [Fp; 16] = [
    Fp::from_raw_unchecked([
        0xeb6c359d47e52b1c,
        0x18ef5f8a10634d60,
        0xddfa71a0889d5b7e,
        0x723e71dcc5fc1323,
        0x52f45700b70d5c69,
        0x0a8b981ee47691f1,
    ]),
    Fp::from_raw_unchecked([
        0x616a3c4f5535b9fb,
        0x6f5f037395dbd911,
        0xf25f4cc5e35c65da,
        0x3e50dffea3c62658,
        0x6a33dca523560776,
        0x0fadeff77b6bfe3e,
    ]),
    Fp::from_raw_unchecked([
        0x2be9b66df470059c,
        0x24a2c159a3d36742,
        0x115dbe7ad10c2a37,
        0xb6634a652ee5884d,
        0x04fe8bb2b8d81af4,
        0x01c2a7a256fe9c41,
    ]),
    Fp::from_raw_unchecked([
        0xf27bf8ef3b75a386,
        0x898b367476c9073f,
        0x24482e6b8c2f4e5f,
        0xc8e0bbd6fe110806,
        0x59b0c17f7631448a,
        0x11037cd58b3dbfbd,
    ]),
    Fp::from_raw_unchecked([
        0x31c7912ea267eec6,
        0x1dbf6f1c5fcdb700,
        0xd30d4fe3ba86fdb1,
        0x3cae528fbee9a2a4,
        0xb1cce69b6aa9ad9a,
        0x044393bb632d94fb,
    ]),
    Fp::from_raw_unchecked([
        0xc66ef6efeeb5c7e8,
        0x9824c289dd72bb55,
        0x71b1a4d2f119981d,
        0x104fc1aafb0919cc,
        0x0e49df01d942a628,
        0x096c3a09773272d4,
    ]),
    Fp::from_raw_unchecked([
        0x9abc11eb5fadeff4,
        0x32dca50a885728f0,
        0xfb1fa3721569734c,
        0xc4b76271ea6506b3,
        0xd466a75599ce728e,
        0x0c81d4645f4cb6ed,
    ]),
    Fp::from_raw_unchecked([
        0x4199f10e5b8be45b,
        0xda64e495b1e87930,
        0xcb353efe9b33e4ff,
        0x9e9efb24aa6424c6,
        0xf08d33680a237465,
        0x0d3378023e4c7406,
    ]),
    Fp::from_raw_unchecked([
        0x7eb4ae92ec74d3a5,
        0xc341b4aa9fac3497,
        0x5be603899e907687,
        0x03bfd9cca75cbdeb,
        0x564c2935a96bfa93,
        0x0ef3c33371e2fdb5,
    ]),
    Fp::from_raw_unchecked([
        0x7ee91fd449f6ac2e,
        0xe5d5bd5cb9357a30,
        0x773a8ca5196b1380,
        0xd0fda172174ed023,
        0x6cb95e0fa776aead,
        0x0d22d5a40cec7cff,
    ]),
    Fp::from_raw_unchecked([
        0xf727e09285fd8519,
        0xdc9d55a83017897b,
        0x7549d8bd057894ae,
        0x178419613d90d8f8,
        0xfce95ebdeb5b490a,
        0x0467ffaef23fc49e,
    ]),
    Fp::from_raw_unchecked([
        0xc1769e6a7c385f1b,
        0x79bc930deac01c03,
        0x5461c75a23ede3b5,
        0x6e20829e5c230c45,
        0x828e0f1e772a53cd,
        0x116aefa749127bff,
    ]),
    Fp::from_raw_unchecked([
        0x101c10bf2744c10a,
        0xbbf18d053a6a3154,
        0xa0ecf39ef026f602,
        0xfc009d4996dc5153,
        0xb9000209d5bd08d3,
        0x189e5fe4470cd73c,
    ]),
    Fp::from_raw_unchecked([
        0x7ebd546ca1575ed2,
        0xe47d5a981d081b55,
        0x57b2b625b6d4ca21,
        0xb0a1ba04228520cc,
        0x98738983c2107ff3,
        0x13dddbc4799d81d6,
    ]),
    Fp::from_raw_unchecked([
        0x09319f2e39834935,
        0x039e952cbdb05c21,
        0x55ba77a9a2f76493,
        0xfd04e3dfc6086467,
        0xfb95832e7d78742e,
        0x0ef9c24eccaf5e0e,
    ]),
    Fp::from_raw_unchecked([
        0x760900000002fffd,
        0xebf4000bc40c0002,
        0x5f48985753c758ba,
        0x77ce585370525745,
        0x5c071a97a256ec6d,
        0x15f65ec3fa80e493,
    ]),
];

impl IsogenyMap for G1Projective {
    fn isogeny_map(&mut self) {
        self.eval_iso([&XNUM[..], &XDEN[..], &YNUM[..], &YDEN[..]]);
    }
}

pub struct G1Compressed([u8; 48]);

impl fmt::Debug for G1Compressed {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G1Compressed {
    fn default() -> Self {
        G1Compressed([0; 48])
    }
}

impl AsRef<[u8]> for G1Compressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G1Compressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

pub struct G1Uncompressed([u8; 96]);

impl fmt::Debug for G1Uncompressed {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G1Uncompressed {
    fn default() -> Self {
        G1Uncompressed([0; 96])
    }
}

impl AsRef<[u8]> for G1Uncompressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G1Uncompressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl Group for G1Projective {
    type Scalar = Scalar;

    fn random(mut rng: impl RngCore) -> Self {
        loop {
            let x = Fp::random(&mut rng);
            let flip_sign = rng.next_u32() % 2 != 0;

            // Obtain the corresponding y-coordinate given x as y = sqrt(x^3 + 4)
            let p = ((x.square() * x) + B).sqrt().map(|y| G1Affine {
                x,
                y: if flip_sign { -y } else { y },
                infinity: 0.into(),
            });

            if p.is_some().into() {
                let p = p.unwrap().to_curve().clear_cofactor();

                if bool::from(!p.is_identity()) {
                    return p;
                }
            }
        }
    }

    fn identity() -> Self {
        Self::identity()
    }

    fn generator() -> Self {
        Self::generator()
    }

    fn is_identity(&self) -> Choice {
        self.is_identity()
    }

    #[must_use]
    fn double(&self) -> Self {
        self.double()
    }
}

#[cfg(feature = "alloc")]
impl WnafGroup for G1Projective {
    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
        const RECOMMENDATIONS: [usize; 12] =
            [1, 3, 7, 20, 43, 120, 273, 563, 1630, 3128, 7933, 62569];

        let mut ret = 4;
        for r in &RECOMMENDATIONS {
            if num_scalars > *r {
                ret += 1;
            } else {
                break;
            }
        }

        ret
    }
}

impl PrimeGroup for G1Projective {}

impl Curve for G1Projective {
    type AffineRepr = G1Affine;

    fn batch_normalize(p: &[Self], q: &mut [Self::AffineRepr]) {
        Self::batch_normalize(p, q);
    }

    fn to_affine(&self) -> Self::AffineRepr {
        self.into()
    }
}

impl PrimeCurve for G1Projective {
    type Affine = G1Affine;
}

impl PrimeCurveAffine for G1Affine {
    type Scalar = Scalar;
    type Curve = G1Projective;

    fn identity() -> Self {
        Self::identity()
    }

    fn generator() -> Self {
        Self::generator()
    }

    fn is_identity(&self) -> Choice {
        self.is_identity()
    }

    fn to_curve(&self) -> Self::Curve {
        self.into()
    }
}

impl GroupEncoding for G1Projective {
    type Repr = G1Compressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        G1Affine::from_bytes(bytes).map(Self::from)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        G1Affine::from_bytes_unchecked(bytes).map(Self::from)
    }

    fn to_bytes(&self) -> Self::Repr {
        G1Affine::from(self).to_bytes()
    }
}

impl GroupEncoding for G1Affine {
    type Repr = G1Compressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_compressed(&bytes.0)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_compressed_unchecked(&bytes.0)
    }

    fn to_bytes(&self) -> Self::Repr {
        G1Compressed(self.to_compressed())
    }
}

impl UncompressedEncoding for G1Affine {
    type Uncompressed = G1Uncompressed;

    fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self> {
        Self::from_uncompressed(&bytes.0)
    }

    fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self> {
        Self::from_uncompressed_unchecked(&bytes.0)
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        G1Uncompressed(self.to_uncompressed())
    }
}

#[test]
fn test_is_on_curve() {
    assert!(bool::from(G1Affine::identity().is_on_curve()));
    assert!(bool::from(G1Affine::generator().is_on_curve()));
    assert!(bool::from(G1Projective::identity().is_on_curve()));
    assert!(bool::from(G1Projective::generator().is_on_curve()));

    let z = Fp::from_raw_unchecked([
        0xba7a_fa1f_9a6f_e250,
        0xfa0f_5b59_5eaf_e731,
        0x3bdc_4776_94c3_06e7,
        0x2149_be4b_3949_fa24,
        0x64aa_6e06_49b2_078c,
        0x12b1_08ac_3364_3c3e,
    ]);

    let gen = G1Affine::generator();
    let mut test = G1Projective {
        x: gen.x * z,
        y: gen.y * z,
        z,
    };

    assert!(bool::from(test.is_on_curve()));

    test.x = z;
    assert!(!bool::from(test.is_on_curve()));
}

#[test]
#[allow(clippy::eq_op)]
fn test_affine_point_equality() {
    let a = G1Affine::generator();
    let b = G1Affine::identity();

    assert!(a == a);
    assert!(b == b);
    assert!(a != b);
    assert!(b != a);
}

#[test]
#[allow(clippy::eq_op)]
fn test_projective_point_equality() {
    let a = G1Projective::generator();
    let b = G1Projective::identity();

    assert!(a == a);
    assert!(b == b);
    assert!(a != b);
    assert!(b != a);

    let z = Fp::from_raw_unchecked([
        0xba7a_fa1f_9a6f_e250,
        0xfa0f_5b59_5eaf_e731,
        0x3bdc_4776_94c3_06e7,
        0x2149_be4b_3949_fa24,
        0x64aa_6e06_49b2_078c,
        0x12b1_08ac_3364_3c3e,
    ]);

    let mut c = G1Projective {
        x: a.x * z,
        y: a.y * z,
        z,
    };
    assert!(bool::from(c.is_on_curve()));

    assert!(a == c);
    assert!(b != c);
    assert!(c == a);
    assert!(c != b);

    c.y = -c.y;
    assert!(bool::from(c.is_on_curve()));

    assert!(a != c);
    assert!(b != c);
    assert!(c != a);
    assert!(c != b);

    c.y = -c.y;
    c.x = z;
    assert!(!bool::from(c.is_on_curve()));
    assert!(a != b);
    assert!(a != c);
    assert!(b != c);
}

#[test]
fn test_conditionally_select_affine() {
    let a = G1Affine::generator();
    let b = G1Affine::identity();

    assert_eq!(G1Affine::conditional_select(&a, &b, Choice::from(0u8)), a);
    assert_eq!(G1Affine::conditional_select(&a, &b, Choice::from(1u8)), b);
}

#[test]
fn test_conditionally_select_projective() {
    let a = G1Projective::generator();
    let b = G1Projective::identity();

    assert_eq!(
        G1Projective::conditional_select(&a, &b, Choice::from(0u8)),
        a
    );
    assert_eq!(
        G1Projective::conditional_select(&a, &b, Choice::from(1u8)),
        b
    );
}

#[test]
fn test_projective_to_affine() {
    let a = G1Projective::generator();
    let b = G1Projective::identity();

    assert!(bool::from(G1Affine::from(a).is_on_curve()));
    assert!(!bool::from(G1Affine::from(a).is_identity()));
    assert!(bool::from(G1Affine::from(b).is_on_curve()));
    assert!(bool::from(G1Affine::from(b).is_identity()));

    let z = Fp::from_raw_unchecked([
        0xba7a_fa1f_9a6f_e250,
        0xfa0f_5b59_5eaf_e731,
        0x3bdc_4776_94c3_06e7,
        0x2149_be4b_3949_fa24,
        0x64aa_6e06_49b2_078c,
        0x12b1_08ac_3364_3c3e,
    ]);

    let c = G1Projective {
        x: a.x * z,
        y: a.y * z,
        z,
    };

    assert_eq!(G1Affine::from(c), G1Affine::generator());
}

#[test]
fn test_affine_to_projective() {
    let a = G1Affine::generator();
    let b = G1Affine::identity();

    assert!(bool::from(G1Projective::from(a).is_on_curve()));
    assert!(!bool::from(G1Projective::from(a).is_identity()));
    assert!(bool::from(G1Projective::from(b).is_on_curve()));
    assert!(bool::from(G1Projective::from(b).is_identity()));
}

#[test]
fn test_doubling() {
    {
        let tmp = G1Projective::identity().double();
        assert!(bool::from(tmp.is_identity()));
        assert!(bool::from(tmp.is_on_curve()));
    }
    {
        let tmp = G1Projective::generator().double();
        assert!(!bool::from(tmp.is_identity()));
        assert!(bool::from(tmp.is_on_curve()));

        assert_eq!(
            G1Affine::from(tmp),
            G1Affine {
                x: Fp::from_raw_unchecked([
                    0x53e9_78ce_58a9_ba3c,
                    0x3ea0_583c_4f3d_65f9,
                    0x4d20_bb47_f001_2960,
                    0xa54c_664a_e5b2_b5d9,
                    0x26b5_52a3_9d7e_b21f,
                    0x0008_895d_26e6_8785,
                ]),
                y: Fp::from_raw_unchecked([
                    0x7011_0b32_9829_3940,
                    0xda33_c539_3f1f_6afc,
                    0xb86e_dfd1_6a5a_a785,
                    0xaec6_d1c9_e7b1_c895,
                    0x25cf_c2b5_22d1_1720,
                    0x0636_1c83_f8d0_9b15,
                ]),
                infinity: Choice::from(0u8)
            }
        );
    }
}

#[test]
fn test_projective_addition() {
    {
        let a = G1Projective::identity();
        let b = G1Projective::identity();
        let c = a + b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
    {
        let a = G1Projective::identity();
        let mut b = G1Projective::generator();
        {
            let z = Fp::from_raw_unchecked([
                0xba7a_fa1f_9a6f_e250,
                0xfa0f_5b59_5eaf_e731,
                0x3bdc_4776_94c3_06e7,
                0x2149_be4b_3949_fa24,
                0x64aa_6e06_49b2_078c,
                0x12b1_08ac_3364_3c3e,
            ]);

            b = G1Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = a + b;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G1Projective::generator());
    }
    {
        let a = G1Projective::identity();
        let mut b = G1Projective::generator();
        {
            let z = Fp::from_raw_unchecked([
                0xba7a_fa1f_9a6f_e250,
                0xfa0f_5b59_5eaf_e731,
                0x3bdc_4776_94c3_06e7,
                0x2149_be4b_3949_fa24,
                0x64aa_6e06_49b2_078c,
                0x12b1_08ac_3364_3c3e,
            ]);

            b = G1Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = b + a;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G1Projective::generator());
    }
    {
        let a = G1Projective::generator().double().double(); // 4P
        let b = G1Projective::generator().double(); // 2P
        let c = a + b;

        let mut d = G1Projective::generator();
        for _ in 0..5 {
            d += G1Projective::generator();
        }
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(!bool::from(d.is_identity()));
        assert!(bool::from(d.is_on_curve()));
        assert_eq!(c, d);
    }

    // Degenerate case
    {
        let beta = Fp::from_raw_unchecked([
            0xcd03_c9e4_8671_f071,
            0x5dab_2246_1fcd_a5d2,
            0x5870_42af_d385_1b95,
            0x8eb6_0ebe_01ba_cb9e,
            0x03f9_7d6e_83d0_50d2,
            0x18f0_2065_5463_8741,
        ]);
        let beta = beta.square();
        let a = G1Projective::generator().double().double();
        let b = G1Projective {
            x: a.x * beta,
            y: -a.y,
            z: a.z,
        };
        assert!(bool::from(a.is_on_curve()));
        assert!(bool::from(b.is_on_curve()));

        let c = a + b;
        assert_eq!(
            G1Affine::from(c),
            G1Affine::from(G1Projective {
                x: Fp::from_raw_unchecked([
                    0x29e1_e987_ef68_f2d0,
                    0xc5f3_ec53_1db0_3233,
                    0xacd6_c4b6_ca19_730f,
                    0x18ad_9e82_7bc2_bab7,
                    0x46e3_b2c5_785c_c7a9,
                    0x07e5_71d4_2d22_ddd6,
                ]),
                y: Fp::from_raw_unchecked([
                    0x94d1_17a7_e5a5_39e7,
                    0x8e17_ef67_3d4b_5d22,
                    0x9d74_6aaf_508a_33ea,
                    0x8c6d_883d_2516_c9a2,
                    0x0bc3_b8d5_fb04_47f7,
                    0x07bf_a4c7_210f_4f44,
                ]),
                z: Fp::one()
            })
        );
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
}

#[test]
fn test_mixed_addition() {
    {
        let a = G1Affine::identity();
        let b = G1Projective::identity();
        let c = a + b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
    {
        let a = G1Affine::identity();
        let mut b = G1Projective::generator();
        {
            let z = Fp::from_raw_unchecked([
                0xba7a_fa1f_9a6f_e250,
                0xfa0f_5b59_5eaf_e731,
                0x3bdc_4776_94c3_06e7,
                0x2149_be4b_3949_fa24,
                0x64aa_6e06_49b2_078c,
                0x12b1_08ac_3364_3c3e,
            ]);

            b = G1Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = a + b;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G1Projective::generator());
    }
    {
        let a = G1Affine::identity();
        let mut b = G1Projective::generator();
        {
            let z = Fp::from_raw_unchecked([
                0xba7a_fa1f_9a6f_e250,
                0xfa0f_5b59_5eaf_e731,
                0x3bdc_4776_94c3_06e7,
                0x2149_be4b_3949_fa24,
                0x64aa_6e06_49b2_078c,
                0x12b1_08ac_3364_3c3e,
            ]);

            b = G1Projective {
                x: b.x * z,
                y: b.y * z,
                z,
            };
        }
        let c = b + a;
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(c == G1Projective::generator());
    }
    {
        let a = G1Projective::generator().double().double(); // 4P
        let b = G1Projective::generator().double(); // 2P
        let c = a + b;

        let mut d = G1Projective::generator();
        for _ in 0..5 {
            d += G1Affine::generator();
        }
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(!bool::from(d.is_identity()));
        assert!(bool::from(d.is_on_curve()));
        assert_eq!(c, d);
    }

    // Degenerate case
    {
        let beta = Fp::from_raw_unchecked([
            0xcd03_c9e4_8671_f071,
            0x5dab_2246_1fcd_a5d2,
            0x5870_42af_d385_1b95,
            0x8eb6_0ebe_01ba_cb9e,
            0x03f9_7d6e_83d0_50d2,
            0x18f0_2065_5463_8741,
        ]);
        let beta = beta.square();
        let a = G1Projective::generator().double().double();
        let b = G1Projective {
            x: a.x * beta,
            y: -a.y,
            z: a.z,
        };
        let a = G1Affine::from(a);
        assert!(bool::from(a.is_on_curve()));
        assert!(bool::from(b.is_on_curve()));

        let c = a + b;
        assert_eq!(
            G1Affine::from(c),
            G1Affine::from(G1Projective {
                x: Fp::from_raw_unchecked([
                    0x29e1_e987_ef68_f2d0,
                    0xc5f3_ec53_1db0_3233,
                    0xacd6_c4b6_ca19_730f,
                    0x18ad_9e82_7bc2_bab7,
                    0x46e3_b2c5_785c_c7a9,
                    0x07e5_71d4_2d22_ddd6,
                ]),
                y: Fp::from_raw_unchecked([
                    0x94d1_17a7_e5a5_39e7,
                    0x8e17_ef67_3d4b_5d22,
                    0x9d74_6aaf_508a_33ea,
                    0x8c6d_883d_2516_c9a2,
                    0x0bc3_b8d5_fb04_47f7,
                    0x07bf_a4c7_210f_4f44,
                ]),
                z: Fp::one()
            })
        );
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
    }
}

#[test]
#[allow(clippy::eq_op)]
fn test_projective_negation_and_subtraction() {
    let a = G1Projective::generator().double();
    assert_eq!(a + (-a), G1Projective::identity());
    assert_eq!(a + (-a), a - a);
}

#[test]
fn test_affine_negation_and_subtraction() {
    let a = G1Affine::generator();
    assert_eq!(G1Projective::from(a) + (-a), G1Projective::identity());
    assert_eq!(G1Projective::from(a) + (-a), G1Projective::from(a) - a);
}

#[test]
fn test_projective_scalar_multiplication() {
    let g = G1Projective::generator();
    let a = Scalar::from_raw([
        0x2b56_8297_a56d_a71c,
        0xd8c3_9ecb_0ef3_75d1,
        0x435c_38da_67bf_bf96,
        0x8088_a050_26b6_59b2,
    ]);
    let b = Scalar::from_raw([
        0x785f_dd9b_26ef_8b85,
        0xc997_f258_3769_5c18,
        0x4c8d_bc39_e7b7_56c1,
        0x70d9_b6cc_6d87_df20,
    ]);
    let c = a * b;

    assert_eq!((g * a) * b, g * c);
}

#[test]
fn test_affine_scalar_multiplication() {
    let g = G1Affine::generator();
    let a = Scalar::from_raw([
        0x2b56_8297_a56d_a71c,
        0xd8c3_9ecb_0ef3_75d1,
        0x435c_38da_67bf_bf96,
        0x8088_a050_26b6_59b2,
    ]);
    let b = Scalar::from_raw([
        0x785f_dd9b_26ef_8b85,
        0xc997_f258_3769_5c18,
        0x4c8d_bc39_e7b7_56c1,
        0x70d9_b6cc_6d87_df20,
    ]);
    let c = a * b;

    assert_eq!(G1Affine::from(g * a) * b, g * c);
}

#[test]
fn test_is_torsion_free() {
    let a = G1Affine {
        x: Fp::from_raw_unchecked([
            0x0aba_f895_b97e_43c8,
            0xba4c_6432_eb9b_61b0,
            0x1250_6f52_adfe_307f,
            0x7502_8c34_3933_6b72,
            0x8474_4f05_b8e9_bd71,
            0x113d_554f_b095_54f7,
        ]),
        y: Fp::from_raw_unchecked([
            0x73e9_0e88_f5cf_01c0,
            0x3700_7b65_dd31_97e2,
            0x5cf9_a199_2f0d_7c78,
            0x4f83_c10b_9eb3_330d,
            0xf6a6_3f6f_07f6_0961,
            0x0c53_b5b9_7e63_4df3,
        ]),
        infinity: Choice::from(0u8),
    };
    assert!(!bool::from(a.is_torsion_free()));

    assert!(bool::from(G1Affine::identity().is_torsion_free()));
    assert!(bool::from(G1Affine::generator().is_torsion_free()));
}

#[test]
fn test_mul_by_x() {
    // multiplying by `x` a point in G1 is the same as multiplying by
    // the equivalent scalar.
    let generator = G1Projective::generator();
    let x = if crate::BLS_X_IS_NEGATIVE {
        -Scalar::from(crate::BLS_X)
    } else {
        Scalar::from(crate::BLS_X)
    };
    assert_eq!(generator.mul_by_x(), generator * x);

    let point = G1Projective::generator() * Scalar::from(42);
    assert_eq!(point.mul_by_x(), point * x);
}

#[test]
fn test_clear_cofactor() {
    // the generator (and the identity) are always on the curve,
    // even after clearing the cofactor
    let generator = G1Projective::generator();
    assert!(bool::from(generator.clear_cofactor().is_on_curve()));
    let id = G1Projective::identity();
    assert!(bool::from(id.clear_cofactor().is_on_curve()));

    let z = Fp::from_raw_unchecked([
        0x3d2d1c670671394e,
        0x0ee3a800a2f7c1ca,
        0x270f4f21da2e5050,
        0xe02840a53f1be768,
        0x55debeb597512690,
        0x08bd25353dc8f791,
    ]);

    let point = G1Projective {
        x: Fp::from_raw_unchecked([
            0x48af5ff540c817f0,
            0xd73893acaf379d5a,
            0xe6c43584e18e023c,
            0x1eda39c30f188b3e,
            0xf618c6d3ccc0f8d8,
            0x0073542cd671e16c,
        ]) * z,
        y: Fp::from_raw_unchecked([
            0x57bf8be79461d0ba,
            0xfc61459cee3547c3,
            0x0d23567df1ef147b,
            0x0ee187bcce1d9b64,
            0xb0c8cfbe9dc8fdc1,
            0x1328661767ef368b,
        ]),
        z: z.square() * z,
    };

    assert!(bool::from(point.is_on_curve()));
    assert!(!bool::from(G1Affine::from(point).is_torsion_free()));
    let cleared_point = point.clear_cofactor();
    assert!(bool::from(cleared_point.is_on_curve()));
    assert!(bool::from(G1Affine::from(cleared_point).is_torsion_free()));

    // in BLS12-381 the cofactor in G1 can be
    // cleared multiplying by (1-x)
    let h_eff = Scalar::from(1) + Scalar::from(crate::BLS_X);
    assert_eq!(point.clear_cofactor(), point * h_eff);
}

#[test]
fn test_batch_normalize() {
    let a = G1Projective::generator().double();
    let b = a.double();
    let c = b.double();

    for a_identity in (0..1).map(|n| n == 1) {
        for b_identity in (0..1).map(|n| n == 1) {
            for c_identity in (0..1).map(|n| n == 1) {
                let mut v = [a, b, c];
                if a_identity {
                    v[0] = G1Projective::identity()
                }
                if b_identity {
                    v[1] = G1Projective::identity()
                }
                if c_identity {
                    v[2] = G1Projective::identity()
                }

                let mut t = [
                    G1Affine::identity(),
                    G1Affine::identity(),
                    G1Affine::identity(),
                ];
                let expected = [
                    G1Affine::from(v[0]),
                    G1Affine::from(v[1]),
                    G1Affine::from(v[2]),
                ];

                G1Projective::batch_normalize(&v[..], &mut t[..]);

                assert_eq!(&t[..], &expected[..]);
            }
        }
    }
}

#[test]
fn test_projective_iso11_zero() {
    let zero = Fp::zero();
    let mut pt = G1Projective {
        x: Fp::zero(),
        y: Fp::zero(),
        z: Fp::zero(),
    };
    pt.isogeny_map();
    assert_eq!(pt.x, zero);
    assert_eq!(pt.y, zero);
    assert_eq!(pt.z, zero);
}

#[test]
fn test_projective_iso11_one() {
    let mut pt = G1Projective {
        x: Fp::one(),
        y: Fp::one(),
        z: Fp::one(),
    };
    pt.isogeny_map();
    assert_eq!(
        pt.x,
        Fp::from_raw_unchecked([
            0xc02cdce8a2f88d1b,
            0x5c5d398a0fffbaae,
            0x8183696afb034ce6,
            0xc711c18a9199f2fc,
            0xbaa27175064212a5,
            0x133080bfaef68b8d
        ])
    );
    assert_eq!(
        pt.y,
        Fp::from_raw_unchecked([
            0x2a9cfaaa906a40ec,
            0xd93a736a56402602,
            0xee9e99207bbc4d29,
            0xdb0aa86255189ad5,
            0xe0cf89b1e6badd9d,
            0x1952f3da38fd3f9
        ])
    );
    assert_eq!(
        pt.z,
        Fp::from_raw_unchecked([
            0x964515fe48c609db,
            0x6735ed08a11adde5,
            0x28a2b43715a40f,
            0x4d8ec57a748671de,
            0xf67c258133b77d96,
            0x1723f4caff98e8d2
        ])
    );
}

#[test]
fn test_projective_iso11_fixed() {
    let xi = Fp::from_raw_unchecked([
        0x62653ef95550b55a,
        0x238de03a1308654,
        0x19623cd2898b0e9b,
        0x484bf2559df8de8a,
        0xe2e93595da9983f5,
        0x39b5d6b0479e04,
    ]);
    let yi = Fp::from_raw_unchecked([
        0xe9fe6726d2d392bb,
        0xee7122a3b30bd837,
        0xcb2283c872cab64b,
        0xf6bd444fc828fb36,
        0x8061113f1da387bd,
        0x88bd49568bdad70,
    ]);
    let zi = Fp::from_raw_unchecked([
        0xc81bb8f8c311a54f,
        0x37d939e4f58ae80b,
        0xf78fc1246067c90b,
        0x1abffdad4be45003,
        0x97bc93b1d4fd1bdf,
        0xbebe6e3794aa5f8,
    ]);

    let mut pt = G1Projective {
        x: xi,
        y: yi,
        z: zi,
    };
    pt.isogeny_map();
    assert_eq!(
        pt.x,
        Fp::from_raw_unchecked([
            0x98fd5875e4a58e20,
            0xf00708855ba24728,
            0xc54b92477ec91cf8,
            0x62147c9f75dc4da0,
            0x784ea4acc57cc997,
            0x13f1dc06973b63c2
        ])
    );
    assert_eq!(
        pt.y,
        Fp::from_raw_unchecked([
            0x641ce4fa319fb527,
            0xb25350f92bb911f0,
            0x60f575e35412341a,
            0x7cc7cdd785a14540,
            0x912d881277ec9fbf,
            0x3dd9bf1f358a671,
        ])
    );
    assert_eq!(
        pt.z,
        Fp::from_raw_unchecked([
            0xc52c49bfd570446c,
            0x2673d927b164af30,
            0x2ab772e0665915f1,
            0xd5dc90b3e7644d3d,
            0xccce39bd5f3ee359,
            0xeaa09526b938f64,
        ])
    );
}
