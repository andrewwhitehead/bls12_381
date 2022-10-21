use crate::fp12::Fp12;
use crate::fp2::Fp2;
use crate::g2::fp2_mul_by_3b;
use crate::gt::Gt;
use crate::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar, BLS_X, BLS_X_IS_NEGATIVE};

use core::ops;

use pairing::{Engine, PairingCurveAffine};
use subtle::{Choice, ConditionallySelectable};

#[cfg(feature = "alloc")]
use alloc::{boxed::Box, vec::Vec};
#[cfg(feature = "alloc")]
use pairing::MultiMillerLoop;

/// Represents results of a Miller loop, one of the most expensive portions
/// of the pairing function. `MillerLoopResult`s cannot be compared with each
/// other until `.final_exponentiation()` is called, which is also expensive.
#[cfg_attr(docsrs, doc(cfg(feature = "pairings")))]
#[derive(Copy, Clone, Debug)]
pub struct MillerLoopResult(pub(crate) Fp12);

impl Default for MillerLoopResult {
    fn default() -> Self {
        MillerLoopResult(Fp12::one())
    }
}

#[cfg(feature = "zeroize")]
impl zeroize::DefaultIsZeroes for MillerLoopResult {}

impl ConditionallySelectable for MillerLoopResult {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        MillerLoopResult(Fp12::conditional_select(&a.0, &b.0, choice))
    }
}

impl MillerLoopResult {
    /// This performs a "final exponentiation" routine to convert the result
    /// of a Miller loop into an element of `Gt` with help of efficient squaring
    /// operation in the so-called `cyclotomic subgroup` of `Fq6` so that
    /// it can be compared with other elements of `Gt`.
    pub fn final_exponentiation(&self) -> Gt {
        Gt(self.0.final_exponentiation())
    }
}

impl<'a, 'b> ops::Add<&'b MillerLoopResult> for &'a MillerLoopResult {
    type Output = MillerLoopResult;

    #[inline]
    fn add(self, rhs: &'b MillerLoopResult) -> MillerLoopResult {
        MillerLoopResult(self.0 * rhs.0)
    }
}

impl_add_binop_specify_output!(MillerLoopResult, MillerLoopResult, MillerLoopResult);

impl ops::AddAssign<MillerLoopResult> for MillerLoopResult {
    #[inline]
    fn add_assign(&mut self, rhs: MillerLoopResult) {
        *self = *self + rhs;
    }
}

impl<'b> ops::AddAssign<&'b MillerLoopResult> for MillerLoopResult {
    #[inline]
    fn add_assign(&mut self, rhs: &'b MillerLoopResult) {
        *self = *self + rhs;
    }
}

#[cfg(feature = "alloc")]
#[cfg_attr(docsrs, doc(cfg(all(feature = "pairings", feature = "alloc"))))]
#[derive(Clone, Debug)]
/// This structure contains cached computations pertaining to a $\mathbb{G}_2$
/// element as part of the pairing function (specifically, the Miller loop) and
/// so should be computed whenever a $\mathbb{G}_2$ element is being used in
/// multiple pairings or is otherwise known in advance. This should be used in
/// conjunction with the [`multi_miller_loop`](crate::multi_miller_loop)
/// function provided by this crate.
///
/// Requires the `alloc` and `pairing` crate features to be enabled.
pub struct G2Prepared {
    infinity: Choice,
    coeffs: Box<[(Fp2, Fp2, Fp2)]>,
}

#[cfg(feature = "alloc")]
impl From<G2Affine> for G2Prepared {
    fn from(q: G2Affine) -> G2Prepared {
        struct Adder {
            cur: G2Projective,
            base: G2Affine,
            coeffs: Vec<(Fp2, Fp2, Fp2)>,
        }

        impl MillerLoopDriver for Adder {
            type Output = ();

            #[inline]
            fn doubling_step(&mut self, _: Self::Output) -> Self::Output {
                let coeffs = doubling_step(&mut self.cur);
                self.coeffs.push(coeffs);
            }
            #[inline]
            fn addition_step(&mut self, _: Self::Output) -> Self::Output {
                let coeffs = addition_step(&mut self.cur, &self.base);
                self.coeffs.push(coeffs);
            }
            #[inline(always)]
            fn square_output(_: Self::Output) -> Self::Output {}
            #[inline(always)]
            fn conjugate(_: Self::Output) -> Self::Output {}
            #[inline(always)]
            fn one() -> Self::Output {}
        }

        let is_identity = q.is_identity();
        let q = G2Affine::conditional_select(&q, &G2Affine::generator(), is_identity);

        let mut adder = Adder {
            cur: G2Projective::from(q),
            base: q,
            coeffs: Vec::new(),
        };
        adder.coeffs.reserve_exact(68);

        miller_loop(&mut adder);

        assert_eq!(adder.coeffs.len(), 68);

        G2Prepared {
            infinity: is_identity,
            coeffs: adder.coeffs.into_boxed_slice(),
        }
    }
}

#[cfg(feature = "alloc")]
#[cfg_attr(docsrs, doc(cfg(all(feature = "pairings", feature = "alloc"))))]
/// Computes $$\sum_{i=1}^n \textbf{ML}(a_i, b_i)$$ given a series of terms
/// $$(a_1, b_1), (a_2, b_2), ..., (a_n, b_n).$$
///
/// Requires the `alloc` and `pairing` crate features to be enabled.
pub fn multi_miller_loop(terms: &[(&G1Affine, &G2Prepared)]) -> MillerLoopResult {
    struct Adder<'a, 'b, 'c> {
        terms: &'c [(&'a G1Affine, &'b G2Prepared)],
        index: usize,
    }

    impl<'a, 'b, 'c> MillerLoopDriver for Adder<'a, 'b, 'c> {
        type Output = Fp12;

        #[inline]
        fn doubling_step(&mut self, mut f: Self::Output) -> Self::Output {
            let index = self.index;
            for &(p, r) in self.terms {
                let either_identity = p.is_identity() | r.infinity;

                let new_f = ell(f, &r.coeffs[index], p);
                f = Fp12::conditional_select(&new_f, &f, either_identity);
            }
            self.index += 1;

            f
        }
        #[inline]
        fn addition_step(&mut self, mut f: Self::Output) -> Self::Output {
            let index = self.index;
            for &(p, r) in self.terms {
                let either_identity = p.is_identity() | r.infinity;

                let new_f = ell(f, &r.coeffs[index], p);
                f = Fp12::conditional_select(&new_f, &f, either_identity);
            }
            self.index += 1;

            f
        }
        #[inline]
        fn square_output(f: Self::Output) -> Self::Output {
            f.square()
        }
        #[inline]
        fn conjugate(f: Self::Output) -> Self::Output {
            f.conjugate()
        }
        #[inline(always)]
        fn one() -> Self::Output {
            Fp12::one()
        }
    }

    let mut adder = Adder { terms, index: 0 };

    let tmp = miller_loop(&mut adder);

    MillerLoopResult(tmp)
}

/// Invoke the pairing function without the use of precomputation and other optimizations.
#[cfg_attr(docsrs, doc(cfg(feature = "pairings")))]
pub fn pairing(p: &G1Affine, q: &G2Affine) -> Gt {
    struct Adder {
        cur: G2Projective,
        base: G2Affine,
        p: G1Affine,
    }

    impl MillerLoopDriver for Adder {
        type Output = Fp12;

        #[inline]
        fn doubling_step(&mut self, f: Self::Output) -> Self::Output {
            let coeffs = doubling_step(&mut self.cur);
            ell(f, &coeffs, &self.p)
        }
        #[inline]
        fn addition_step(&mut self, f: Self::Output) -> Self::Output {
            let coeffs = addition_step(&mut self.cur, &self.base);
            ell(f, &coeffs, &self.p)
        }
        #[inline]
        fn square_output(f: Self::Output) -> Self::Output {
            f.square()
        }
        #[inline]
        fn conjugate(f: Self::Output) -> Self::Output {
            f.conjugate()
        }
        #[inline(always)]
        fn one() -> Self::Output {
            Fp12::one()
        }
    }

    let either_identity = p.is_identity() | q.is_identity();
    let p = G1Affine::conditional_select(p, &G1Affine::generator(), either_identity);
    let q = G2Affine::conditional_select(q, &G2Affine::generator(), either_identity);

    let mut adder = Adder {
        cur: G2Projective::from(q),
        base: q,
        p,
    };

    let tmp = miller_loop(&mut adder);
    let tmp = MillerLoopResult(Fp12::conditional_select(
        &tmp,
        &Fp12::one(),
        either_identity,
    ));
    tmp.final_exponentiation()
}

/// Invoke the pairing function without the use of precomputation and other optimizations.
#[cfg_attr(docsrs, doc(cfg(feature = "pairings")))]
pub fn multi_pairing<const T: usize>(terms: &[(&G1Affine, &G2Affine); T]) -> Gt {
    struct Adder<'a, const U: usize> {
        terms: &'a [(&'a G1Affine, &'a G2Affine)],
        cur: [G2Projective; U],
    }

    impl<'a, const U: usize> MillerLoopDriver for Adder<'a, U> {
        type Output = Fp12;

        #[inline]
        fn doubling_step(&mut self, mut f: Self::Output) -> Self::Output {
            for (idx, &(p, r)) in self.terms.iter().enumerate() {
                let either_identity = p.is_identity() | r.is_identity();

                let coeffs = doubling_step(&mut self.cur[idx]);
                let new_f = ell(f, &coeffs, p);
                f = Fp12::conditional_select(&new_f, &f, either_identity);
            }
            f
        }
        #[inline]
        fn addition_step(&mut self, mut f: Self::Output) -> Self::Output {
            for (idx, &(p, r)) in self.terms.iter().enumerate() {
                let either_identity = p.is_identity() | r.is_identity();

                let coeffs = addition_step(&mut self.cur[idx], r);
                let new_f = ell(f, &coeffs, p);
                f = Fp12::conditional_select(&new_f, &f, either_identity);
            }
            f
        }
        #[inline]
        fn square_output(f: Self::Output) -> Self::Output {
            f.square()
        }
        #[inline]
        fn conjugate(f: Self::Output) -> Self::Output {
            f.conjugate()
        }
        #[inline(always)]
        fn one() -> Self::Output {
            Fp12::one()
        }
    }

    let mut adder = Adder {
        terms,
        cur: [G2Projective::identity(); T],
    };
    for idx in 0..T {
        adder.cur[idx] = G2Projective::from(terms[idx].1);
    }

    let tmp = miller_loop(&mut adder);
    MillerLoopResult(tmp).final_exponentiation()
}

trait MillerLoopDriver {
    type Output;

    fn doubling_step(&mut self, f: Self::Output) -> Self::Output;
    fn addition_step(&mut self, f: Self::Output) -> Self::Output;
    fn square_output(f: Self::Output) -> Self::Output;
    fn conjugate(f: Self::Output) -> Self::Output;
    fn one() -> Self::Output;
}

/// This is a "generic" implementation of the Miller loop to avoid duplicating code
/// structure elsewhere; instead, we'll write concrete instantiations of
/// `MillerLoopDriver` for whatever purposes we need (such as caching modes).
#[inline]
fn miller_loop<D: MillerLoopDriver>(driver: &mut D) -> D::Output {
    let mut acc = D::one();
    let mut i = 1u64 << 62;
    loop {
        acc = driver.doubling_step(acc);
        if BLS_X & i != 0 {
            acc = driver.addition_step(acc);
        }
        i >>= 1;
        if i == 0 {
            break;
        }
        acc = D::square_output(acc);
    }
    if BLS_X_IS_NEGATIVE {
        D::conjugate(acc)
    } else {
        acc
    }
}

#[inline]
fn ell(f: Fp12, coeffs: &(Fp2, Fp2, Fp2), p: &G1Affine) -> Fp12 {
    let (mut c4, mut c1, ref c0) = coeffs;

    c1.c0 *= p.x;
    c1.c1 *= p.x;

    c4.c0 *= p.y;
    c4.c1 *= p.y;

    f.mul_by_014(c0, &c1, &c4)
}

fn doubling_step(r: &mut G2Projective) -> (Fp2, Fp2, Fp2) {
    // Adaptation of Section 5 of https://eprint.iacr.org/2009/615.pdf

    // A = X1^2, B = Y1^2, C = Z1^2
    let a = r.x.square();
    let b = r.y.square();
    let c = r.z.square();

    // D = 3b'C, E = (X1 + Y1)^2 - A - B
    let d = fp2_mul_by_3b(c);
    let e = (r.x * r.y).double();

    // F = (Y1 + Z1)^2 - B - C, G = 3D
    let f = (r.y * r.z).double();
    let g = d.double() + d;

    // X3 = E * (B - G), Y3 = (B + G)^2 - 12D^2, Z3 = 4BF
    r.x = e * (b - g);
    let t0 = d.double().square();
    r.y = (b + g).square() - t0.double() - t0;
    r.z = b * f.double().double();

    // l00 = D - B, l11 = -F, l01 = 3A
    let l00 = d - b;
    let l11 = -f;
    let l01 = a.double() + a;

    (l11, l01, l00)
}

fn addition_step(r: &mut G2Projective, q: &G2Affine) -> (Fp2, Fp2, Fp2) {
    // Adaptation of Section 4.3, https://eprint.iacr.org/2013/722.pdf

    // A = Y1 - Z1 * Y2
    let a = r.y - r.z * q.y;

    // B = X1 - Z1 * X2
    let b = r.x - r.z * q.x;

    // C = A^2, D = B^2, E = B^3, F = E + Z1 * C
    let c = a.square();
    let d = b.square();
    let e = b * d;
    let f = e + r.z * c;

    // G = X1 * D, H = F - 2 * G
    let g = r.x * d;
    let h = f - g.double();

    // X3 = B * H, Y3 = A * (G - H) - Y1 * E, Z3 = Z1 * E
    r.x = b * h;
    r.y = Fp2::lin_comb(&a, &(g - h), &-r.y, &e);
    r.z = r.z * e;

    // l00 = A * X2 - B * Y2, l11 = B, l01 = -A
    let l00 = Fp2::lin_comb(&a, &q.x, &-b, &q.y);
    let l01 = -a;
    let l11 = b;

    (l11, l01, l00)
}

impl PairingCurveAffine for G1Affine {
    type Pair = G2Affine;
    type PairingResult = Gt;

    fn pairing_with(&self, other: &Self::Pair) -> Self::PairingResult {
        pairing(self, other)
    }
}

impl PairingCurveAffine for G2Affine {
    type Pair = G1Affine;
    type PairingResult = Gt;

    fn pairing_with(&self, other: &Self::Pair) -> Self::PairingResult {
        pairing(other, self)
    }
}

/// A [`pairing::Engine`] for BLS12-381 pairing operations.
#[cfg_attr(docsrs, doc(cfg(feature = "pairings")))]
#[derive(Copy, Clone, Debug)]
pub struct Bls12;

impl Engine for Bls12 {
    type Fr = Scalar;
    type G1 = G1Projective;
    type G1Affine = G1Affine;
    type G2 = G2Projective;
    type G2Affine = G2Affine;
    type Gt = Gt;

    fn pairing(p: &Self::G1Affine, q: &Self::G2Affine) -> Self::Gt {
        pairing(p, q)
    }
}

impl pairing::MillerLoopResult for MillerLoopResult {
    type Gt = Gt;

    fn final_exponentiation(&self) -> Self::Gt {
        self.final_exponentiation()
    }
}

#[cfg(feature = "alloc")]
impl MultiMillerLoop for Bls12 {
    type G2Prepared = G2Prepared;
    type Result = MillerLoopResult;

    fn multi_miller_loop(terms: &[(&Self::G1Affine, &Self::G2Prepared)]) -> Self::Result {
        multi_miller_loop(terms)
    }
}

#[test]
fn test_gt_generator() {
    assert_eq!(
        Gt::generator(),
        pairing(&G1Affine::generator(), &G2Affine::generator())
    );
}

#[test]
fn test_bilinearity() {
    use crate::Scalar;

    let a = Scalar::from_raw([1, 2, 3, 4]).invert().unwrap().square();
    let b = Scalar::from_raw([5, 6, 7, 8]).invert().unwrap().square();
    let c = a * b;

    let g = G1Affine::from(G1Affine::generator() * a);
    let h = G2Affine::from(G2Affine::generator() * b);
    let p = pairing(&g, &h);

    assert!(p != Gt::identity());

    let expected = G1Affine::from(G1Affine::generator() * c);

    assert_eq!(p, pairing(&expected, &G2Affine::generator()));
    assert_eq!(
        p,
        pairing(&G1Affine::generator(), &G2Affine::generator()) * c
    );
}

#[test]
fn test_unitary() {
    let g = G1Affine::generator();
    let h = G2Affine::generator();
    let p = -pairing(&g, &h);
    let q = pairing(&g, &-h);
    let r = pairing(&-g, &h);

    assert_eq!(p, q);
    assert_eq!(q, r);
}

#[cfg(feature = "alloc")]
#[test]
fn test_multi_miller_loop() {
    let a1 = G1Affine::generator();
    let b1 = G2Affine::generator();

    let a2 = G1Affine::from(
        G1Affine::generator() * Scalar::from_raw([1, 2, 3, 4]).invert().unwrap().square(),
    );
    let b2 = G2Affine::from(
        G2Affine::generator() * Scalar::from_raw([4, 2, 2, 4]).invert().unwrap().square(),
    );

    let a3 = G1Affine::identity();
    let b3 = G2Affine::from(
        G2Affine::generator() * Scalar::from_raw([9, 2, 2, 4]).invert().unwrap().square(),
    );

    let a4 = G1Affine::from(
        G1Affine::generator() * Scalar::from_raw([5, 5, 5, 5]).invert().unwrap().square(),
    );
    let b4 = G2Affine::identity();

    let a5 = G1Affine::from(
        G1Affine::generator() * Scalar::from_raw([323, 32, 3, 1]).invert().unwrap().square(),
    );
    let b5 = G2Affine::from(
        G2Affine::generator() * Scalar::from_raw([4, 2, 2, 9099]).invert().unwrap().square(),
    );

    let b1_prepared = G2Prepared::from(b1);
    let b2_prepared = G2Prepared::from(b2);
    let b3_prepared = G2Prepared::from(b3);
    let b4_prepared = G2Prepared::from(b4);
    let b5_prepared = G2Prepared::from(b5);

    let expected = pairing(&a1, &b1)
        + pairing(&a2, &b2)
        + pairing(&a3, &b3)
        + pairing(&a4, &b4)
        + pairing(&a5, &b5);

    let test = multi_miller_loop(&[
        (&a1, &b1_prepared),
        (&a2, &b2_prepared),
        (&a3, &b3_prepared),
        (&a4, &b4_prepared),
        (&a5, &b5_prepared),
    ])
    .final_exponentiation();

    assert_eq!(expected, test);

    let test2 = multi_pairing(&[(&a1, &b1), (&a2, &b2), (&a3, &b3), (&a4, &b4), (&a5, &b5)]);
    assert_eq!(expected, test2);
}

#[test]
fn test_miller_loop_result_default() {
    assert_eq!(
        MillerLoopResult::default().final_exponentiation(),
        Gt::identity(),
    );
}

#[cfg(feature = "zeroize")]
#[test]
fn test_miller_loop_result_zeroize() {
    use zeroize::Zeroize;

    let mut m = multi_miller_loop(&[
        (&G1Affine::generator(), &G2Affine::generator().into()),
        (&-G1Affine::generator(), &G2Affine::generator().into()),
    ]);
    m.zeroize();
    assert_eq!(m.0, MillerLoopResult::default().0);
}

#[test]
fn tricking_miller_loop_result() {
    assert_eq!(
        multi_miller_loop(&[(&G1Affine::identity(), &G2Affine::generator().into())]).0,
        Fp12::one()
    );
    assert_eq!(
        multi_miller_loop(&[(&G1Affine::generator(), &G2Affine::identity().into())]).0,
        Fp12::one()
    );
    assert_ne!(
        multi_miller_loop(&[
            (&G1Affine::generator(), &G2Affine::generator().into()),
            (&-G1Affine::generator(), &G2Affine::generator().into())
        ])
        .0,
        Fp12::one()
    );
    assert_eq!(
        multi_miller_loop(&[
            (&G1Affine::generator(), &G2Affine::generator().into()),
            (&-G1Affine::generator(), &G2Affine::generator().into())
        ])
        .final_exponentiation(),
        Gt::identity()
    );
}

#[test]
fn tricking_multi_pairing_result() {
    assert_eq!(
        multi_miller_loop(&[(&G1Affine::identity(), &G2Affine::generator().into())]).0,
        Fp12::one()
    );
    assert_eq!(
        multi_miller_loop(&[(&G1Affine::generator(), &G2Affine::identity().into())]).0,
        Fp12::one()
    );
    assert_ne!(
        multi_miller_loop(&[
            (&G1Affine::generator(), &G2Affine::generator().into()),
            (&-G1Affine::generator(), &G2Affine::generator().into())
        ])
        .0,
        Fp12::one()
    );
    assert_eq!(
        multi_miller_loop(&[
            (&G1Affine::generator(), &G2Affine::generator().into()),
            (&-G1Affine::generator(), &G2Affine::generator().into())
        ])
        .final_exponentiation(),
        Gt::identity()
    );
}
