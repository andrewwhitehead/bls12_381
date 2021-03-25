/// Compute a + b + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
pub const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a + (b * c) + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

macro_rules! impl_add_binop_specify_output {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl<'b> Add<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn add(self, rhs: &'b $rhs) -> $output {
                &self + rhs
            }
        }

        impl<'a> Add<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn add(self, rhs: $rhs) -> $output {
                self + &rhs
            }
        }

        impl Add<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn add(self, rhs: $rhs) -> $output {
                &self + &rhs
            }
        }
    };
}

macro_rules! impl_sub_binop_specify_output {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl<'b> Sub<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: &'b $rhs) -> $output {
                &self - rhs
            }
        }

        impl<'a> Sub<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: $rhs) -> $output {
                self - &rhs
            }
        }

        impl Sub<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: $rhs) -> $output {
                &self - &rhs
            }
        }
    };
}

macro_rules! impl_binops_additive_specify_output {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl_add_binop_specify_output!($lhs, $rhs, $output);
        impl_sub_binop_specify_output!($lhs, $rhs, $output);
    };
}

macro_rules! impl_binops_multiplicative_mixed {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl<'b> Mul<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: &'b $rhs) -> $output {
                &self * rhs
            }
        }

        impl<'a> Mul<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: $rhs) -> $output {
                self * &rhs
            }
        }

        impl Mul<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: $rhs) -> $output {
                &self * &rhs
            }
        }
    };
}

macro_rules! impl_binops_additive {
    ($lhs:ident, $rhs:ident) => {
        impl_binops_additive_specify_output!($lhs, $rhs, $lhs);

        impl SubAssign<$rhs> for $lhs {
            #[inline]
            fn sub_assign(&mut self, rhs: $rhs) {
                *self = &*self - &rhs;
            }
        }

        impl AddAssign<$rhs> for $lhs {
            #[inline]
            fn add_assign(&mut self, rhs: $rhs) {
                *self = &*self + &rhs;
            }
        }

        impl<'b> SubAssign<&'b $rhs> for $lhs {
            #[inline]
            fn sub_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self - rhs;
            }
        }

        impl<'b> AddAssign<&'b $rhs> for $lhs {
            #[inline]
            fn add_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self + rhs;
            }
        }
    };
}

macro_rules! impl_binops_multiplicative {
    ($lhs:ident, $rhs:ident) => {
        impl_binops_multiplicative_mixed!($lhs, $rhs, $lhs);

        impl MulAssign<$rhs> for $lhs {
            #[inline]
            fn mul_assign(&mut self, rhs: $rhs) {
                *self = &*self * &rhs;
            }
        }

        impl<'b> MulAssign<&'b $rhs> for $lhs {
            #[inline]
            fn mul_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self * rhs;
            }
        }
    };
}

macro_rules! osswu_helper {
    ($F:ident, $u:expr, $xi:expr, $ellp_a:expr, $ellp_b:expr) => {{
        use core::ops::{AddAssign, MulAssign};
        use subtle::ConditionallySelectable;

        let usq = ($u).square();

        let (nd_common, xi_usq, xi2_u4) = {
            let mut tmp = usq;
            tmp.mul_assign($xi); // xi * u^2
            let tmp2 = tmp;
            tmp = tmp.square(); // xi^2 * u^4
            let tmp3 = tmp;
            tmp.add_assign(&tmp2); // xi^2 * u^4 + xi * u^2
            (tmp, tmp2, tmp3)
        };

        // nd_common = xi^2 * u^4 + xi * u^2

        let x0_num = {
            let mut tmp = nd_common;
            tmp.add_assign(&$F::one()); // 1 + nd_common
            tmp.mul_assign($ellp_b); // B * (1 + nd_common)
            tmp
        };

        let x0_den = {
            let mut ifnz = *($ellp_a);
            let mut ifz = ifnz;
            ifz.mul_assign($xi);
            ifnz.mul_assign(&nd_common);
            ifnz = ifnz.neg();
            $F::conditional_select(&ifnz, &ifz, nd_common.is_zero())
        };
        let x0_densq = x0_den.square(); // x0_den^2

        // compute g(X0(u))
        let gx0_den = {
            let mut tmp = x0_densq; // x0_den^2
            tmp.mul_assign(&x0_den);
            tmp // x0_den ^ 3
        };

        let gx0_num = {
            let mut tmp1 = gx0_den;
            tmp1.mul_assign($ellp_b); // B * x0_den^3
            let mut tmp2 = x0_densq; // x0_den^2
            tmp2.mul_assign(&x0_num); // x0_num * x0_den^2
            tmp2.mul_assign($ellp_a); // A * x0_num * x0_den^2
            tmp1.add_assign(&tmp2); // ^^^ + B * x0_den^3
            tmp2 = x0_num.square(); // x0_num^2
            tmp2.mul_assign(&x0_num); // x0_num^3
            tmp1.add_assign(&tmp2); // x0_num^3 + A * x0_num * x0_den^2 + B * x0_den^3
            tmp1
        };

        [usq, xi_usq, xi2_u4, x0_num, x0_den, gx0_num, gx0_den]
    }};
}

// #[cfg(test)]
// /// Check that the point (X : Y : Z)==(X/Z^2, Y/Z^3) is on E: y^2 = x^3 + ELLP_A * x + ELLP_B.
// macro_rules! check_g_prime {
//     ($x:expr, $y:expr, $z:expr, $a:expr, $b:expr) => {
//         use std::ops::{AddAssign, MulAssign};

//         let lhs = ($y).square(); // y^2
//         let rhs = {
//             // x^3 + A x z^4 + B z^6
//             let zsq = ($z).square();
//             let z4 = zsq.square();

//             let mut tmp1 = ($x).square();
//             tmp1.mul_assign($x); // x^3

//             let mut tmp2 = *($x);
//             tmp2.mul_assign(&z4);
//             tmp2.mul_assign($a);
//             tmp1.add_assign(&tmp2); // + A x z^4

//             tmp2 = z4;
//             tmp2.mul_assign(&zsq);
//             tmp2.mul_assign($b);
//             tmp1.add_assign(&tmp2); // + B z^6

//             tmp1
//         };

//         assert_eq!(lhs, rhs, "not prime");
//     };
// }
