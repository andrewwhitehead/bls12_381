use crate::fp::*;
use crate::fp2::*;
use crate::fp6::*;

use core::fmt;
use core::ops;

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "pairings")]
use rand_core::RngCore;

pub(crate) const FP12_FROBENIUS_COEFFS: [Fp2; 3] = [
    Fp2 {
        c0: Fp::from_hex_unchecked(
            "1904d3bf02bb0667\
            c231beb4202c0d1f\
            0fd603fd3cbd5f4f\
            7b2443d784bab9c4\
            f67ea53d63e7813d\
            8d0775ed92235fb8",
        ),
        c1: Fp::from_hex_unchecked(
            "00fc3e2b36c4e032\
             88e9e902231f9fb8\
             54a14787b6c7b36f\
             ec0c8ec971f63c5f\
             282d5ac14d6c7ec2\
             2cf78a126ddc4af3",
        ),
    },
    Fp2 {
        c0: Fp::from_hex_unchecked(
            "135203e60180a68e\
            e2e9c448d77a2cd9\
            1c3dedd930b1cf60\
            ef396489f61eb45e\
            304466cf3e67fa0a\
            f1ee7b04121bdea2",
        ),
        c1: Fp::from_hex_unchecked(
            "06af0e0437ff400b\
             6831e36d6bd17ffe\
             48395dabc2d3435e\
             77f76e17009241c5\
             ee67992f72ec05f4\
             c81084fbede3cc09",
        ),
    },
    Fp2 {
        c0: Fp::from_hex_unchecked(
            "144e4211384586c1\
            6bd3ad4afa99cc91\
            70df3560e77982d0\
            db45f3536814f0bd\
            5871c1908bd478cd\
            1ee605167ff82995",
        ),
        c1: Fp::from_hex_unchecked(
            "05b2cfd9013a5fd8\
             df47fa6b48b1e045\
             f39816240c0b8fee\
             8beadf4d8e9c0566\
             c63a3e6e257f8732\
             9b18fae980078116",
        ),
    },
];

/// This represents an element $c_0 + c_1 w$ of $\mathbb{F}_{p^12} = \mathbb{F}_{p^6} / w^2 - v$.
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl From<Fp> for Fp12 {
    fn from(f: Fp) -> Fp12 {
        Fp12 {
            c0: Fp6::from(f),
            c1: Fp6::zero(),
        }
    }
}

impl From<Fp2> for Fp12 {
    fn from(f: Fp2) -> Fp12 {
        Fp12 {
            c0: Fp6::from(f),
            c1: Fp6::zero(),
        }
    }
}

impl From<Fp6> for Fp12 {
    fn from(f: Fp6) -> Fp12 {
        Fp12 {
            c0: f,
            c1: Fp6::zero(),
        }
    }
}

impl PartialEq for Fp12 {
    fn eq(&self, other: &Fp12) -> bool {
        self.ct_eq(other).into()
    }
}

impl Copy for Fp12 {}
impl Clone for Fp12 {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl Default for Fp12 {
    fn default() -> Self {
        Fp12::zero()
    }
}

#[cfg(feature = "zeroize")]
impl zeroize::DefaultIsZeroes for Fp12 {}

impl fmt::Debug for Fp12 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*w", self.c0, self.c1)
    }
}

impl ConditionallySelectable for Fp12 {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fp12 {
            c0: Fp6::conditional_select(&a.c0, &b.c0, choice),
            c1: Fp6::conditional_select(&a.c1, &b.c1, choice),
        }
    }
}

impl ConstantTimeEq for Fp12 {
    #[inline(always)]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1)
    }
}

impl Fp12 {
    #[inline]
    pub const fn zero() -> Self {
        Fp12 {
            c0: Fp6::zero(),
            c1: Fp6::zero(),
        }
    }

    #[inline]
    pub const fn one() -> Self {
        Fp12 {
            c0: Fp6::one(),
            c1: Fp6::zero(),
        }
    }

    #[cfg(feature = "pairings")]
    pub(crate) fn random(mut rng: impl RngCore) -> Self {
        Fp12 {
            c0: Fp6::random(&mut rng),
            c1: Fp6::random(&mut rng),
        }
    }

    #[inline(always)]
    pub fn mul_by_014(&self, c0: &Fp2, c1: &Fp2, c4: &Fp2) -> Fp12 {
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let c1 = (self.c1 + self.c0).mul_by_01(c0, &(c1 + c4)) - (aa + bb);
        let c0 = bb.mul_by_nonresidue() + aa;

        Fp12 { c0, c1 }
    }

    #[inline(always)]
    pub fn is_zero(&self) -> Choice {
        self.c0.is_zero() & self.c1.is_zero()
    }

    #[inline(always)]
    pub const fn conjugate(&self) -> Self {
        Fp12 {
            c0: self.c0,
            c1: self.c1.neg(),
        }
    }

    /// Raises this element to p^n.
    #[inline]
    pub fn frobenius_map(&self, n: usize) -> Self {
        let c0 = self.c0.frobenius_map(n);
        let c1 = self.c1.frobenius_map(n);

        // gamma = (u + 1)^((p^n - 1) / 6)
        // this is simply a roundabout way of reducing the number
        // of duplicate precomputed values.
        let gamma = if n % 2 == 0 {
            Fp2 {
                c0: FP6_FROBENIUS_COEFFS_2[5 - ((n / 2 + 5) % 6)],
                c1: Fp::zero(),
            }
        } else {
            let nb = ((n - 1) / 2) % 12;
            let coeff = FP12_FROBENIUS_COEFFS[nb % 3];
            if nb >= 3 {
                Fp2 {
                    c0: coeff.c1,
                    c1: coeff.c0,
                }
            } else {
                coeff
            }
        };

        let c1 = Fp6 {
            c0: c1.c0.mul(&gamma),
            c1: c1.c1.mul(&gamma),
            c2: c1.c2.mul(&gamma),
        };

        Fp12 { c0, c1 }
    }

    #[inline]
    pub fn square(&self) -> Self {
        let ab = self.c0 * self.c1;
        let c0 = self.c0 + self.c1.mul_by_nonresidue();
        let c0 = c0 * (self.c0 + self.c1);
        let c0 = c0 - (ab + ab.mul_by_nonresidue());
        let c1 = ab.double();

        Fp12 { c0, c1 }
    }

    #[inline]
    pub fn invert(&self) -> CtOption<Self> {
        (self.c0.square() - self.c1.square().mul_by_nonresidue())
            .invert()
            .map(|t| Fp12 {
                c0: self.c0 * t,
                c1: self.c1 * -t,
            })
    }
}

impl<'a, 'b> ops::Mul<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn mul(self, other: &'b Fp12) -> Self::Output {
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;
        let c0 = bb.mul_by_nonresidue();
        let c0 = c0 + aa;
        let o = other.c0 + other.c1;
        let c1 = self.c1 + self.c0;
        let c1 = c1 * o;
        let c1 = c1 - (aa + bb);

        Fp12 { c0, c1 }
    }
}

impl<'a, 'b> ops::Add<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn add(self, rhs: &'b Fp12) -> Self::Output {
        Fp12 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl<'a> ops::Neg for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Self::Output {
        Fp12 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

impl ops::Neg for Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<'a, 'b> ops::Sub<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn sub(self, rhs: &'b Fp12) -> Self::Output {
        Fp12 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl_binops_additive!(Fp12, Fp12);
impl_binops_multiplicative!(Fp12, Fp12);

#[test]
fn test_arithmetic() {
    use crate::fp::*;
    use crate::fp2::*;

    let a = Fp12 {
        c0: Fp6 {
            c0: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x47f9_cb98_b1b8_2d58,
                    0x5fe9_11eb_a3aa_1d9d,
                    0x96bf_1b5f_4dd8_1db3,
                    0x8100_d27c_c925_9f5b,
                    0xafa2_0b96_7464_0eab,
                    0x09bb_cea7_d8d9_497d,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x0303_cb98_b166_2daa,
                    0xd931_10aa_0a62_1d5a,
                    0xbfa9_820c_5be4_a468,
                    0x0ba3_643e_cb05_a348,
                    0xdc35_34bb_1f1c_25a6,
                    0x06c3_05bb_19c0_e1c1,
                ]),
            },
            c1: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x46f9_cb98_b162_d858,
                    0x0be9_109c_f7aa_1d57,
                    0xc791_bc55_fece_41d2,
                    0xf84c_5770_4e38_5ec2,
                    0xcb49_c1d9_c010_e60f,
                    0x0acd_b8e1_58bf_e3c8,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x8aef_cb98_b15f_8306,
                    0x3ea1_108f_e4f2_1d54,
                    0xcf79_f69f_a1b7_df3b,
                    0xe4f5_4aa1_d16b_1a3c,
                    0xba5e_4ef8_6105_a679,
                    0x0ed8_6c07_97be_e5cf,
                ]),
            },
            c2: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xcee5_cb98_b15c_2db4,
                    0x7159_1082_d23a_1d51,
                    0xd762_30e9_44a1_7ca4,
                    0xd19e_3dd3_549d_d5b6,
                    0xa972_dc17_01fa_66e3,
                    0x12e3_1f2d_d6bd_e7d6,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xad2a_cb98_b173_2d9d,
                    0x2cfd_10dd_0696_1d64,
                    0x0739_6b86_c6ef_24e8,
                    0xbd76_e2fd_b1bf_c820,
                    0x6afe_a7f6_de94_d0d5,
                    0x1099_4b0c_5744_c040,
                ]),
            },
        },
        c1: Fp6 {
            c0: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x47f9_cb98_b1b8_2d58,
                    0x5fe9_11eb_a3aa_1d9d,
                    0x96bf_1b5f_4dd8_1db3,
                    0x8100_d27c_c925_9f5b,
                    0xafa2_0b96_7464_0eab,
                    0x09bb_cea7_d8d9_497d,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x0303_cb98_b166_2daa,
                    0xd931_10aa_0a62_1d5a,
                    0xbfa9_820c_5be4_a468,
                    0x0ba3_643e_cb05_a348,
                    0xdc35_34bb_1f1c_25a6,
                    0x06c3_05bb_19c0_e1c1,
                ]),
            },
            c1: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x46f9_cb98_b162_d858,
                    0x0be9_109c_f7aa_1d57,
                    0xc791_bc55_fece_41d2,
                    0xf84c_5770_4e38_5ec2,
                    0xcb49_c1d9_c010_e60f,
                    0x0acd_b8e1_58bf_e3c8,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x8aef_cb98_b15f_8306,
                    0x3ea1_108f_e4f2_1d54,
                    0xcf79_f69f_a1b7_df3b,
                    0xe4f5_4aa1_d16b_1a3c,
                    0xba5e_4ef8_6105_a679,
                    0x0ed8_6c07_97be_e5cf,
                ]),
            },
            c2: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xcee5_cb98_b15c_2db4,
                    0x7159_1082_d23a_1d51,
                    0xd762_30e9_44a1_7ca4,
                    0xd19e_3dd3_549d_d5b6,
                    0xa972_dc17_01fa_66e3,
                    0x12e3_1f2d_d6bd_e7d6,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xad2a_cb98_b173_2d9d,
                    0x2cfd_10dd_0696_1d64,
                    0x0739_6b86_c6ef_24e8,
                    0xbd76_e2fd_b1bf_c820,
                    0x6afe_a7f6_de94_d0d5,
                    0x1099_4b0c_5744_c040,
                ]),
            },
        },
    };

    let b = Fp12 {
        c0: Fp6 {
            c0: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x47f9_cb98_b1b8_2d58,
                    0x5fe9_11eb_a3aa_1d9d,
                    0x96bf_1b5f_4dd8_1db3,
                    0x8100_d272_c925_9f5b,
                    0xafa2_0b96_7464_0eab,
                    0x09bb_cea7_d8d9_497d,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x0303_cb98_b166_2daa,
                    0xd931_10aa_0a62_1d5a,
                    0xbfa9_820c_5be4_a468,
                    0x0ba3_643e_cb05_a348,
                    0xdc35_34bb_1f1c_25a6,
                    0x06c3_05bb_19c0_e1c1,
                ]),
            },
            c1: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x46f9_cb98_b162_d858,
                    0x0be9_109c_f7aa_1d57,
                    0xc791_bc55_fece_41d2,
                    0xf84c_5770_4e38_5ec2,
                    0xcb49_c1d9_c010_e60f,
                    0x0acd_b8e1_58bf_e348,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x8aef_cb98_b15f_8306,
                    0x3ea1_108f_e4f2_1d54,
                    0xcf79_f69f_a1b7_df3b,
                    0xe4f5_4aa1_d16b_1a3c,
                    0xba5e_4ef8_6105_a679,
                    0x0ed8_6c07_97be_e5cf,
                ]),
            },
            c2: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xcee5_cb98_b15c_2db4,
                    0x7159_1082_d23a_1d51,
                    0xd762_30e9_44a1_7ca4,
                    0xd19e_3dd3_549d_d5b6,
                    0xa972_dc17_01fa_66e3,
                    0x12e3_1f2d_d6bd_e7d6,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xad2a_cb98_b173_2d9d,
                    0x2cfd_10dd_0696_1d64,
                    0x0739_6b86_c6ef_24e8,
                    0xbd76_e2fd_b1bf_c820,
                    0x6afe_a7f6_de94_d0d5,
                    0x1099_4b0c_5744_c040,
                ]),
            },
        },
        c1: Fp6 {
            c0: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x47f9_cb98_b1b8_2d58,
                    0x5fe9_11eb_a3aa_1d9d,
                    0x96bf_1b5f_4dd2_1db3,
                    0x8100_d27c_c925_9f5b,
                    0xafa2_0b96_7464_0eab,
                    0x09bb_cea7_d8d9_497d,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x0303_cb98_b166_2daa,
                    0xd931_10aa_0a62_1d5a,
                    0xbfa9_820c_5be4_a468,
                    0x0ba3_643e_cb05_a348,
                    0xdc35_34bb_1f1c_25a6,
                    0x06c3_05bb_19c0_e1c1,
                ]),
            },
            c1: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x46f9_cb98_b162_d858,
                    0x0be9_109c_f7aa_1d57,
                    0xc791_bc55_fece_41d2,
                    0xf84c_5770_4e38_5ec2,
                    0xcb49_c1d9_c010_e60f,
                    0x0acd_b8e1_58bf_e3c8,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x8aef_cb98_b15f_8306,
                    0x3ea1_108f_e4f2_1d54,
                    0xcf79_f69f_a117_df3b,
                    0xe4f5_4aa1_d16b_1a3c,
                    0xba5e_4ef8_6105_a679,
                    0x0ed8_6c07_97be_e5cf,
                ]),
            },
            c2: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xcee5_cb98_b15c_2db4,
                    0x7159_1082_d23a_1d51,
                    0xd762_30e9_44a1_7ca4,
                    0xd19e_3dd3_549d_d5b6,
                    0xa972_dc17_01fa_66e3,
                    0x12e3_1f2d_d6bd_e7d6,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xad2a_cb98_b173_2d9d,
                    0x2cfd_10dd_0696_1d64,
                    0x0739_6b86_c6ef_24e8,
                    0xbd76_e2fd_b1bf_c820,
                    0x6afe_a7f6_de94_d0d5,
                    0x1099_4b0c_5744_c040,
                ]),
            },
        },
    };

    let c = Fp12 {
        c0: Fp6 {
            c0: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x47f9_cb98_71b8_2d58,
                    0x5fe9_11eb_a3aa_1d9d,
                    0x96bf_1b5f_4dd8_1db3,
                    0x8100_d27c_c925_9f5b,
                    0xafa2_0b96_7464_0eab,
                    0x09bb_cea7_d8d9_497d,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x0303_cb98_b166_2daa,
                    0xd931_10aa_0a62_1d5a,
                    0xbfa9_820c_5be4_a468,
                    0x0ba3_643e_cb05_a348,
                    0xdc35_34bb_1f1c_25a6,
                    0x06c3_05bb_19c0_e1c1,
                ]),
            },
            c1: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x46f9_cb98_b162_d858,
                    0x0be9_109c_f7aa_1d57,
                    0x7791_bc55_fece_41d2,
                    0xf84c_5770_4e38_5ec2,
                    0xcb49_c1d9_c010_e60f,
                    0x0acd_b8e1_58bf_e3c8,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x8aef_cb98_b15f_8306,
                    0x3ea1_108f_e4f2_1d54,
                    0xcf79_f69f_a1b7_df3b,
                    0xe4f5_4aa1_d16b_133c,
                    0xba5e_4ef8_6105_a679,
                    0x0ed8_6c07_97be_e5cf,
                ]),
            },
            c2: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xcee5_cb98_b15c_2db4,
                    0x7159_1082_d23a_1d51,
                    0xd762_40e9_44a1_7ca4,
                    0xd19e_3dd3_549d_d5b6,
                    0xa972_dc17_01fa_66e3,
                    0x12e3_1f2d_d6bd_e7d6,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xad2a_cb98_b173_2d9d,
                    0x2cfd_10dd_0696_1d64,
                    0x0739_6b86_c6ef_24e8,
                    0xbd76_e2fd_b1bf_c820,
                    0x6afe_a7f6_de94_d0d5,
                    0x1099_4b0c_1744_c040,
                ]),
            },
        },
        c1: Fp6 {
            c0: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x47f9_cb98_b1b8_2d58,
                    0x5fe9_11eb_a3aa_1d9d,
                    0x96bf_1b5f_4dd8_1db3,
                    0x8100_d27c_c925_9f5b,
                    0xafa2_0b96_7464_0eab,
                    0x09bb_cea7_d8d9_497d,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x0303_cb98_b166_2daa,
                    0xd931_10aa_0a62_1d5a,
                    0xbfa9_820c_5be4_a468,
                    0x0ba3_643e_cb05_a348,
                    0xdc35_34bb_1f1c_25a6,
                    0x06c3_05bb_19c0_e1c1,
                ]),
            },
            c1: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0x46f9_cb98_b162_d858,
                    0x0be9_109c_f7aa_1d57,
                    0xc791_bc55_fece_41d2,
                    0xf84c_5770_4e38_5ec2,
                    0xcb49_c1d3_c010_e60f,
                    0x0acd_b8e1_58bf_e3c8,
                ]),
                c1: Fp::from_raw_unchecked([
                    0x8aef_cb98_b15f_8306,
                    0x3ea1_108f_e4f2_1d54,
                    0xcf79_f69f_a1b7_df3b,
                    0xe4f5_4aa1_d16b_1a3c,
                    0xba5e_4ef8_6105_a679,
                    0x0ed8_6c07_97be_e5cf,
                ]),
            },
            c2: Fp2 {
                c0: Fp::from_raw_unchecked([
                    0xcee5_cb98_b15c_2db4,
                    0x7159_1082_d23a_1d51,
                    0xd762_30e9_44a1_7ca4,
                    0xd19e_3dd3_549d_d5b6,
                    0xa972_dc17_01fa_66e3,
                    0x12e3_1f2d_d6bd_e7d6,
                ]),
                c1: Fp::from_raw_unchecked([
                    0xad2a_cb98_b173_2d9d,
                    0x2cfd_10dd_0696_1d64,
                    0x0739_6b86_c6ef_24e8,
                    0xbd76_e2fd_b1bf_c820,
                    0x6afe_a7f6_de94_d0d5,
                    0x1099_4b0c_5744_1040,
                ]),
            },
        },
    };

    // because a and b and c are similar to each other and
    // I was lazy, this is just some arbitrary way to make
    // them a little more different
    let a = a.square().invert().unwrap().square() + c;
    let b = b.square().invert().unwrap().square() + a;
    let c = c.square().invert().unwrap().square() + b;

    assert_eq!(a.square(), a * a);
    assert_eq!(b.square(), b * b);
    assert_eq!(c.square(), c * c);

    assert_eq!((a + b) * c.square(), (c * c * a) + (c * c * b));

    assert_eq!(
        a.invert().unwrap() * b.invert().unwrap(),
        (a * b).invert().unwrap()
    );
    assert_eq!(a.invert().unwrap() * a, Fp12::one());

    assert_ne!(a, a.frobenius_map(1));
    for i in 0..12 {
        assert_eq!(a.frobenius_map(i).frobenius_map(1), a.frobenius_map(i + 1));
    }
    assert_eq!(a, a.frobenius_map(12));

    let d = Fp12 {
        c0: Fp6 {
            c0: c.c0.c0,
            c1: c.c0.c1,
            c2: Fp2::zero(),
        },
        c1: Fp6 {
            c0: Fp2::zero(),
            c1: c.c1.c1,
            c2: Fp2::zero(),
        },
    };

    assert_eq!(a.mul_by_014(&d.c0.c0, &d.c0.c1, &d.c1.c1), a * d);
}

#[cfg(feature = "zeroize")]
#[test]
fn test_zeroize() {
    use zeroize::Zeroize;

    let mut a = Fp12::one();
    a.zeroize();
    assert!(bool::from(a.is_zero()));
}
