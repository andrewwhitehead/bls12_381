use subtle::{Choice, ConditionallyNegatable, ConditionallySelectable, ConstantTimeEq};

use super::chain::chain_p2m9div16;
use super::MessageToField;
use crate::{fp::Fp, fp2::Fp2, g2::G2Projective};
use crate::{
    generic_array::{
        typenum::{U128, U64},
        GenericArray,
    },
    G1Projective,
};

/// Coefficients of the 3-isogeny x map's numerator
const ISO3_XNUM: [Fp2; 4] = [
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x47f671c71ce05e62,
            0x06dd57071206393e,
            0x7c80cd2af3fd71a2,
            0x048103ea9e6cd062,
            0xc54516acc8d037f6,
            0x13808f550920ea41,
        ]),
        c1: Fp::from_raw_unchecked([
            0x47f671c71ce05e62,
            0x06dd57071206393e,
            0x7c80cd2af3fd71a2,
            0x048103ea9e6cd062,
            0xc54516acc8d037f6,
            0x13808f550920ea41,
        ]),
    },
    Fp2 {
        c0: Fp::zero(),
        c1: Fp::from_raw_unchecked([
            0x5fe55555554c71d0,
            0x873fffdd236aaaa3,
            0x6a6b4619b26ef918,
            0x21c2888408874945,
            0x2836cda7028cabc5,
            0x0ac73310a7fd5abd,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x0a0c5555555971c3,
            0xdb0c00101f9eaaae,
            0xb1fb2f941d797997,
            0xd3960742ef416e1c,
            0xb70040e2c20556f4,
            0x149d7861e581393b,
        ]),
        c1: Fp::from_raw_unchecked([
            0xaff2aaaaaaa638e8,
            0x439fffee91b55551,
            0xb535a30cd9377c8c,
            0x90e144420443a4a2,
            0x941b66d3814655e2,
            0x0563998853fead5e,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x40aac71c71c725ed,
            0x190955557a84e38e,
            0xd817050a8f41abc3,
            0xd86485d4c87f6fb1,
            0x696eb479f885d059,
            0x198e1a74328002d2,
        ]),
        c1: Fp::zero(),
    },
];

/// Coefficients of the 3-isogeny x map's denominator
const ISO3_XDEN: [Fp2; 3] = [
    Fp2 {
        c0: Fp::zero(),
        c1: Fp::from_raw_unchecked([
            0x1f3affffff13ab97,
            0xf25bfc611da3ff3e,
            0xca3757cb3819b208,
            0x3e6427366f8cec18,
            0x03977bc86095b089,
            0x04f69db13f39a952,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x447600000027552e,
            0xdcb8009a43480020,
            0x6f7ee9ce4a6e8b59,
            0xb10330b7c0a95bc6,
            0x6140b1fcfb1e54b7,
            0x0381be097f0bb4e1,
        ]),
        c1: Fp::from_raw_unchecked([
            0x7588ffffffd8557d,
            0x41f3ff646e0bffdf,
            0xf7b1e8d2ac426aca,
            0xb3741acd32dbb6f8,
            0xe9daf5b9482d581f,
            0x167f53e0ba7431b8,
        ]),
    },
    Fp2::one(),
];

/// Coefficients of the 3-isogeny y map's numerator
const ISO3_YNUM: [Fp2; 4] = [
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x96d8f684bdfc77be,
            0xb530e4f43b66d0e2,
            0x184a88ff379652fd,
            0x57cb23ecfae804e1,
            0x0fd2e39eada3eba9,
            0x08c8055e31c5d5c3,
        ]),
        c1: Fp::from_raw_unchecked([
            0x96d8f684bdfc77be,
            0xb530e4f43b66d0e2,
            0x184a88ff379652fd,
            0x57cb23ecfae804e1,
            0x0fd2e39eada3eba9,
            0x08c8055e31c5d5c3,
        ]),
    },
    Fp2 {
        c0: Fp::zero(),
        c1: Fp::from_raw_unchecked([
            0xbf0a71c71c91b406,
            0x4d6d55d28b7638fd,
            0x9d82f98e5f205aee,
            0xa27aa27b1d1a18d5,
            0x02c3b2b2d2938e86,
            0x0c7d13420b09807f,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0xd7f9555555531c74,
            0x21cffff748daaaa8,
            0x5a9ad1866c9bbe46,
            0x4870a2210221d251,
            0x4a0db369c0a32af1,
            0x02b1ccc429ff56af,
        ]),
        c1: Fp::from_raw_unchecked([
            0xe205aaaaaaac8e37,
            0xfcdc000768795556,
            0x0c96011a8a1537dd,
            0x1c06a963f163406e,
            0x010df44c82a881e6,
            0x174f45260f808feb,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0xa470bda12f67f35c,
            0xc0fe38e23327b425,
            0xc9d3d0f2c6f0678d,
            0x1c55c9935b5a982e,
            0x27f6c0e2f0746764,
            0x117c5e6e28aa9054,
        ]),
        c1: Fp::zero(),
    },
];

/// Coefficients of the 3-isogeny y map's denominator
const ISO3_YDEN: [Fp2; 4] = [
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x0162fffffa765adf,
            0x8f7bea480083fb75,
            0x561b3c2259e93611,
            0x11e19fc1a9c875d5,
            0xca713efc00367660,
            0x03c6a03d41da1151,
        ]),
        c1: Fp::from_raw_unchecked([
            0x0162fffffa765adf,
            0x8f7bea480083fb75,
            0x561b3c2259e93611,
            0x11e19fc1a9c875d5,
            0xca713efc00367660,
            0x03c6a03d41da1151,
        ]),
    },
    Fp2 {
        c0: Fp::zero(),
        c1: Fp::from_raw_unchecked([
            0x5db0fffffd3b02c5,
            0xd713f52358ebfdba,
            0x5ea60761a84d161a,
            0xbb2c75a34ea6c44a,
            0x0ac6735921c1119b,
            0x0ee3d913bdacfbf6,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x66b10000003affc5,
            0xcb1400e764ec0030,
            0xa73e5eb56fa5d106,
            0x8984c913a0fe09a9,
            0x11e10afb78ad7f13,
            0x05429d0e3e918f52,
        ]),
        c1: Fp::from_raw_unchecked([
            0x534dffffffc4aae6,
            0x5397ff174c67ffcf,
            0xbff273eb870b251d,
            0xdaf2827152870915,
            0x393a9cbaca9e2dc3,
            0x14be74dbfaee5748,
        ]),
    },
    Fp2::one(),
];

const SSWU_ELLP_A: Fp2 = Fp2 {
    c0: Fp::zero(),
    c1: Fp::from_raw_unchecked([
        0xe53a000003135242,
        0x01080c0fdef80285,
        0xe7889edbe340f6bd,
        0x0b51375126310601,
        0x02d6985717c744ab,
        0x1220b4e979ea5467,
    ]),
};

const SSWU_ELLP_B: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0x22ea00000cf89db2,
        0x6ec832df71380aa4,
        0x6e1b94403db5a66e,
        0x75bf3c53a79473ba,
        0x3dd3a569412c0a34,
        0x125cdb5e74dc4fd1,
    ]),
    c1: Fp::from_raw_unchecked([
        0x22ea00000cf89db2,
        0x6ec832df71380aa4,
        0x6e1b94403db5a66e,
        0x75bf3c53a79473ba,
        0x3dd3a569412c0a34,
        0x125cdb5e74dc4fd1,
    ]),
};

const SSWU_XI: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0x87ebfffffff9555c,
        0x656fffe5da8ffffa,
        0x0fd0749345d33ad2,
        0xd951e663066576f4,
        0xde291a3d41e980d3,
        0x0815664c7dfe040d,
    ]),
    c1: Fp::from_raw_unchecked([
        0x43f5fffffffcaaae,
        0x32b7fff2ed47fffd,
        0x07e83a49a2e99d69,
        0xeca8f3318332bb7a,
        0xef148d1ea0f4c069,
        0x040ab3263eff0206,
    ]),
};

const SSWU_ETAS: [Fp2; 4] = [
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x05e514668ac736d2,
            0x9089b4d6b84f3ea5,
            0x603c384c224a8b32,
            0xf3257909536afea6,
            0x5c5cdbabae656d81,
            0x075bfa0863c987e9,
        ]),
        c1: Fp::from_raw_unchecked([
            0x338d9bfe08087330,
            0x7b8e48b2bd83cefe,
            0x530dad5d306b5be7,
            0x5a4d7e8e6c408b6d,
            0x6258f7a6232cab9b,
            0x0b985811cce14db5,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x86716401f7f7377b,
            0xa31db74bf3d03101,
            0x14232543c6459a3c,
            0x0a29ccf687448752,
            0xe8c2b010201f013c,
            0x0e68b9d86c9e98e4,
        ]),
        c1: Fp::from_raw_unchecked([
            0x05e514668ac736d2,
            0x9089b4d6b84f3ea5,
            0x603c384c224a8b32,
            0xf3257909536afea6,
            0x5c5cdbabae656d81,
            0x075bfa0863c987e9,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x718fdad24ee1d90f,
            0xa58c025bed8276af,
            0x0c3a10230ab7976f,
            0xf0c54df5c8f275e1,
            0x4ec2478c28baf465,
            0x1129373a90c508e6,
        ]),
        c1: Fp::from_raw_unchecked([
            0x019af5f980a3680c,
            0x4ed7da0e66063afa,
            0x600354723b5d9972,
            0x8b2f958b20d09d72,
            0x0474938f02d461db,
            0x0dcf8b9e0684ab1c,
        ]),
    },
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0xb8640a067f5c429f,
            0xcfd425f04b4dc505,
            0x072d7e2ebb535cb1,
            0xd947b5f9d2b4754d,
            0x46a7142740774afb,
            0x0c31864c32fb3b7e,
        ]),
        c1: Fp::from_raw_unchecked([
            0x718fdad24ee1d90f,
            0xa58c025bed8276af,
            0x0c3a10230ab7976f,
            0xf0c54df5c8f275e1,
            0x4ec2478c28baf465,
            0x1129373a90c508e6,
        ]),
    },
];

const SSWU_RV1: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0x7bcfa7a25aa30fda,
        0xdc17dec12a927e7c,
        0x2f088dd86b4ebef1,
        0xd1ca2087da74d4a7,
        0x2da2596696cebc1d,
        0x0e2b7eedbbfd87d2,
    ]),
    c1: Fp::from_raw_unchecked([
        0x7bcfa7a25aa30fda,
        0xdc17dec12a927e7c,
        0x2f088dd86b4ebef1,
        0xd1ca2087da74d4a7,
        0x2da2596696cebc1d,
        0x0e2b7eedbbfd87d2,
    ]),
};

impl MessageToField for G2Projective {
    type InputLength = U128;
    type Pt = Fp2;

    fn input_okm(okm: &GenericArray<u8, U128>) -> Fp2 {
        let c0 = <G1Projective as MessageToField>::input_okm(GenericArray::<u8, U64>::from_slice(
            &okm[..64],
        ));
        let c1 = <G1Projective as MessageToField>::input_okm(GenericArray::<u8, U64>::from_slice(
            &okm[64..],
        ));
        Fp2 { c0, c1 }
    }

    fn isogeny_map(&mut self) {
        const COEFFS: [&[Fp2]; 4] = [&ISO3_XNUM, &ISO3_XDEN, &ISO3_YNUM, &ISO3_YDEN];

        // xnum, xden, ynum, yden
        let mut mapvals = [Fp2::zero(); 4];

        // unpack input point
        let G2Projective { x, y, z } = *self;

        // compute powers of z
        let zsq = z.square();
        let zpows = [z, zsq, zsq * z];

        // compute map value by Horner's rule
        for idx in 0..4 {
            let coeff = COEFFS[idx];
            let clast = coeff.len() - 1;
            mapvals[idx] = coeff[clast];
            for jdx in 0..clast {
                mapvals[idx] = mapvals[idx] * x + zpows[jdx] * coeff[clast - 1 - jdx];
            }
        }

        // x denominator is order 1 less than x numerator, so we need an extra factor of z
        mapvals[1] *= z;

        // multiply result of Y map by the y-coord, y / z
        mapvals[2] *= y;
        mapvals[3] *= z;

        // projective coordinates of resulting point
        *self = G2Projective {
            x: mapvals[0] * mapvals[3], // xnum * yden,
            y: mapvals[2] * mapvals[1], // ynum * xden,
            z: mapvals[1] * mapvals[3], // xden * yden
        }
    }

    fn osswu_map(u: &Fp2) -> Self {
        let usq = u.square();
        let xi_usq = SSWU_XI * usq;
        let xisq_u4 = xi_usq.square();
        let nd_common = xisq_u4 + xi_usq; // XI^2 * u^4 + XI * u^2
        let x_den =
            SSWU_ELLP_A * Fp2::conditional_select(&(-nd_common), &SSWU_XI, nd_common.is_zero());
        let x0_num = SSWU_ELLP_B * (Fp2::one() + nd_common); // B * (1 + (XI^2 * u^4 + XI * u^2))

        // compute g(x0(u))
        let x_densq = x_den.square();
        let gx_den = x_densq * x_den;
        // x0_num^3 + A * x0_num * x_den^2 + B * x_den^3
        let gx0_num = (x0_num.square() + SSWU_ELLP_A * x_densq) * x0_num + SSWU_ELLP_B * gx_den;

        // compute g(x0(u)) ^ ((p^2 - 9) // 16)
        let sqrt_candidate = {
            let vsq = gx_den.square(); // v^2
            let v_3 = vsq * gx_den; // v^3
            let v_4 = vsq.square(); // v^4
            let uv_7 = gx0_num * v_3 * v_4; // u v^7
            let uv_15 = uv_7 * v_4.square(); // u v^15
            uv_7 * chain_p2m9div16(&uv_15) // u v^7 (u v^15) ^ ((p^2 - 9) // 16)
        };

        // set y = sqrt_candidate * Fp2::one(), check candidate against other roots of unity
        let mut y = sqrt_candidate;
        // check Fp2(0, 1)
        let tmp = Fp2 {
            c0: -sqrt_candidate.c1,
            c1: sqrt_candidate.c0,
        };
        y.conditional_assign(&tmp, (tmp.square() * gx_den).ct_eq(&gx0_num));
        // check Fp2(RV1, RV1)
        let tmp = sqrt_candidate * SSWU_RV1;
        y.conditional_assign(&tmp, (tmp.square() * gx_den).ct_eq(&gx0_num));
        // check Fp2(RV1, -RV1)
        let tmp = Fp2 {
            c0: tmp.c1,
            c1: -tmp.c0,
        };
        y.conditional_assign(&tmp, (tmp.square() * gx_den).ct_eq(&gx0_num));

        // compute g(x1(u)) = g(x0(u)) * XI^3 * u^6
        let gx1_num = gx0_num * xi_usq * xisq_u4;
        // compute g(x1(u)) * u^3
        let sqrt_candidate = sqrt_candidate * usq * u;
        let mut eta_found = Choice::from(0u8);
        for eta in &SSWU_ETAS[..] {
            let tmp = sqrt_candidate * eta;
            let found = (tmp.square() * gx_den).ct_eq(&gx1_num);
            y.conditional_assign(&tmp, found);
            eta_found |= found;
        }

        let x_num = Fp2::conditional_select(&x0_num, &(x0_num * xi_usq), eta_found);
        // ensure sign of y and sign of u agree
        y.conditional_negate(u.sgn0() ^ y.sgn0()); // fix sign of y

        G2Projective {
            x: x_num,
            y: y * x_den,
            z: x_den,
        }
    }

    fn clear_h(&mut self) {
        *self = self.clear_cofactor();
    }
}

#[cfg(test)]
fn check_g2_prime(x: &Fp2, y: &Fp2, z: &Fp2) -> subtle::Choice {
    // (X : Y : Z)==(X/Z, Y/Z) is on E: y^2 = x^3 + ELLP_A * x + ELLP_B.
    // y^2 z = (x^3) + A (x z^2) + B z^3
    let zsq = z.square();
    (y.square() * z).ct_eq(&(x.square() * x + SSWU_ELLP_A * x * zsq + SSWU_ELLP_B * zsq * z))
}

#[test]
fn test_osswu_semirandom() {
    use group::Group;
    use rand_core::SeedableRng;
    let mut rng = rand_xorshift::XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    for _ in 0..32 {
        let input = Fp2::random(&mut rng);
        let p = G2Projective::osswu_map(&input);
        let G2Projective { x, y, z } = &p;
        check_g2_prime(x, y, z);

        // FIXME remove
        let p = G2Projective::random(&mut rng);
        assert_eq!(p + p - p, p);
        // assert_eq!(p.double(), p + p);
    }
}

#[test]
fn test_encode_to_curve_07() {
    use crate::{
        g2::G2Affine,
        hash_to_curve::{ExpandMsgXmd, HashToCurve},
    };

    struct TestCase {
        msg: &'static [u8],
        expected: [&'static str; 4],
    }
    impl TestCase {
        fn expected(&self) -> String {
            self.expected[0].to_string() + self.expected[1] + self.expected[2] + self.expected[3]
        }
    }

    const DOMAIN: &[u8] = b"BLS12381G2_XMD:SHA-256_SSWU_NU_TESTGEN";

    let cases = vec![
        TestCase {
            msg: b"",
            expected: [
        "0d4333b77becbf9f9dfa3ca928002233d1ecc854b1447e5a71f751c9042d000f42db91c1d6649a5e0ad22bd7bf7398b8",
        "027e4bfada0b47f9f07e04aec463c7371e68f2fd0c738cd517932ea3801a35acf09db018deda57387b0f270f7a219e4d",
        "0cc76dc777ea0d447e02a41004f37a0a7b1fafb6746884e8d9fc276716ccf47e4e0899548a2ec71c2bdf1a2a50e876db",
        "053674cba9ef516ddc218fedb37324e6c47de27f88ab7ef123b006127d738293c0277187f7e2f80a299a24d84ed03da7",
            ],
        },
        TestCase {
            msg: b"abc",
            expected: [
        "18f0f87b40af67c056915dbaf48534c592524e82c1c2b50c3734d02c0172c80df780a60b5683759298a3303c5d942778",
        "09349f1cb5b2e55489dcd45a38545343451cc30a1681c57acd4fb0a6db125f8352c09f4a67eb7d1d8242cb7d3405f97b",
        "10a2ba341bc689ab947b7941ce6ef39be17acaab067bd32bd652b471ab0792c53a2bd03bdac47f96aaafe96e441f63c0",
        "02f2d9deb2c7742512f5b8230bf0fd83ea42279d7d39779543c1a43b61c885982b611f6a7a24b514995e8a098496b811",
            ],
        },
        TestCase {
            msg: b"abcdef0123456789",
            expected: [
        "19808ec5930a53c7cf5912ccce1cc33f1b3dcff24a53ce1cc4cba41fd6996dbed4843ccdd2eaf6a0cd801e562718d163",
        "149fe43777d34f0d25430dea463889bd9393bdfb4932946db23671727081c629ebb98a89604f3433fba1c67d356a4af7",
        "04783e391c30c83f805ca271e353582fdf19d159f6a4c39b73acbb637a9b8ac820cfbe2738d683368a7c07ad020e3e33",
        "04c0d6793a766233b2982087b5f4a254f261003ccb3262ea7c50903eecef3e871d1502c293f9e063d7d293f6384f4551",
            ]
        },
        TestCase {
            msg: b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            expected: [
        "0b8e0094c886487870372eb6264613a6a087c7eb9804fab789be4e47a57b29eb19b1983a51165a1b5eb025865e9fc63a",
        "0804152cbf8474669ad7d1796ab92d7ca21f32d8bed70898a748ed4e4e0ec557069003732fc86866d938538a2ae95552",
        "14c80f068ece15a3936bb00c3c883966f75b4e8d9ddde809c11f781ab92d23a2d1d103ad48f6f3bb158bf3e3a4063449",
        "09e5c8242dd7281ad32c03fe4af3f19167770016255fb25ad9b67ec51d62fade31a1af101e8f6172ec2ee8857662be3a",
            ]
        }
    ];

    for case in cases {
        let g = <G2Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::encode_to_curve(
            &case.msg, DOMAIN,
        );
        let g_uncompressed = G2Affine::from(g).to_uncompressed();

        assert_eq!(case.expected(), hex::encode(&g_uncompressed[..]));
    }
}

#[test]
fn test_hash_to_curve_07() {
    use crate::{
        g2::G2Affine,
        hash_to_curve::{ExpandMsgXmd, HashToCurve},
    };

    struct TestCase {
        msg: &'static [u8],
        expected: [&'static str; 4],
    }
    impl TestCase {
        fn expected(&self) -> String {
            self.expected[0].to_string() + self.expected[1] + self.expected[2] + self.expected[3]
        }
    }

    const DOMAIN: &[u8] = b"BLS12381G2_XMD:SHA-256_SSWU_RO_TESTGEN";

    let cases = vec![
        TestCase {
            msg: b"",
            expected: [
        "0fbdae26f9f9586a46d4b0b70390d09064ef2afe5c99348438a3c7d9756471e015cb534204c1b6824617a85024c772dc",
        "0a650bd36ae7455cb3fe5d8bb1310594551456f5c6593aec9ee0c03d2f6cb693bd2c5e99d4e23cbaec767609314f51d3",
        "02e5cf8f9b7348428cc9e66b9a9b36fe45ba0b0a146290c3a68d92895b1af0e1f2d9f889fb412670ae8478d8abd4c5aa",
        "0d8d49e7737d8f9fc5cef7c4b8817633103faf2613016cb86a1f3fc29968fe2413e232d9208d2d74a89bf7a48ac36f83",
            ],
        },
        TestCase {
            msg: b"abc",
            expected: [
        "03578447618463deb106b60e609c6f7cc446dc6035f84a72801ba17c94cd800583b493b948eff0033f09086fdd7f6175",
        "1953ce6d4267939c7360756d9cca8eb34aac4633ef35369a7dc249445069888e7d1b3f9d2e75fbd468fbcbba7110ea02",
        "0184d26779ae9d4670aca9b267dbd4d3b30443ad05b8546d36a195686e1ccc3a59194aea05ed5bce7c3144a29ec047c4",
        "0882ab045b8fe4d7d557ebb59a63a35ac9f3d312581b509af0f8eaa2960cbc5e1e36bb969b6e22980b5cbdd0787fcf4e",
            ],
        },
        TestCase {
            msg: b"abcdef0123456789",
            expected: [
        "195fad48982e186ce3c5c82133aefc9b26d55979b6f530992a8849d4263ec5d57f7a181553c8799bcc83da44847bdc8d",
        "17b461fc3b96a30c2408958cbfa5f5927b6063a8ad199d5ebf2d7cdeffa9c20c85487204804fab53f950b2f87db365aa",
        "005cdf3d984e3391e7e969276fb4bc02323c5924a4449af167030d855acc2600cf3d4fab025432c6d868c79571a95bef",
        "174a3473a3af2d0302b9065e895ca4adba4ece6ce0b41148ba597001abb152f852dd9a96fb45c9de0a43d944746f833e",
            ]
        },
        TestCase {
            msg: b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            expected: [
        "123b6bd9feeba26dd4ad00f8bfda2718c9700dc093ea5287d7711844644eb981848316d3f3f57d5d3a652c6cdc816aca",
        "0a162306f3b0f2bb326f0c4fb0e1fea020019c3af796dcd1d7264f50ddae94cacf3cade74603834d44b9ab3d5d0a6c98",
        "05483f3b96d9252dd4fc0868344dfaf3c9d145e3387db23fa8e449304fab6a7b6ec9c15f05c0a1ea66ff0efcc03e001a",
        "15c1d4f1a685bb63ee67ca1fd96155e3d091e852a684b78d085fd34f6091e5249ddddbdcf2e7ec82ce6c04c63647eeb7",
            ]
        }
    ];

    for case in cases {
        let g = <G2Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(
            &case.msg, DOMAIN,
        );
        let g_uncompressed = G2Affine::from(g).to_uncompressed();

        assert_eq!(case.expected(), hex::encode(&g_uncompressed[..]));
    }
}
