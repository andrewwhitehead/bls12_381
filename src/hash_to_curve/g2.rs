use subtle::{ConditionallyNegatable, ConditionallySelectable, ConstantTimeEq};

use super::chain::{chain_h2_eff, chain_p2m9div16};
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
const XNUM: [Fp2; 4] = [
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
        c0: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
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
        c1: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
    },
];

/// Coefficients of the 3-isogeny x map's denominator
const XDEN: [Fp2; 3] = [
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
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
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x760900000002fffd,
            0xebf4000bc40c0002,
            0x5f48985753c758ba,
            0x77ce585370525745,
            0x5c071a97a256ec6d,
            0x15f65ec3fa80e493,
        ]),
        c1: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
    },
];

/// Coefficients of the 3-isogeny y map's numerator
const YNUM: [Fp2; 4] = [
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
        c0: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
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
        c1: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
    },
];

/// Coefficients of the 3-isogeny y map's denominator
const YDEN: [Fp2; 4] = [
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
        c0: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
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
    Fp2 {
        c0: Fp::from_raw_unchecked([
            0x760900000002fffd,
            0xebf4000bc40c0002,
            0x5f48985753c758ba,
            0x77ce585370525745,
            0x5c071a97a256ec6d,
            0x15f65ec3fa80e493,
        ]),
        c1: Fp::from_raw_unchecked([
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        ]),
    },
];

const ELLP_A: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]),
    c1: Fp::from_raw_unchecked([
        0xe53a000003135242u64,
        0x01080c0fdef80285u64,
        0xe7889edbe340f6bdu64,
        0x0b51375126310601u64,
        0x02d6985717c744abu64,
        0x1220b4e979ea5467u64,
    ]),
};

const ELLP_B: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0x22ea00000cf89db2u64,
        0x6ec832df71380aa4u64,
        0x6e1b94403db5a66eu64,
        0x75bf3c53a79473bau64,
        0x3dd3a569412c0a34u64,
        0x125cdb5e74dc4fd1u64,
    ]),
    c1: Fp::from_raw_unchecked([
        0x22ea00000cf89db2u64,
        0x6ec832df71380aa4u64,
        0x6e1b94403db5a66eu64,
        0x75bf3c53a79473bau64,
        0x3dd3a569412c0a34u64,
        0x125cdb5e74dc4fd1u64,
    ]),
};

const XI: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0x87ebfffffff9555cu64,
        0x656fffe5da8ffffau64,
        0xfd0749345d33ad2u64,
        0xd951e663066576f4u64,
        0xde291a3d41e980d3u64,
        0x815664c7dfe040du64,
    ]),
    c1: Fp::from_raw_unchecked([
        0x43f5fffffffcaaaeu64,
        0x32b7fff2ed47fffdu64,
        0x7e83a49a2e99d69u64,
        0xeca8f3318332bb7au64,
        0xef148d1ea0f4c069u64,
        0x40ab3263eff0206u64,
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

    /// Default isogeny evaluation function.
    fn isogeny_map(&mut self) {
        // self.eval_iso([&XNUM[..], &XDEN[..], &YNUM[..], &YDEN[..]]);
        *self = self.map_isogeny([&XNUM[..], &XDEN[..], &YNUM[..], &YDEN[..]]);
    }

    fn osswu_map(u: &Fp2) -> Self {
        // sqrt(-1)
        const C2: Fp2 = Fp2 {
            c0: Fp::from_raw_unchecked([0x0, 0x0, 0x0, 0x0, 0x0, 0x0]),
            c1: Fp::from_raw_unchecked([
                0x43f5fffffffcaaae,
                0x32b7fff2ed47fffd,
                0x7e83a49a2e99d69,
                0xeca8f3318332bb7a,
                0xef148d1ea0f4c069,
                0x40ab3263eff0206,
            ]),
        };
        // sqrt(c2)
        const C3: Fp2 = Fp2 {
            c0: Fp::from_raw_unchecked([
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
                0x2da2596696cebc1d,
                0xe2b7eedbbfd87d2,
            ]),
            c1: Fp::from_raw_unchecked([
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
                0x2da2596696cebc1d,
                0xe2b7eedbbfd87d2,
            ]),
        };
        // sqrt(xi^3 / c3)
        const C4: Fp2 = Fp2 {
            c0: Fp::from_raw_unchecked([
                0x486f252db11dd19c,
                0x791ffda2c3d18950,
                0x5af6c27debf95eb4,
                0x73b1fd8f2a929cde,
                0xfc59602a1a90b871,
                0x8d7daafa8baddb3,
            ]),
            c1: Fp::from_raw_unchecked([
                0xb8640a067f5c429f,
                0xcfd425f04b4dc505,
                0x72d7e2ebb535cb1,
                0xd947b5f9d2b4754d,
                0x46a7142740774afb,
                0xc31864c32fb3b7e,
            ]),
        };
        // sqrt(xi^3 / (c2 * c3))
        const C5: Fp2 = Fp2 {
            c0: Fp::from_raw_unchecked([
                0x338d9bfe08087330,
                0x7b8e48b2bd83cefe,
                0x530dad5d306b5be7,
                0x5a4d7e8e6c408b6d,
                0x6258f7a6232cab9b,
                0xb985811cce14db5,
            ]),
            c1: Fp::from_raw_unchecked([
                0xb419eb99753873d9,
                0x8e224b27f904c15a,
                0x6f49a54d4666af1,
                0x7151d27ba01a1419,
                0xeebecc0a94e63f55,
                0x12a517e1d5b65eb0,
            ]),
        };

        let mut tmp1 = u.square();
        let tmp3 = XI * tmp1;
        let mut tmp5 = tmp3.square();
        let mut xd = tmp5 + tmp3;
        let x1n = (xd + Fp2::one()) * ELLP_B;
        xd *= ELLP_A;
        xd = xd.neg();
        xd.conditional_assign(&(XI * ELLP_A), xd.is_zero()); // If xd == 0, set xd = Z * A
        let mut tmp2 = xd.square();
        let gxd = tmp2 * xd; // xd^3
        tmp2 *= ELLP_A;
        let gx1 = (x1n.square() + tmp2) * x1n + ELLP_B * gxd; // x1n^3 + A * x1n * xd^2 + B * xd^3
        let mut tmp4 = gxd.square();
        tmp2 = tmp4 * gxd; // gxd^3
        tmp4 = tmp4.square(); // gxd^4
        tmp2 = tmp2 * tmp4; // gxd^7
        tmp2 = tmp2 * gx1; // gx1 * gxd^7
        tmp4 = tmp2 * tmp4.square(); // gx1 * gxd^15
        let mut y = chain_p2m9div16(&tmp4);
        y *= tmp2; // This is almost sqrt(gx1)

        // Check the four possible square roots
        tmp4 = y * C2;
        tmp2 = tmp4.square() * gxd;
        y.conditional_assign(&tmp4, tmp2.ct_eq(&gx1));
        tmp4 = y * C3;
        tmp2 = tmp4.square() * gxd;
        y.conditional_assign(&tmp4, tmp2.ct_eq(&gx1));
        tmp4 = tmp4 * C2;
        tmp2 = tmp4.square() * gxd;
        // if x1 is square, this is its sqrt
        y.conditional_assign(&tmp4, tmp2.ct_eq(&gx1));

        let gx2 = gx1 * tmp5 * tmp3; // gx1 * XI^3 * u^6
        tmp5 = y * tmp1 * u; // this is almost sqrt(gx2)

        // Check the four possible square roots
        tmp1 = tmp5 * C4;
        tmp4 = tmp1 * C2;
        tmp2 = tmp4.square() * gxd;
        tmp1.conditional_assign(&tmp4, tmp2.ct_eq(&gx2));
        tmp4 = tmp5 * C5;
        tmp2 = tmp4.square() * gxd;
        tmp1.conditional_assign(&tmp4, tmp2.ct_eq(&gx2));
        tmp4 = tmp4 * C4;
        tmp2 = tmp4.square() * gxd;
        // if x2 is square, this is its sqrt
        tmp1.conditional_assign(&tmp4, tmp2.ct_eq(&gx2));

        tmp2 = y.square() * gxd;
        let e8 = tmp2.ct_eq(&gx1);
        y.conditional_assign(&tmp1, !e8); // choose correct y coordinate
        tmp2 = tmp3 * x1n; // x2n = x2n / xd = Z * u^2 * x1n / xd
        let xn = Fp2::conditional_select(&tmp2, &x1n, e8); // choose correct x-coordinate
        y.conditional_negate(u.sgn0() ^ y.sgn0()); // fix sign of y

        G2Projective {
            x: xn,
            y: y * xd,
            z: xd,
        }

        // let p = ((XI * XI * XI) * (c2 * c3).invert().unwrap()).sqrt();
        // println!("c5: {:#x?}", p);
        // assert!(false);

        // compute x0 and g(x0)
        // let [usq, xi_usq, xi2_u4, x0_num, x0_den, gx0_num, gx0_den] =
        //     osswu_helper!(Fp2, u, &XI, &ELLP_A, &ELLP_B);

        // // compute g(x0(u)) ^ ((p - 9) // 16)
        // let sqrt_candidate = {
        //     let mut tmp1 = gx0_den.square(); // v^2
        //     let mut tmp2 = tmp1;
        //     tmp1 = tmp1.square(); // v^4
        //     tmp2.mul_assign(&tmp1); // v^6
        //     tmp2.mul_assign(&gx0_den); // v^7
        //     tmp2.mul_assign(&gx0_num); // u v^7
        //     tmp1 = tmp1.square(); // v^8
        //     tmp1.mul_assign(&tmp2); // u v^15
        //     let tmp3 = tmp1;
        //     tmp1 = chain_p2m9div16(&tmp3); // (u v^15) ^ ((p - 9) // 16)
        //     tmp1.mul_assign(&tmp2); // u v^7 (u v^15) ^ ((p - 9) // 16)
        //     tmp1
        // };

        // let mut result = G2Projective::identity();
        // let mut done = Choice::from(0u8);

        // for root in &ROOTS_OF_UNITY[..] {
        //     let mut y0 = *root;
        //     y0.mul_assign(&sqrt_candidate);

        //     let mut tmp = y0.square();
        //     tmp.mul_assign(&gx0_den);

        //     let found = Choice::from((tmp == gx0_num) as u8) & !done;

        //     let sgn0_y_xor_u = y0.sgn0() ^ u.sgn0();
        //     y0.conditional_negate(sgn0_y_xor_u);
        //     // y0.mul_assign(&gx0_den); // y * x0_den^3 / x0_den^3 = y
        //     y0.mul_assign(&x0_den); // y * x0_den / x0_den = y

        //     // tmp = x0_num;
        //     // tmp.mul_assign(&x0_den); // x0_num * x0_den / x0_den^2 = x0_num / x0_den

        //     let found_e = G2Projective {
        //         x: x0_num, //tmp,
        //         y: y0,
        //         z: x0_den,
        //     };
        //     result.conditional_assign(&found_e, found);
        //     done |= found;
        // }

        // // If we've gotten here, g(X0(u)) is not square. Use X1 instead.
        // let x1_num = {
        //     let mut tmp = x0_num;
        //     tmp.mul_assign(&xi_usq);
        //     tmp
        // };
        // let gx1_num = {
        //     let mut tmp = xi2_u4;
        //     tmp.mul_assign(&xi_usq); // xi^3 u^6
        //     tmp.mul_assign(&gx0_num);
        //     tmp
        // };
        // let sqrt_candidate = {
        //     let mut tmp = sqrt_candidate;
        //     tmp.mul_assign(&usq);
        //     tmp.mul_assign(u);
        //     tmp
        // };
        // for eta in &ETAS[..] {
        //     let mut y1 = *eta;
        //     y1.mul_assign(&sqrt_candidate);

        //     let mut tmp = y1.square();
        //     tmp.mul_assign(&gx0_den);

        //     let found = Choice::from((tmp == gx1_num) as u8) & !done;

        //     let sgn0_y_xor_u = y1.sgn0() ^ u.sgn0();
        //     y1.conditional_negate(sgn0_y_xor_u);
        //     // y1.mul_assign(&gx0_den); // y * x0_den^3 / x0_den^3 = y
        //     y1.mul_assign(&x0_den); // y * x0_den / x0_den = y

        //     // tmp = x1_num;
        //     // tmp.mul_assign(&x0_den); // x1_num * x0_den / x0_den^2 = x1_num / x0_den

        //     let found_e = G2Projective {
        //         x: x1_num, // tmp,
        //         y: y1,
        //         z: x0_den,
        //     };
        //     result.conditional_assign(&found_e, found);
        //     done |= found;
        // }

        // if !bool::from(done) {
        //     panic!("Failed to find square root in G2 osswu_map");
        // }
        // result
    }

    fn clear_h(&mut self) {
        //*self = self.clear_cofactor();
        *self = chain_h2_eff(&self);
    }
}

#[test]
fn test_projective_iso3_zero() {
    let zero = Fp2::zero();
    let mut pt = G2Projective {
        x: Fp2::zero(),
        y: Fp2::zero(),
        z: Fp2::zero(),
    };
    pt.isogeny_map();
    assert_eq!(pt.x, zero);
    assert_eq!(pt.y, zero);
    assert_eq!(pt.z, zero);
}

#[test]
fn test_projective_iso3_one() {
    let mut pt = G2Projective {
        x: Fp2::one(),
        y: Fp2::one(),
        z: Fp2::one(),
    };
    pt.isogeny_map();
    println!("{:#x?}", pt);
    println!("{:#x?}", crate::g2::G2Affine::from(pt));
    println!(
        "{}",
        hex::encode(&crate::g2::G2Affine::from(pt).to_compressed()[..])
    );
    let ptb = pt + pt;
    pt = pt.double_safe();
    println!(
        "{}",
        hex::encode(&crate::g2::G2Affine::from(pt).to_compressed()[..])
    );
    println!(
        "{}",
        hex::encode(&crate::g2::G2Affine::from(ptb).to_compressed()[..])
    );
    assert!(bool::from(pt.is_on_curve()));

    let x_expect = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x57c6555579807bcau64,
            0xc285c71b6d7a38e3u64,
            0xde7b4e7d31a614c6u64,
            0x31b21e4af64b0e94u64,
            0x8fc02d1bfb73bf52u64,
            0x1439b899baf1b35bu64,
        ]),
        c1: Fp::from_raw_unchecked([
            0xf58daaab358a307bu64,
            0x665f8e3829a071c6u64,
            0x55c5ca596c9b3369u64,
            0xfeecf110f9110a6au64,
            0xd464b281b39bd1ccu64,
            0x0e725f493c63801cu64,
        ]),
    };
    let y_expect = Fp2 {
        c0: Fp::from_raw_unchecked([
            0xa72f3db7cb8405a4u64,
            0x221fda12b88ad097u64,
            0x71ec98c879891123u64,
            0x54f9a5b05305ae23u64,
            0xf176e62b3bde9b44u64,
            0x04d0ca6dbecbd55eu64,
        ]),
        c1: Fp::from_raw_unchecked([
            0xe1b3626ab65e39a9u64,
            0x4e79097a56dc4bd9u64,
            0xb0e977c69aa27452u64,
            0x761b0f37a1e26286u64,
            0xfbf7043de3811ad0u64,
            0x124c9ad43b6cf79bu64,
        ]),
    };
    let z_expect = Fp2 {
        c0: Fp::from_raw_unchecked([
            0xb9fefffffffebb2au64,
            0x1eabfffeb153ffffu64,
            0x6730d2a0f6b0f624u64,
            0x64774b84f38512bfu64,
            0x4b1ba7b6434bacd7u64,
            0x1a0111ea397fe69au64,
        ]),
        c1: Fp::from_raw_unchecked([
            0x00000000000065b2u64,
            0x0000000000000000u64,
            0x0000000000000000u64,
            0x0000000000000000u64,
            0x0000000000000000u64,
            0x0000000000000000u64,
        ]),
    };
    let x_expect = x_expect * z_expect;
    let y_expect = y_expect * z_expect * z_expect;
    // assert_eq!(pt.x, x_expect);
    // assert_eq!(pt.y, y_expect);
    assert_eq!(pt.z, z_expect);
}

#[test]
fn test_projective_iso3_fixed() {
    let xi = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x2b4f1b0418ec2ab9,
            0x8ccfc3bd38b8b1bd,
            0x160d21c60264b158,
            0x44d11146d827540,
            0xe9a1ff8efbfa3e55,
            0x2790421956b94db,
        ]),
        c1: Fp::from_raw_unchecked([
            0xf53130c2f21627e9,
            0x8bf41fb299f22777,
            0x22de3385337eef77,
            0xa3d650b28238d936,
            0xc26ac36ef74be788,
            0x6b0de5801bcdb55,
        ]),
    };
    let yi = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x8c3478d244e73a1e,
            0x9ad80eaea2847ef,
            0xa009d2c6a19c4d7c,
            0xad14ea37d3caa2e6,
            0xf4d4bfd1cd09cd5e,
            0x633ac7191e6e1cd,
        ]),
        c1: Fp::from_raw_unchecked([
            0xa242b42d478776e5,
            0x3474f7bb939b1fde,
            0xca9317dcdf327fbb,
            0x92c8b4def6629d23,
            0xf3db74fd5208df9e,
            0xc3ff42813809dbe,
        ]),
    };
    let zi = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x396e3a171d3682eb,
            0x9d8a8a66679bed76,
            0xd149d6008a42ad3c,
            0x642a4f268fc07724,
            0xfe94535d55e01ead,
            0x17ee1dc7f35478ab,
        ]),
        c1: Fp::from_raw_unchecked([
            0xaec2bbaef501dcc,
            0xf647afeae70603a1,
            0x268aeee6dc5d1d38,
            0xc9a8fe231d8fec6b,
            0xbbfb6687b11e9507,
            0x8f1b25b6b4c4d0e,
        ]),
    };
    let mut pt = G2Projective {
        x: xi,
        y: yi,
        z: zi,
    };
    pt.isogeny_map();

    let x_expect = Fp2 {
        c0: Fp::from_raw_unchecked([
            0xaa51c0e1180b0d57,
            0x28b0e686761afc8c,
            0x19ff407a3a484438,
            0x63ea8ff9abf6fddf,
            0x5c6aa531d2636d17,
            0x185badd0900f1073,
        ]),
        c1: Fp::from_raw_unchecked([
            0x6ac4e087d0ba1981,
            0x37118ac6f0aa47f6,
            0x1ee51ba02ed89f1f,
            0xb7d4b12c8096b32,
            0x27bb5e93c6038b03,
            0x184c84a00727706,
        ]),
    };
    let y_expect = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x8c76ec95e02f5d6a,
            0x5d31fde09b3d4d4e,
            0xa4e17e369860ac14,
            0xda389b35e29b699a,
            0xe67fa8059a237afa,
            0x1220f324257645f7,
        ]),
        c1: Fp::from_raw_unchecked([
            0x2fdab2cdbc61a883,
            0x5dcf5ada26fc9714,
            0x995cd083b242fec6,
            0xa0079db2fc99f60b,
            0xe01e053554e727bb,
            0x55a9d8217bab54f,
        ]),
    };
    let z_expect = Fp2 {
        c0: Fp::from_raw_unchecked([
            0x74daa669f740f9f2,
            0x283671a4e6a1c8f,
            0x14230a1fecab7ad9,
            0xb34dce82bd2eb9dd,
            0x866e706af63aef04,
            0x252ef31770c147d,
        ]),
        c1: Fp::from_raw_unchecked([
            0xc1b8083edd204658,
            0xa9c56c859e3d920c,
            0x1e0b5d8a59b68688,
            0x60d10c39c8c72e4e,
            0x948f508d12a11262,
            0xb7e31697882943e,
        ]),
    };
    assert_eq!(pt.x, x_expect);
    assert_eq!(pt.y, y_expect);
    assert_eq!(pt.z, z_expect);
}

#[cfg(test)]
fn check_g2_prime(x: &Fp2, y: &Fp2, z: &Fp2) -> subtle::Choice {
    // (X : Y : Z)==(X/Z, Y/Z) is on E: y^2 = x^3 + ELLP_A * x + ELLP_B.
    // y^2 z = (x^3) + A (x z^2) + B z^3
    let zsq = z.square();
    (y.square() * z).ct_eq(&(x.square() * x + ELLP_A * x * zsq + ELLP_B * zsq * z))
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
    use crate::g2::G2Affine;
    use crate::hash_to_curve::{ExpandMsgXmd, HashToCurve};

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

// #[test]
// fn test_hash_to_curve_07() {
//     use crate::{ExpandMsgXmd, HashToCurve};

//     struct TestCase {
//         msg: &'static [u8],
//         expected: [&'static str; 4],
//     }
//     impl TestCase {
//         fn expected(&self) -> String {
//             self.expected[0].to_string() + self.expected[1] + self.expected[2] + self.expected[3]
//         }
//     }

//     const DOMAIN: &[u8] = b"BLS12381G2_XMD:SHA-256_SSWU_RO_TESTGEN";

//     let cases = vec![
//         TestCase {
//             msg: b"",
//             expected: [
//         "0fbdae26f9f9586a46d4b0b70390d09064ef2afe5c99348438a3c7d9756471e015cb534204c1b6824617a85024c772dc",
//         "0a650bd36ae7455cb3fe5d8bb1310594551456f5c6593aec9ee0c03d2f6cb693bd2c5e99d4e23cbaec767609314f51d3",
//         "02e5cf8f9b7348428cc9e66b9a9b36fe45ba0b0a146290c3a68d92895b1af0e1f2d9f889fb412670ae8478d8abd4c5aa",
//         "0d8d49e7737d8f9fc5cef7c4b8817633103faf2613016cb86a1f3fc29968fe2413e232d9208d2d74a89bf7a48ac36f83",
//             ],
//         },
//         TestCase {
//             msg: b"abc",
//             expected: [
//         "03578447618463deb106b60e609c6f7cc446dc6035f84a72801ba17c94cd800583b493b948eff0033f09086fdd7f6175",
//         "1953ce6d4267939c7360756d9cca8eb34aac4633ef35369a7dc249445069888e7d1b3f9d2e75fbd468fbcbba7110ea02",
//         "0184d26779ae9d4670aca9b267dbd4d3b30443ad05b8546d36a195686e1ccc3a59194aea05ed5bce7c3144a29ec047c4",
//         "0882ab045b8fe4d7d557ebb59a63a35ac9f3d312581b509af0f8eaa2960cbc5e1e36bb969b6e22980b5cbdd0787fcf4e",
//             ],
//         },
//         TestCase {
//             msg: b"abcdef0123456789",
//             expected: [
//         "195fad48982e186ce3c5c82133aefc9b26d55979b6f530992a8849d4263ec5d57f7a181553c8799bcc83da44847bdc8d",
//         "17b461fc3b96a30c2408958cbfa5f5927b6063a8ad199d5ebf2d7cdeffa9c20c85487204804fab53f950b2f87db365aa",
//         "005cdf3d984e3391e7e969276fb4bc02323c5924a4449af167030d855acc2600cf3d4fab025432c6d868c79571a95bef",
//         "174a3473a3af2d0302b9065e895ca4adba4ece6ce0b41148ba597001abb152f852dd9a96fb45c9de0a43d944746f833e",
//             ]
//         },
//         TestCase {
//             msg: b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
//                    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
//             expected: [
//         "123b6bd9feeba26dd4ad00f8bfda2718c9700dc093ea5287d7711844644eb981848316d3f3f57d5d3a652c6cdc816aca",
//         "0a162306f3b0f2bb326f0c4fb0e1fea020019c3af796dcd1d7264f50ddae94cacf3cade74603834d44b9ab3d5d0a6c98",
//         "05483f3b96d9252dd4fc0868344dfaf3c9d145e3387db23fa8e449304fab6a7b6ec9c15f05c0a1ea66ff0efcc03e001a",
//         "15c1d4f1a685bb63ee67ca1fd96155e3d091e852a684b78d085fd34f6091e5249ddddbdcf2e7ec82ce6c04c63647eeb7",
//             ]
//         }
//     ];

//     for case in cases {
//         let g = <G2 as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(&case.msg, DOMAIN);
//         let g_uncompressed = g.into_affine().into_uncompressed();

//         assert_eq!(case.expected(), hex::encode(&g_uncompressed.0[..]));
//     }
// }
