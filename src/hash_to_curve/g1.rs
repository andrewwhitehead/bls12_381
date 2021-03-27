use subtle::{ConditionallyNegatable, ConditionallySelectable, ConstantTimeEq};

use super::chain::chain_pm3div4;
use super::MessageToField;
use crate::fp::Fp;
use crate::g1::G1Projective;
use crate::generic_array::{typenum::U64, GenericArray};

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

const ELLP_A: Fp = Fp::from_raw_unchecked([
    0x2f65aa0e9af5aa51,
    0x86464c2d1e8416c3,
    0xb85ce591b7bd31e2,
    0x27e11c91b5f24e7c,
    0x28376eda6bfc1835,
    0x155455c3e5071d85,
]);

const ELLP_B: Fp = Fp::from_raw_unchecked([
    0xfb996971fe22a1e0,
    0x9aa93eb35b742d6f,
    0x8c476013de99c5c4,
    0x873e27c3a221e571,
    0xca72b5e45a52d888,
    0x06824061418a386b,
]);

const XI: Fp = Fp::from_raw_unchecked([
    0x886c00000023ffdc,
    0x0f70008d3090001d,
    0x77672417ed5828c3,
    0x9dac23e943dc1740,
    0x50553f1b9c131521,
    0x078c712fbe0ab6e8,
]);

const SQRT_M_XI_CUBED: Fp = Fp::from_raw_unchecked([
    0x43b571cad3215f1f,
    0xccb460ef1c702dc2,
    0x742d884f4f97100b,
    0xdb2c3e3238a3382b,
    0xe40f3fa13fce8f88,
    0x0073a2af9892a2ff,
]);

impl MessageToField for G1Projective {
    type InputLength = U64;
    type Pt = Fp;

    fn input_okm(okm: &GenericArray<u8, U64>) -> Fp {
        const F_2_256: Fp = Fp::from_raw_unchecked([
            0x75b3cd7c5ce820f,
            0x3ec6ba621c3edb0b,
            0x168a13d82bff6bce,
            0x87663c4bf8c449d2,
            0x15f34c83ddc8d830,
            0xf9628b49caa2e85,
        ]);

        let mut bs = [0u8; 48];
        bs[16..].copy_from_slice(&okm[..32]);
        let db = Fp::from_bytes(&bs).unwrap();
        // println!("db: {:?} {:#?}", db, db);
        bs[16..].copy_from_slice(&okm[32..]);
        let da = Fp::from_bytes(&bs).unwrap();

        db * F_2_256 + da
    }

    fn isogeny_map(&mut self) {
        *self = self.map_isogeny([&XNUM[..], &XDEN[..], &YNUM[..], &YDEN[..]]);
    }

    fn osswu_map(u: &Fp) -> G1Projective {
        // compute x0 and g(x0)
        // let [usq, xi_usq, _, x0_num, x0_den, gx0_num, gx0_den] =
        //     osswu_helper!(Fp, u, &XI, &ELLP_A, &ELLP_B);

        let tmp1 = u.square();
        let tmp3 = tmp1 * XI;
        let mut tmp2 = tmp3.square();
        let mut xd = tmp3 + tmp2;
        let x1n = (xd + Fp::one()) * ELLP_B;
        xd *= ELLP_A;
        xd = xd.neg();
        xd.conditional_assign(&(XI * ELLP_A), xd.is_zero()); // If xd == 0, set xd = Z * A
        tmp2 = xd.square();
        let gxd = tmp2 * xd; // xd^3
        tmp2 *= ELLP_A;
        let gx1 = (x1n.square() + tmp2) * x1n + ELLP_B * gxd; // x1n^3 + A * x1n * xd^2 + B * xd^3
        tmp2 = gx1 * gxd;
        let tmp4 = tmp2 * gxd.square(); // gx1 * gxd^3
        let mut y1 = chain_pm3div4(&tmp4); // (gx1 * gxd^3)^((p - 3) / 4)
        y1 *= tmp2; // gx1 * gxd * (gx1 * gxd^3)^((p - 3) / 4)
        let x2n = tmp3 * x1n; // x2 = x2n / xd = Z * u^2 * x1n / xd
        let y2 = y1 * SQRT_M_XI_CUBED * tmp1 * u;
        tmp2 = y1.square() * gxd;
        let e2 = tmp2.ct_eq(&gx1);
        let xn = Fp::conditional_select(&x2n, &x1n, e2);
        let mut y = Fp::conditional_select(&y2, &y1, e2);
        let e3 = y.sgn0() ^ u.sgn0();
        y.conditional_negate(e3);

        G1Projective {
            x: xn,
            y: y * xd,
            z: xd,
        }

        // // compute g(X0(u)) ^ ((p - 3) // 4)
        // let sqrt_candidate = {
        //     let mut tmp1 = gx0_num;
        //     tmp1.mul_assign(&gx0_den); // u * v
        //     let mut tmp2 = gx0_den.square(); // v^2
        //     tmp2.mul_assign(&tmp1); // u * v^3
        //     let tmp3 = tmp2;
        //     chain_pm3div4(&mut tmp2, &tmp3); // (u v^3) ^ ((p - 3) // 4)
        //     tmp2.mul_assign(&tmp1); // u v (u v^3) ^ ((p - 3) // 4)
        //     tmp2
        // };

        // // select correct values for y and for x numerator
        // let (mut x_num, mut y) = {
        //     let mut test_cand = sqrt_candidate.square();
        //     test_cand.mul_assign(&gx0_den);
        //     let mut x1_num = x0_num;
        //     x1_num.mul_assign(&xi_usq); // x1 = xi u^2 g(x0)
        //     let mut y1 = usq; // y1 = sqrt(-xi**3) * u^3 g(x0) ^ ((p - 1) // 4)
        //     y1.mul_assign(u);
        //     y1.mul_assign(&sqrt_candidate);
        //     y1.mul_assign(&SQRT_M_XI_CUBED);
        //     let gx0_square = test_cand.ct_eq(&gx0_num); // g(x0) is square
        //     (
        //         Fp::conditional_select(&x1_num, &x0_num, gx0_square),
        //         Fp::conditional_select(&y1, &sqrt_candidate, gx0_square),
        //     )
        // };

        // // make sure sign of y and sign of u agree
        // let sgn0_y_xor_u = y.sgn0() ^ u.sgn0();
        // y.conditional_negate(sgn0_y_xor_u);

        // G1Projective {
        //     x: x_num,
        //     y: y * x0_den,
        //     z: x0_den,
        // }
    }

    fn clear_h(&mut self) {
        *self = self.clear_cofactor();
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
    println!("from {}", crate::g1::G1Affine::from(pt));
    println!(
        "from {}",
        hex::encode(&crate::g1::G1Affine::from(pt).to_compressed()[..])
    );

    pt.isogeny_map();
    assert!(bool::from(pt.is_on_curve()));
    let xo = Fp::from_raw_unchecked([
        0xb129fab9bef88edd,
        0x1c5429e2f4b8bc35,
        0xcaab8cc9ec4893f2,
        0x9e9c31f30a607c8b,
        0x9661fcf22bedfddb,
        0x10fc4a3ba5f48e07,
    ]);
    let yo = Fp::from_raw_unchecked([
        0xaf52c5fbd490f370u64,
        0x1533c0f27b46c02fu64,
        0xc8890dd0987b134fu64,
        0x43e2d5f172257d50u64,
        0x538ebef63fb145beu64,
        0x11eab1145b95cb9fu64,
    ]);
    let zo = Fp::from_raw_unchecked([
        0x7441c43513e11f49u64,
        0x620b0af2483ad30fu64,
        0x678c5bf3ad4090b4u64,
        0xc75152c6f387d070u64,
        0x5f3cc0ed1bd3f0eeu64,
        0x12514e630a486abbu64,
    ]);

    let xo = xo * zo; // * zo.invert().unwrap();
    let yo = yo * zo; // * zo; //zo.invert().unwrap() * zo.invert().unwrap();
    println!(
        "expect {:#x?}",
        G1Projective {
            x: xo,
            y: yo,
            z: zo
        }
    );
    println!(
        "?? {}",
        crate::g1::G1Affine::from(G1Projective {
            x: xo,
            y: yo,
            z: zo
        })
    );
    println!(
        "?? {}",
        hex::encode(
            &crate::g1::G1Affine::from(G1Projective {
                x: xo,
                y: yo,
                z: zo
            })
            .to_compressed()[..]
        )
    );

    println!("got {:#x?}", pt);

    assert_eq!(&pt.x, &xo);
    assert_eq!(&pt.y, &yo);
    assert_eq!(&pt.z, &zo);
}

#[test]
fn test_projective_iso11_fixed() {
    let xi = Fp::from_raw_unchecked([
        0xf6adc4118ae592abu64,
        0xa384a7ab165def35u64,
        0x2365b1fb1c8a73bfu64,
        0xc40dc338ca285231u64,
        0x47ff3364428c59b3u64,
        0x1789051238d025e3u64,
    ]);
    let yi = Fp::from_raw_unchecked([
        0x1a635634e9cced27u64,
        0x03f604e47bc51aa9u64,
        0x06f6ff472fa7276eu64,
        0x0459ed10f1f8abb1u64,
        0x8e76c82bd4a29d21u64,
        0x088cb5712bf81924u64,
    ]);
    let zi = Fp::from_raw_unchecked([
        0x0416411fe2e97d06u64,
        0xaced7fec7a63fe65u64,
        0x683295bcaed54202u64,
        0xbdc3405df9ff0a3bu64,
        0xf9698f57510273fbu64,
        0x064bb4b501466b2au64,
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

#[cfg(test)]
fn check_g1_prime(x: &Fp, y: &Fp, z: &Fp) -> subtle::Choice {
    // (X : Y : Z)==(X/Z, Y/Z) is on E: y^2 = x^3 + ELLP_A * x + ELLP_B.
    // y^2 z = (x^3) + A (x z^2) + B z^3
    let zsq = z.square();
    (y.square() * z).ct_eq(&(x.square() * x + ELLP_A * x * zsq + ELLP_B * zsq * z))
}

#[test]
fn test_osswu_expected() {
    // exceptional case: zero
    let p = G1Projective::osswu_map(&Fp::zero());
    let G1Projective { x, y, z } = &p;
    let xo = Fp::from_raw_unchecked([
        0xfb996971fe22a1e0,
        0x9aa93eb35b742d6f,
        0x8c476013de99c5c4,
        0x873e27c3a221e571,
        0xca72b5e45a52d888,
        0x6824061418a386b,
    ]);
    let yo = Fp::from_raw_unchecked([
        0xfd6fced87a7f11a3,
        0x9a6b314b03c8db31,
        0x41f85416e0eab593,
        0xfeeb089f7e6ec4d7,
        0x85a134c37ed1278f,
        0x575c525bb9f74bb,
    ]);
    let zo = Fp::from_raw_unchecked([
        0x7f674ea0a8915178,
        0xb0f945fc13b8fa65,
        0x4b46759a38e87d76,
        0x2e7a929641bbb6a1,
        0x1668ddfa462bf6b6,
        0x960e2ed1cf294c,
    ]);
    println!("proj {:?}", crate::g1::G1Affine::from(p));
    println!("proj {:#?}", p);
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    assert!(bool::from(check_g1_prime(x, y, z)));

    // exceptional case: sqrt(-1/XI) (positive)
    let excp = Fp::from_raw_unchecked([
        0xf3d0477e91edbf,
        0x8d6621e4ca8dc69,
        0xb9cf7927b19b9726,
        0xba133c996cafa2ec,
        0xed2a5ccd5ca7bb68,
        0x19cb022f8ee9d73b,
    ]);
    let p = G1Projective::osswu_map(&excp);
    let G1Projective { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    assert!(bool::from(check_g1_prime(x, y, z)));

    // exceptional case: sqrt(-1/XI) (negative)
    let excp = Fp::from_raw_unchecked([
        0xb90b2fb8816dbcec,
        0x15d59de064ab2396,
        0xad61597945155efe,
        0xaa640eeb86d56fd2,
        0x5df14ae8e6a3f16e,
        0x360fbaaa960f5e,
    ]);
    let p = G1Projective::osswu_map(&excp);
    let G1Projective { x, y, z } = &p;
    let myo = -yo;
    assert_eq!(x, &xo);
    assert_eq!(y, &myo);
    assert_eq!(z, &zo);
    assert!(bool::from(check_g1_prime(x, y, z)));

    let u = Fp::from_raw_unchecked([
        0xa618fa19f7e2eadc,
        0x93c7f1fc876ba245,
        0xe2ed4cc47b5c0ae0,
        0xd49efa74e4a8d000,
        0xa0b23ba692b5431c,
        0xd1551f2d7d8d193,
    ]);
    let xo = Fp::from_raw_unchecked([
        0x2197ca55fab3ba48,
        0x591deb39f434949a,
        0xf9df7fb4f1fa6a08,
        0x59e3c16a9dfa8fa5,
        0xe5929b194aad5f7a,
        0x130a46a4c61b44ed,
    ]);
    let yo = Fp::from_raw_unchecked([
        0xf7215b58c7200ad0,
        0x890516313a4e66bf,
        0xc9031acc8a3619a8,
        0xea1f9978fde3ffec,
        0x548f02d6cfbf472,
        0x169375573529163f,
    ]);
    let zo = Fp::from_raw_unchecked([
        0xf36feb2e1128ade0,
        0x42e22214250bcd94,
        0xb94f6ba2dddf62d6,
        0xf56d4392782bf0a2,
        0xb2d7ce1ec26309e7,
        0x182b57ed6b99f0a1,
    ]);
    let p = G1Projective::osswu_map(&u);
    let G1Projective { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    assert!(bool::from(check_g1_prime(x, y, z)));
}

#[test]
fn test_osswu_semirandom() {
    use rand_core::SeedableRng;
    let mut rng = rand_xorshift::XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    for _ in 0..32 {
        let input = Fp::random(&mut rng);
        let p = G1Projective::osswu_map(&input);
        let G1Projective { x, y, z } = &p;
        assert!(bool::from(check_g1_prime(x, y, z)));
    }
}

#[cfg(feature = "alloc")]
#[test]
fn test_encode_to_curve_07() {
    use alloc::string::{String, ToString};

    use crate::g1::G1Affine;
    use crate::hash_to_curve::{ExpandMsgXmd, HashToCurve};

    struct TestCase {
        msg: &'static [u8],
        expected: [&'static str; 2],
    }
    impl TestCase {
        fn expected(&self) -> String {
            self.expected[0].to_string() + self.expected[1]
        }
    }

    const DOMAIN: &[u8] = b"BLS12381G1_XMD:SHA-256_SSWU_NU_TESTGEN";

    let cases = vec![
        TestCase {
            msg: b"",
            expected: [
        "1223effdbb2d38152495a864d78eee14cb0992d89a241707abb03819a91a6d2fd65854ab9a69e9aacb0cbebfd490732c",
        "0f925d61e0b235ecd945cbf0309291878df0d06e5d80d6b84aa4ff3e00633b26f9a7cb3523ef737d90e6d71e8b98b2d5",
            ],
        },
        TestCase {
            msg: b"abc",
            expected: [
        "179d3fd0b4fb1da43aad06cea1fb3f828806ddb1b1fa9424b1e3944dfdbab6e763c42636404017da03099af0dcca0fd6",
        "0d037cb1c6d495c0f5f22b061d23f1be3d7fe64d3c6820cfcd99b6b36fa69f7b4c1f4addba2ae7aa46fb25901ab483e4",

            ],
        },
        TestCase {
            msg: b"abcdef0123456789",
            expected: [
        "15aa66c77eded1209db694e8b1ba49daf8b686733afaa7b68c683d0b01788dfb0617a2e2d04c0856db4981921d3004af",
        "0952bb2f61739dd1d201dd0a79d74cda3285403d47655ee886afe860593a8a4e51c5b77a22d2133e3a4280eaaaa8b788",
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
        "06328ce5106e837935e8da84bd9af473422e62492930aa5f460369baad9545defa468d9399854c23a75495d2a80487ee",
        "094bfdfe3e552447433b5a00967498a3f1314b86ce7a7164c8a8f4131f99333b30a574607e301d5f774172c627fd0bca",
            ]
        }
    ];

    for case in cases {
        let g = <G1Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::encode_to_curve(
            &case.msg, DOMAIN,
        );
        let aff = G1Affine::from(g);
        println!("aff {}", hex::encode(&aff.to_compressed()[..]));
        let g_uncompressed = aff.to_uncompressed();

        assert_eq!(case.expected(), hex::encode(&g_uncompressed[..]));
    }
}

#[cfg(feature = "alloc")]
#[test]
fn test_hash_to_curve_07() {
    use alloc::string::{String, ToString};

    use crate::g1::G1Affine;
    use crate::hash_to_curve::{ExpandMsgXmd, HashToCurve};

    struct TestCase {
        msg: &'static [u8],
        expected: [&'static str; 2],
    }
    impl TestCase {
        fn expected(&self) -> String {
            self.expected[0].to_string() + self.expected[1]
        }
    }

    const DOMAIN: &[u8] = b"BLS12381G1_XMD:SHA-256_SSWU_RO_TESTGEN";

    let cases = vec![
        TestCase {
            msg: b"",
            expected: [
                "0576730ab036cbac1d95b38dca905586f28d0a59048db4e8778782d89bff856ddef89277ead5a21e2975c4a6e3d8c79e",
                "1273e568bebf1864393c517f999b87c1eaa1b8432f95aea8160cd981b5b05d8cd4a7cf00103b6ef87f728e4b547dd7ae",
            ],
        },
        TestCase {
            msg: b"abc",
            expected: [
                "061daf0cc00d8912dac1d4cf5a7c32fca97f8b3bf3f805121888e5eb89f77f9a9f406569027ac6d0e61b1229f42c43d6",
                "0de1601e5ba02cb637c1d35266f5700acee9850796dc88e860d022d7b9e7e3dce5950952e97861e5bb16d215c87f030d"
            ],
        },
        TestCase {
            msg: b"abcdef0123456789",
            expected: [
                "0fb3455436843e76079c7cf3dfef75e5a104dfe257a29a850c145568d500ad31ccfe79be9ae0ea31a722548070cf98cd",
                "177989f7e2c751658df1b26943ee829d3ebcf131d8f805571712f3a7527ee5334ecff8a97fc2a50cea86f5e6212e9a57"
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
                "0514af2137c1ae1d78d5cb97ee606ea142824c199f0f25ac463a0c78200de57640d34686521d3e9cf6b3721834f8a038",
                "047a85d6898416a0899e26219bca7c4f0fa682717199de196b02b95eaf9fb55456ac3b810e78571a1b7f5692b7c58ab6"
            ]
        }
    ];

    for case in cases {
        let g = <G1Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(
            &case.msg, DOMAIN,
        );
        let g_uncompressed = G1Affine::from(g).to_uncompressed();

        assert_eq!(case.expected(), hex::encode(&g_uncompressed[..]));
    }
}