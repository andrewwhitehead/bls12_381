#[macro_use]
extern crate criterion;

extern crate bls12_381;
use bls12_381::*;

use criterion::{black_box, Criterion};

fn criterion_benchmark(c: &mut Criterion) {
    // G1Affine
    {
        let name = "G1Affine";
        let a = G1Affine::generator();
        let s = Scalar::from_raw([1, 2, 3, 4]);
        let compressed = [0u8; 48];
        let uncompressed = [0u8; 96];

        use rand_core::SeedableRng;
        let mut rng = rand_xorshift::XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        const M: usize = 256;
        let v = vec![G1Affine::generator(); M];
        let ss: Vec<Scalar> = (0..M).map(|_| Scalar::random(&mut rng)).collect();

        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} doubling", name), move |b| {
            b.iter(|| black_box(a).double())
        });
        c.bench_function(&format!("{} addition", name), move |b| {
            b.iter(|| black_box(a).add(&a))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} sum of products n={}", name, M), |b| {
            b.iter(|| G1Affine::sum_of_products(black_box(&v[..M]), black_box(&ss[..M])))
        });
        c.bench_function(&format!("{} sum of products vartime n={}", name, M), |b| {
            b.iter(|| G1Affine::sum_of_products_vartime(black_box(&v[..M]), black_box(&ss[..M])))
        });
        c.bench_function(&format!("{} subgroup check", name), move |b| {
            b.iter(|| black_box(a).is_torsion_free())
        });
        c.bench_function(
            &format!("{} deserialize compressed point", name),
            move |b| b.iter(|| G1Affine::from_compressed(black_box(&compressed))),
        );
        c.bench_function(
            &format!("{} deserialize uncompressed point", name),
            move |b| b.iter(|| G1Affine::from_uncompressed(black_box(&uncompressed))),
        );
    }

    // G1Projective
    {
        let name = "G1Projective";
        let a = G1Projective::generator();
        let a_affine = G1Affine::generator();
        let s = Scalar::from_raw([1, 2, 3, 4]);

        use rand_core::SeedableRng;
        let mut rng = rand_xorshift::XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        const N: usize = 10000;
        const M: usize = 256;
        let v = vec![G1Projective::generator(); N];
        let mut q = vec![G1Affine::identity(); N];
        let ss: Vec<Scalar> = (0..M).map(|_| Scalar::random(&mut rng)).collect();

        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} to affine", name), move |b| {
            b.iter(|| G1Affine::from(black_box(a)))
        });
        c.bench_function(&format!("{} doubling", name), move |b| {
            b.iter(|| black_box(a).double())
        });
        c.bench_function(&format!("{} addition", name), move |b| {
            b.iter(|| black_box(a).add(&a))
        });
        c.bench_function(&format!("{} mixed addition", name), move |b| {
            b.iter(|| black_box(a).add_mixed(&a_affine))
        });
        c.bench_function(&format!("{} scalar multiplication", name), |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} sum of products n={}", name, M), |b| {
            b.iter(|| G1Projective::sum_of_products(black_box(&v[..M]), black_box(&ss[..M])))
        });
        c.bench_function(&format!("{} sum of products vartime n={}", name, M), |b| {
            b.iter(|| {
                G1Projective::sum_of_products_vartime(black_box(&v[..M]), black_box(&ss[..M]))
            })
        });
        c.bench_function(&format!("{} batch to affine n={}", name, N), move |b| {
            b.iter(|| {
                G1Projective::batch_normalize(black_box(&v), black_box(&mut q));
                black_box(&q)[0]
            })
        });
    }

    // G2Affine
    {
        let name = "G2Affine";
        let a = G2Affine::generator();
        let s = Scalar::from_raw([1, 2, 3, 4]);
        let compressed = [0u8; 96];
        let uncompressed = [0u8; 192];

        use rand_core::SeedableRng;
        let mut rng = rand_xorshift::XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        const M: usize = 256;
        let v = vec![G2Affine::generator(); M];
        let ss: Vec<Scalar> = (0..M).map(|_| Scalar::random(&mut rng)).collect();

        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} sum of products n={}", name, M), |b| {
            b.iter(|| G2Affine::sum_of_products(black_box(&v[..M]), black_box(&ss[..M])))
        });
        c.bench_function(&format!("{} sum of products vartime n={}", name, M), |b| {
            b.iter(|| G2Affine::sum_of_products_vartime(black_box(&v[..M]), black_box(&ss[..M])))
        });
        c.bench_function(&format!("{} subgroup check", name), move |b| {
            b.iter(|| black_box(a).is_torsion_free())
        });
        c.bench_function(
            &format!("{} deserialize compressed point", name),
            move |b| b.iter(|| G2Affine::from_compressed(black_box(&compressed))),
        );
        c.bench_function(
            &format!("{} deserialize uncompressed point", name),
            move |b| b.iter(|| G2Affine::from_uncompressed(black_box(&uncompressed))),
        );
    }

    // G2Projective
    {
        let name = "G2Projective";
        let a = G2Projective::generator();
        let a_affine = G2Affine::generator();
        let s = Scalar::from_raw([1, 2, 3, 4]);

        use rand_core::SeedableRng;
        let mut rng = rand_xorshift::XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        const N: usize = 10000;
        const M: usize = 256;
        let v = vec![G2Projective::generator(); N];
        let mut q = vec![G2Affine::identity(); N];
        let ss: Vec<Scalar> = (0..M).map(|_| Scalar::random(&mut rng)).collect();

        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("{} to affine", name), move |b| {
            b.iter(|| G2Affine::from(black_box(a)))
        });
        c.bench_function(&format!("{} doubling", name), move |b| {
            b.iter(|| black_box(a).double())
        });
        c.bench_function(&format!("{} addition", name), move |b| {
            b.iter(|| black_box(a).add(&a))
        });
        c.bench_function(&format!("{} mixed addition", name), move |b| {
            b.iter(|| black_box(a).add_mixed(&a_affine))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("{} sum of products n={}", name, M), |b| {
            b.iter(|| G2Projective::sum_of_products(black_box(&v[..M]), black_box(&ss[..M])))
        });
        c.bench_function(&format!("{} sum of products vartime n={}", name, M), |b| {
            b.iter(|| {
                G2Projective::sum_of_products_vartime(black_box(&v[..M]), black_box(&ss[..M]))
            })
        });
        c.bench_function(&format!("{} batch to affine n={}", name, N), move |b| {
            b.iter(|| {
                G2Projective::batch_normalize(black_box(&v), black_box(&mut q));
                black_box(&q)[0]
            })
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
