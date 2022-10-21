#[macro_use]
extern crate criterion;

extern crate bls12_381;
use bls12_381::*;

use criterion::{black_box, Criterion};

fn criterion_benchmark(c: &mut Criterion) {
    // Pairings
    {
        let g = G1Affine::generator();
        let h = G2Affine::generator();
        c.bench_function("full pairing", move |b| {
            b.iter(|| pairing(black_box(&g), black_box(&h)))
        });
        c.bench_function("G2 preparation for pairing", move |b| {
            b.iter(|| G2Prepared::from(h))
        });
        let prep = G2Prepared::from(h);
        c.bench_function("miller loop for pairing", move |b| {
            b.iter(|| multi_miller_loop(&[(&g, &prep)]))
        });
        let prep = G2Prepared::from(h);
        let r = multi_miller_loop(&[(&g, &prep)]);
        c.bench_function("final exponentiation for pairing", move |b| {
            b.iter(|| r.final_exponentiation())
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
