#[macro_use]
extern crate criterion;

extern crate bls12_381;
use bls12_381::*;

use criterion::{black_box, Criterion};
use subtle::ConditionallySelectable;

fn criterion_benchmark(c: &mut Criterion) {
    const BATCH: usize = 1024;
    let x = Scalar::from_raw([1, 2, 3, 4]);
    let y = Scalar::from_raw([5, 6, 7, 8]);
    let narrow = b"01234567801234567801234567801234";
    let wide = b"0123456780123456780123456780123456780123456780123456780123456780";

    c.bench_function("Scalar add", move |b| {
        b.iter(|| black_box(x) + black_box(y))
    });
    c.bench_function("Scalar sub", move |b| {
        b.iter(|| black_box(x) - black_box(y))
    });
    c.bench_function("Scalar double", move |b| b.iter(|| black_box(x).double()));
    c.bench_function("Scalar negate", move |b| b.iter(|| -black_box(x)));
    c.bench_function("Scalar mul", move |b| {
        b.iter(|| black_box(x) * black_box(y))
    });
    c.bench_function("Scalar square", move |b| b.iter(|| black_box(x).square()));
    c.bench_function("Scalar sqrt", move |b| b.iter(|| black_box(x).sqrt()));
    c.bench_function("Scalar invert", move |b| b.iter(|| black_box(x).invert()));
    c.bench_function("Scalar invert vartime", move |b| {
        b.iter(|| black_box(x).invert_vartime())
    });
    c.bench_function(&format!("Scalar batch invert n={}", BATCH), {
        let batch = vec![x; BATCH];
        let mut out = vec![Scalar::zero(); BATCH];
        move |b| b.iter(|| Scalar::batch_invert(black_box(&batch), black_box(&mut out)))
    });
    c.bench_function("Scalar pow", move |b| {
        b.iter(|| black_box(x).pow(&[1111111, 22222222, 33333333, 44444444]))
    });
    c.bench_function("Scalar pow_vartime", move |b| {
        b.iter(|| black_box(x).pow_vartime(&[1111111, 22222222, 33333333, 44444444]))
    });
    c.bench_function("Scalar from_bytes", move |b| {
        b.iter(|| Scalar::from_bytes(black_box(narrow)))
    });
    c.bench_function("Scalar from_bytes_vartime", move |b| {
        b.iter(|| Scalar::from_bytes_vartime(black_box(narrow)))
    });
    c.bench_function("Scalar from_bytes_wide", move |b| {
        b.iter(|| Scalar::from_bytes_wide(black_box(&wide)))
    });
    c.bench_function("Scalar to_bytes", move |b| {
        b.iter(|| black_box(x).to_bytes())
    });
    c.bench_function("Scalar conditional_select", move |b| {
        use subtle::Choice;
        b.iter(|| Scalar::conditional_select(&x, &y, Choice::from(0)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
