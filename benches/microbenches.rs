use std::array;

use ark_ff::{AdditiveGroup, Field};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use polynomial_proving::{interpolate_polynomial, measure_multiplications_for_exclusion, pad_walk, polynomials::{
    EqFixedPoint, HypercubeEvalPoly, HypercubePoint, PiopPolynomial, piop_polynomials::{Mask, P, Q, QConstructionMode}
}};
use util::algebra::field::{arkfield::Fp2256, p434, p503, p610, p751};

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut assignment = [Fp2256::ZERO; 8];
    for (i, a) in assignment.iter_mut().enumerate() {
        *a = Fp2256::from(i as u64);
    }
    let q = Q::<{ 1 << 8 }, _>::new(&assignment);

    let mut point = [Fp2256::ZERO; 2 * 8 + 1];
    for (i, a) in point.iter_mut().enumerate() {
        *a = Fp2256::from(i as u64);
    }

    c.bench_function("Q eval", |b| {
        b.iter(|| q.eval(black_box(&point)));
    });
    let f = Fp2256::from(3423423);
    let hypercube_point = HypercubePoint::from_usize(0b1101001_1_10110011);
    c.bench_function("Q eval field then hypercube", |b| {
        b.iter(|| q.eval_field_then_hypercube(black_box(f), black_box(hypercube_point)));
    });

    let eq = EqFixedPoint::new(HypercubePoint::from_usize(0b1011_0111), 8, Fp2256::ONE);
    c.bench_function("eq eval", |b| {
        b.iter(|| eq.eval(black_box(&point)));
    });

    let vals: [Fp2256; 1 << 9] = array::from_fn(|i| Fp2256::from(i as u32));
    let multilinear_extension = HypercubeEvalPoly::new(&vals, 9);
    c.bench_function("multilinear extension eval", |b| {
        b.iter(|| multilinear_extension.eval(black_box(&point)));
    });
    let point = Fp2256::from(2354234);
    let hypercube = HypercubePoint::from_usize(0b1101_0101);
    c.bench_function("multilinear extension eval field then hypercube", |b| {
        b.iter(|| {
            multilinear_extension.eval_field_then_hypercube(black_box(point), black_box(hypercube))
        });
    });

    let e = [Fp2256::from(2), Fp2256::from(3)];
    let k = [Fp2256::from(25235); 8];
    let b_evals_0 = [Fp2256::from(32423); 2 * (1 << 8)];
    let b_evals_1 = [Fp2256::from(32423); 2 * (1 << 8)];
    let poly_0 = HypercubeEvalPoly::new(&b_evals_0, 1 + 8);
    let poly_1 = HypercubeEvalPoly::new(&b_evals_1, 1 + 8);
    c.bench_function("p new", |b| {
        b.iter(|| {
            P::<17, 9, 8, { 1 << 8 }, { 1 << 9 }, 16, _>::new(
                black_box(e),
                black_box(k),
                black_box(&poly_0),
                black_box(&poly_1),
                QConstructionMode::GrayCodes,
            )
        });
    });

    let p: P<17, 9, 8, { 1 << 8 }, { 2 * (1 << 8) }, 16, Fp2256> = P::new(
        black_box(e),
        black_box(k),
        black_box(&poly_0),
        black_box(&poly_1),
        QConstructionMode::GrayCodes,
    );
    let point = [Fp2256::from(52348); 17];
    c.bench_function("p eval", |b| {
        b.iter(|| p.eval(black_box(&point)));
    });

    let point = HypercubePoint::from_usize(0b10111101_1_10110101);
    c.bench_function("p eval hypercube", |b| {
        b.iter(|| p.eval_hypercube(black_box(point)));
    });
    let point = HypercubePoint::from_usize(0b1011101_1_10110101);
    c.bench_function("p eval f then hypercube", |b| {
        b.iter(|| p.eval_field_then_hypercube(black_box(Fp2256::from(343)), black_box(point)));
    });
    c.bench_function("p new then fix f", |b| {
        b.iter(|| {
            P::<17, 9, 8, { 1 << 8 }, { 2 * (1 << 8) }, 16, Fp2256>::new(
                black_box(e),
                black_box(k),
                black_box(&HypercubeEvalPoly::new(&b_evals_0, 1 + 8)),
                black_box(&HypercubeEvalPoly::new(&b_evals_1, 1 + 8)),
                QConstructionMode::GrayCodes,
            )
            .fix_variable(black_box(Fp2256::from(343)))
        });
    });
    c.bench_function("p eval f then sum hypercube", |b| {
        b.iter(|| {
            P::<17, 9, 8, { 1 << 8 }, { 2 * (1 << 8) }, 16, Fp2256>::new(
                black_box(e),
                black_box(k),
                black_box(&HypercubeEvalPoly::new(&b_evals_0, 1 + 8)),
                black_box(&HypercubeEvalPoly::new(&b_evals_1, 1 + 8)),
                QConstructionMode::GrayCodes,
            )
            .eval_field_then_sum_hypercube(black_box(Fp2256::from(343)))
        });
    });

    c.bench_function("pad walk 705", |b| {
        b.iter(|| pad_walk::<p434::Fp2434>(black_box(1024-705)));
    });
    c.bench_function("pad walk 774", |b| {
        b.iter(|| pad_walk::<p503::Fp2503>(black_box(1024-774)));
    });
    c.bench_function("pad walk 1010", |b| {
        b.iter(|| pad_walk::<p610::Fp2610>(black_box(1024-1010)));
    });
    c.bench_function("pad walk 1280", |b| {
        b.iter(|| pad_walk::<p751::Fp2751>(black_box(2048-1280)));
    });

    c.bench_function("exclude walk 705", |b| {
        b.iter(|| measure_multiplications_for_exclusion(black_box(1024-705), black_box(p434::Fp2434::from(4123))));
    });
    c.bench_function("exclude walk 774", |b| {
        b.iter(|| measure_multiplications_for_exclusion(black_box(1024-774), black_box(p503::Fp2503::from(4123))));
    });
    c.bench_function("exclude walk 1010", |b| {
        b.iter(|| measure_multiplications_for_exclusion(black_box(1024-1010), black_box(p610::Fp2610::from(4123))));
    });
    c.bench_function("exclude walk 1280", |b| {
        b.iter(|| measure_multiplications_for_exclusion(black_box(2048-1280), black_box(p751::Fp2751::from(4123))));
    });

    // Mask
    let mask_evals: [_; 1 + 3 * 17] = array::from_fn(|i| Fp2256::from(i as u64));
    let mask = Mask::<17, _>::new(&mask_evals, &mask_evals);
    let point = [Fp2256::from(24234); 17];
    c.bench_function("mask eval", |b| {
        b.iter(|| mask.eval(black_box(&point)));
    });
    let point: HypercubePoint = HypercubePoint::from_usize(0b1011101_1_10110101);
    c.bench_function("mask eval hypercube", |b| {
        b.iter(|| mask.eval_hypercube(black_box(point)));
    });
    let f = Fp2256::from(2352);
    let point: HypercubePoint = HypercubePoint::from_usize(0b011101_1_10110101);
    c.bench_function("mask eval field then hypercube", |b| {
        b.iter(|| mask.eval_field_then_hypercube(black_box(f), black_box(point)));
    });

    let polynomial_interpolations: Vec<_> = b_evals_0.iter().chain(b_evals_1.iter()).copied().collect();
    c.bench_function("interpolate polynomial", |b| {
        b.iter(|| interpolate_polynomial(black_box(&polynomial_interpolations)));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
