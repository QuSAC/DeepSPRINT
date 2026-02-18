use ark_std::rand::{rngs::StdRng, SeedableRng};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use polynomial_proving::{
    do_sumcheck_pok, expand_keys, verify_pok, CompressedPrivateKey, RunForParamsConfig,
};
use util::{
    algebra::field::{p434, p503, p610, p751, sqisign, FftField},
    random_oracle::RandomOracle,
    CODE_RATE,
};

fn test_prove_verify<
    const VARIABLE_COUNT: usize,
    const PATH_LENGTH: usize,
    const PATH_LENGTH_DIV_64: usize,
    const PATH_LENGTH_TIMES_TWO: usize,
    const PATH_LENGTH_TIMES_FOUR: usize,
    const LOG_2_PATH_LENGTH: usize,
    const LOG_2_PATH_LENGTH_PLUS_ONE: usize,
    const SECURITY_BITS: usize,
    const COMMITMENT_SIZE: usize,
    const Q_VARIABLE_COUNT: usize,
    const FINAL_ROUND_EVALUATIONS: usize,
    F: FftField,
>(
    c: &mut Criterion,
    name: &str,
) {
    // Generate the public and private key
    let compressed_private_key =
        CompressedPrivateKey::<PATH_LENGTH_DIV_64>::rand(StdRng::from_entropy());
    let (public_key, private_key) = expand_keys(&compressed_private_key);
    assert!(private_key.check(&public_key));

    let keygen_id = format!("{name}: keygen ({} bits security)", SECURITY_BITS);
    c.bench_function(&keygen_id, |b| {
        b.iter(|| {
            let compressed_private_key =
                CompressedPrivateKey::<PATH_LENGTH_DIV_64>::rand(StdRng::from_entropy());
            expand_keys::<PATH_LENGTH_DIV_64, F>(&compressed_private_key)
        })
    });

    let random_oracle = RandomOracle::new(
        &mut StdRng::from_seed([0; 32]),
        LOG_2_PATH_LENGTH + 2,
        SECURITY_BITS / CODE_RATE,
    );

    let proof_id = format!("{name}: prove ({} bits security)", SECURITY_BITS);
    let verify_id = format!("{name}: verify ({} bits security)", SECURITY_BITS);

    c.bench_function(&proof_id, |b| {
        b.iter(|| {
            do_sumcheck_pok::<
                VARIABLE_COUNT,
                PATH_LENGTH,
                PATH_LENGTH_TIMES_TWO,
                PATH_LENGTH_TIMES_FOUR,
                LOG_2_PATH_LENGTH,
                LOG_2_PATH_LENGTH_PLUS_ONE,
                Q_VARIABLE_COUNT,
                COMMITMENT_SIZE,
                FINAL_ROUND_EVALUATIONS,
                F,
            >(
                black_box(&private_key),
                black_box(&random_oracle),
                polynomial_proving::MaskCheckMode::Additional,
            )
        })
    });

    let proof = do_sumcheck_pok::<
        VARIABLE_COUNT,
        PATH_LENGTH,
        PATH_LENGTH_TIMES_TWO,
        PATH_LENGTH_TIMES_FOUR,
        LOG_2_PATH_LENGTH,
        LOG_2_PATH_LENGTH_PLUS_ONE,
        Q_VARIABLE_COUNT,
        COMMITMENT_SIZE,
        FINAL_ROUND_EVALUATIONS,
        F,
    >(
        black_box(&private_key),
        black_box(&random_oracle),
        polynomial_proving::MaskCheckMode::Additional,
    )
    .unwrap();

    c.bench_function(&verify_id, |b| {
        b.iter(|| {
            verify_pok::<
                VARIABLE_COUNT,
                PATH_LENGTH,
                PATH_LENGTH_TIMES_TWO,
                PATH_LENGTH_TIMES_FOUR,
                LOG_2_PATH_LENGTH,
                LOG_2_PATH_LENGTH_PLUS_ONE,
                Q_VARIABLE_COUNT,
                COMMITMENT_SIZE,
                FINAL_ROUND_EVALUATIONS,
                F,
            >(
                black_box(&public_key),
                black_box(proof.clone()),
                black_box(&random_oracle),
                polynomial_proving::MaskCheckMode::Additional,
            )
        })
    });
}

pub fn criterion_benchmark(c: &mut Criterion) {
    {
        const CFG_SQISIGN_I: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 8,
            security_bits: 128,
            commitment_size: 96,
        };
        test_prove_verify::<
            { CFG_SQISIGN_I.variable_count() },
            { CFG_SQISIGN_I.path_length() },
            { CFG_SQISIGN_I.path_length_div_64() },
            { CFG_SQISIGN_I.path_length_times_two() },
            { CFG_SQISIGN_I.path_length_times_four() },
            { CFG_SQISIGN_I.log_path_length() },
            { CFG_SQISIGN_I.log_path_length_plus_one() },
            { CFG_SQISIGN_I.security_bits() },
            { CFG_SQISIGN_I.commitment_size },
            { CFG_SQISIGN_I.q_variable_count() },
            { CFG_SQISIGN_I.final_round_evaluations() },
            sqisign::level_i::Fp2251,
        >(c, "Fp2251 (SQISign I)");
    }
    {
        const CFG_SQISIGN_III: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 9,
            security_bits: 192,
            commitment_size: 128,
        };
        test_prove_verify::<
            { CFG_SQISIGN_III.variable_count() },
            { CFG_SQISIGN_III.path_length() },
            { CFG_SQISIGN_III.path_length_div_64() },
            { CFG_SQISIGN_III.path_length_times_two() },
            { CFG_SQISIGN_III.path_length_times_four() },
            { CFG_SQISIGN_III.log_path_length() },
            { CFG_SQISIGN_III.log_path_length_plus_one() },
            { CFG_SQISIGN_III.security_bits() },
            { CFG_SQISIGN_III.commitment_size },
            { CFG_SQISIGN_III.q_variable_count() },
            { CFG_SQISIGN_III.final_round_evaluations() },
            sqisign::level_iii::Fp2383,
        >(c, "Fp2383 (SQISign III)");
    }
    {
        const CFG_SQISIGN_V: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 9,
            security_bits: 256,
            commitment_size: 160,
        };
        test_prove_verify::<
            { CFG_SQISIGN_V.variable_count() },
            { CFG_SQISIGN_V.path_length() },
            { CFG_SQISIGN_V.path_length_div_64() },
            { CFG_SQISIGN_V.path_length_times_two() },
            { CFG_SQISIGN_V.path_length_times_four() },
            { CFG_SQISIGN_V.log_path_length() },
            { CFG_SQISIGN_V.log_path_length_plus_one() },
            { CFG_SQISIGN_V.security_bits() },
            { CFG_SQISIGN_V.commitment_size },
            { CFG_SQISIGN_V.q_variable_count() },
            { CFG_SQISIGN_V.final_round_evaluations() },
            sqisign::level_v::Fp2505,
        >(c, "Fp2505 (SQISign V)");
    }
    {
        const CFG_P434: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 10,
            security_bits: 128,
            commitment_size: 142,
        };
        test_prove_verify::<
            { CFG_P434.variable_count() },
            { CFG_P434.path_length() },
            { CFG_P434.path_length_div_64() },
            { CFG_P434.path_length_times_two() },
            { CFG_P434.path_length_times_four() },
            { CFG_P434.log_path_length() },
            { CFG_P434.log_path_length_plus_one() },
            { CFG_P434.security_bits() },
            { CFG_P434.commitment_size },
            { CFG_P434.q_variable_count() },
            { CFG_P434.final_round_evaluations() },
            p434::Fp2434,
        >(c, "Fp2434 (path length 1024, λ=128)");
    }
    {
        const CFG_P503: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 10,
            security_bits: 128,
            commitment_size: 158,
        };
        test_prove_verify::<
            { CFG_P503.variable_count() },
            { CFG_P503.path_length() },
            { CFG_P503.path_length_div_64() },
            { CFG_P503.path_length_times_two() },
            { CFG_P503.path_length_times_four() },
            { CFG_P503.log_path_length() },
            { CFG_P503.log_path_length_plus_one() },
            { CFG_P503.security_bits() },
            { CFG_P503.commitment_size },
            { CFG_P503.q_variable_count() },
            { CFG_P503.final_round_evaluations() },
            p503::Fp2503,
        >(c, "Fp2503 (path length 1024, λ=128)");
    }
    {
        const CFG_P610: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 10,
            security_bits: 192,
            commitment_size: 186,
        };
        test_prove_verify::<
            { CFG_P610.variable_count() },
            { CFG_P610.path_length() },
            { CFG_P610.path_length_div_64() },
            { CFG_P610.path_length_times_two() },
            { CFG_P610.path_length_times_four() },
            { CFG_P610.log_path_length() },
            { CFG_P610.log_path_length_plus_one() },
            { CFG_P610.security_bits() },
            { CFG_P610.commitment_size },
            { CFG_P610.q_variable_count() },
            { CFG_P610.final_round_evaluations() },
            p610::Fp2610,
        >(c, "Fp2610 (path length 1024, λ=192)");
    }
    {
        const CFG_P751: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 11,
            security_bits: 256,
            commitment_size: 220,
        };
        test_prove_verify::<
            { CFG_P751.variable_count() },
            { CFG_P751.path_length() },
            { CFG_P751.path_length_div_64() },
            { CFG_P751.path_length_times_two() },
            { CFG_P751.path_length_times_four() },
            { CFG_P751.log_path_length() },
            { CFG_P751.log_path_length_plus_one() },
            { CFG_P751.security_bits() },
            { CFG_P751.commitment_size },
            { CFG_P751.q_variable_count() },
            { CFG_P751.final_round_evaluations() },
            p751::Fp2751,
        >(c, "Fp2751 (path length 2048, λ=256)");
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
