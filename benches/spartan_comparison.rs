use ark_ff::Field;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use polynomial_proving::{
    CompressedPrivateKey, JInvariantPublicKey, PublicKey, RadicalPublicKey, RunForParamsConfig, expand_keys, field::Fp2253, interpolate_cosets, prove, verify_pok,
};
use util::{
    algebra::field::{p434, p503, p610, p751, sqisign, FftField},
    random_oracle::RandomOracle,
    CODE_RATE,
};

#[cfg(feature = "rayon")]
const PARALLEL: &str = "(parallel)";
#[cfg(not(feature = "rayon"))]
const PARALLEL: &str = "(single thread)";

trait FieldName {
    const FIELD_NAME: &str;
}

impl FieldName for sqisign::level_i::Fp2251 {
    const FIELD_NAME: &str = "5·2²⁴⁸-1 (SQIsign I)";
}
impl FieldName for sqisign::level_iii::Fp2383 {
    const FIELD_NAME: &str = "65·2³⁷⁶-1 (SQIsign III)";
}
impl FieldName for sqisign::level_v::Fp2505 {
    const FIELD_NAME: &str = "27·2⁵⁰⁰-1 (SQIsign V)";
}
impl FieldName for p434::Fp2434 {
    const FIELD_NAME: &str = "2²¹⁶·3¹³⁷-1 (p434)";
}
impl FieldName for p503::Fp2503 {
    const FIELD_NAME: &str = "2²⁵⁰·3¹⁵⁹-1 (p503)";
}
impl FieldName for p610::Fp2610 {
    const FIELD_NAME: &str = "2³⁰⁵·3¹⁹²-1 (p610)";
}
impl FieldName for p751::Fp2751 {
    const FIELD_NAME: &str = "2³⁷²·3²³⁹-1 (p751)";
}
impl FieldName for Fp2253 {
    const FIELD_NAME: &str = "curve25519 Field";
}

trait PublicKeyName {
    const NAME: &str;
}
impl<F: Field> PublicKeyName for RadicalPublicKey<F> {
    const NAME: &str = "pk=(A, C)";
}
impl<F: Field> PublicKeyName for JInvariantPublicKey<F> {
    const NAME: &str = "pk=j";
}

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
    F: FftField + FieldName,
    PK: PublicKey<F> + PublicKeyName,
>(
    c: &mut Criterion,
) {
    // Generate the public and private key
    let compressed_private_key =
        CompressedPrivateKey::<PATH_LENGTH_DIV_64>::rand(StdRng::from_entropy());
    let (public_key, private_key) = expand_keys(&compressed_private_key);
    assert!(private_key.check(&public_key));

    let id_for_function = |func| format!(
        "{}, e={}, λ={}, {}, {}: {func}",
        F::FIELD_NAME,
        PATH_LENGTH,
        SECURITY_BITS,
        PK::NAME,
        PARALLEL
    );

    let keygen_id = id_for_function("keygen");
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

    let prove_id = id_for_function("prove");
    let verify_id = id_for_function("verify");

    let cosets = interpolate_cosets(LOG_2_PATH_LENGTH + 2);

    c.bench_function(&prove_id, |b| {
        b.iter(|| {
            prove::<
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
                RadicalPublicKey<F>,
            >(
                black_box(&public_key),
                black_box(&private_key),
                black_box(&random_oracle),
                polynomial_proving::MaskCheckMode::Additional,
                &cosets
            )
        })
    });

    let proof = prove::<
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
        RadicalPublicKey<F>,
    >(
        black_box(&public_key),
        black_box(&private_key),
        black_box(&random_oracle),
        polynomial_proving::MaskCheckMode::Additional,
        &cosets
    )
    .unwrap();

    
    let param_gen_id = id_for_function("public parameter generation");
    c.bench_function(&param_gen_id, |b| {
        b.iter(|| interpolate_cosets::<F>(LOG_2_PATH_LENGTH + 2))
    });

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
                RadicalPublicKey<F>,
            >(
                black_box(&public_key),
                black_box(proof.clone()),
                black_box(&random_oracle),
                polynomial_proving::MaskCheckMode::Additional,
                &cosets
            )
        })
    });
}

pub fn criterion_benchmark(c: &mut Criterion) {
    {
        const CFG_SQISIGN_8: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 8,
            security_bits: 128,
            commitment_size: 96,
        };
        test_prove_verify::<
            { CFG_SQISIGN_8.variable_count() },
            { CFG_SQISIGN_8.path_length() },
            { CFG_SQISIGN_8.path_length_div_64() },
            { CFG_SQISIGN_8.path_length_times_two() },
            { CFG_SQISIGN_8.path_length_times_four() },
            { CFG_SQISIGN_8.log_path_length() },
            { CFG_SQISIGN_8.log_path_length_plus_one() },
            { CFG_SQISIGN_8.security_bits() },
            { CFG_SQISIGN_8.commitment_size },
            { CFG_SQISIGN_8.q_variable_count() },
            { CFG_SQISIGN_8.final_round_evaluations() },
            sqisign::level_i::Fp2251,
            JInvariantPublicKey<_>,
        >(c);
    }
    {
        const CFG_SQISIGN_9: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 9,
            security_bits: 128,
            commitment_size: 96,
        };
        test_prove_verify::<
            { CFG_SQISIGN_9.variable_count() },
            { CFG_SQISIGN_9.path_length() },
            { CFG_SQISIGN_9.path_length_div_64() },
            { CFG_SQISIGN_9.path_length_times_two() },
            { CFG_SQISIGN_9.path_length_times_four() },
            { CFG_SQISIGN_9.log_path_length() },
            { CFG_SQISIGN_9.log_path_length_plus_one() },
            { CFG_SQISIGN_9.security_bits() },
            { CFG_SQISIGN_9.commitment_size },
            { CFG_SQISIGN_9.q_variable_count() },
            { CFG_SQISIGN_9.final_round_evaluations() },
            sqisign::level_i::Fp2251,
            JInvariantPublicKey<_>,
        >(c);
    }
    {
        const CFG_SQISIGN_10: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 10,
            security_bits: 128,
            commitment_size: 96,
        };
        test_prove_verify::<
            { CFG_SQISIGN_10.variable_count() },
            { CFG_SQISIGN_10.path_length() },
            { CFG_SQISIGN_10.path_length_div_64() },
            { CFG_SQISIGN_10.path_length_times_two() },
            { CFG_SQISIGN_10.path_length_times_four() },
            { CFG_SQISIGN_10.log_path_length() },
            { CFG_SQISIGN_10.log_path_length_plus_one() },
            { CFG_SQISIGN_10.security_bits() },
            { CFG_SQISIGN_10.commitment_size },
            { CFG_SQISIGN_10.q_variable_count() },
            { CFG_SQISIGN_10.final_round_evaluations() },
            sqisign::level_i::Fp2251,
            JInvariantPublicKey<_>,
        >(c);
    }
    {
        const CFG_SQISIGN_11: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 11,
            security_bits: 128,
            commitment_size: 96,
        };
        test_prove_verify::<
            { CFG_SQISIGN_11.variable_count() },
            { CFG_SQISIGN_11.path_length() },
            { CFG_SQISIGN_11.path_length_div_64() },
            { CFG_SQISIGN_11.path_length_times_two() },
            { CFG_SQISIGN_11.path_length_times_four() },
            { CFG_SQISIGN_11.log_path_length() },
            { CFG_SQISIGN_11.log_path_length_plus_one() },
            { CFG_SQISIGN_11.security_bits() },
            { CFG_SQISIGN_11.commitment_size },
            { CFG_SQISIGN_11.q_variable_count() },
            { CFG_SQISIGN_11.final_round_evaluations() },
            sqisign::level_i::Fp2251,
            JInvariantPublicKey<_>,
        >(c);
    }
    {
        const CFG_SQISIGN_12: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 12,
            security_bits: 128,
            commitment_size: 96,
        };
        test_prove_verify::<
            { CFG_SQISIGN_12.variable_count() },
            { CFG_SQISIGN_12.path_length() },
            { CFG_SQISIGN_12.path_length_div_64() },
            { CFG_SQISIGN_12.path_length_times_two() },
            { CFG_SQISIGN_12.path_length_times_four() },
            { CFG_SQISIGN_12.log_path_length() },
            { CFG_SQISIGN_12.log_path_length_plus_one() },
            { CFG_SQISIGN_12.security_bits() },
            { CFG_SQISIGN_12.commitment_size },
            { CFG_SQISIGN_12.q_variable_count() },
            { CFG_SQISIGN_12.final_round_evaluations() },
            sqisign::level_i::Fp2251,
            JInvariantPublicKey<_>,
        >(c);
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
