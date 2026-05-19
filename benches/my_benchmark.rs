use ark_ff::Field;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use polynomial_proving::{
    expand_keys, prove, verify_pok, CompressedPrivateKey, JInvariantPublicKey, PublicKey,
    RadicalPublicKey, RunForParamsConfig,
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
                RadicalPublicKey<F>,
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
            RadicalPublicKey<_>,
        >(c);
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
            JInvariantPublicKey<_>,
        >(c);
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
            RadicalPublicKey<_>,
        >(c);
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
            JInvariantPublicKey<_>,
        >(c);
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
            RadicalPublicKey<_>,
        >(c);
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
            JInvariantPublicKey<_>,
        >(c);
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
            RadicalPublicKey<_>,
        >(c);
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
            JInvariantPublicKey<_>,
        >(c);
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
            RadicalPublicKey<_>,
        >(c);
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
            JInvariantPublicKey<_>,
        >(c);
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
            RadicalPublicKey<_>,
        >(c);
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
            JInvariantPublicKey<_>,
        >(c);
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
            RadicalPublicKey<_>,
        >(c);
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
            JInvariantPublicKey<_>,
        >(c);
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
