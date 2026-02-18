use std::error::Error;

use polynomial_proving::{run_for_params, RunForParamsConfig};
use util::algebra::field::{arkfield::Fp2256, p434};

#[test]
fn prove_and_verify() -> Result<(), Box<dyn Error>> {
    run_for_params::<
        17,
        { 1 << 8 },
        { (1 << 8) / 64 },
        { 1 << 9 },
        { 1 << 10 },
        8,
        9,
        128,
        96,
        16,
        { 4 + 6 * 8 + 3 },
        Fp2256,
    >(polynomial_proving::MaskCheckMode::Additional)?;
    Ok(())
}

#[test]
fn fp434() -> Result<(), Box<dyn Error>> {
    const CFG: RunForParamsConfig = RunForParamsConfig {
        log_path_length: 10,
        security_bits: 128,
        commitment_size: 142,
    };
    run_for_params::<
        { CFG.variable_count() },
        { CFG.path_length() },
        { CFG.path_length_div_64() },
        { CFG.path_length_times_two() },
        { CFG.path_length_times_four() },
        { CFG.log_path_length() },
        { CFG.log_path_length_plus_one() },
        { CFG.security_bits() },
        { CFG.commitment_size },
        { CFG.q_variable_count() },
        { CFG.final_round_evaluations() },
        p434::Fp2434,
    >(polynomial_proving::MaskCheckMode::Additional)?;
    Ok(())
}
