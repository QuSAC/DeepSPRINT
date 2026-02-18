use std::error::Error;

use polynomial_proving::{run_for_params, MaskCheckMode, RunForParamsConfig};
use util::algebra::field::{p434, p503, p610, p751, sqisign};

fn main() -> Result<(), Box<dyn Error>> {
    let mode = MaskCheckMode::Additional;

    println!("SQISign I:");
    const CFG_SQISIGN_I: RunForParamsConfig = RunForParamsConfig {
        log_path_length: 8,
        security_bits: 128,
        commitment_size: 96,
    };
    run_for_params::<
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
    >(mode)?;

    println!("SQISign III:");
    const CFG_SQISIGN_III: RunForParamsConfig = RunForParamsConfig {
        log_path_length: 9,
        security_bits: 192,
        commitment_size: 128,
    };
    run_for_params::<
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
    >(mode)?;

    println!("SQISign V:");
    const CFG_SQISIGN_V: RunForParamsConfig = RunForParamsConfig {
        log_path_length: 9,
        security_bits: 256,
        commitment_size: 160,
    };
    run_for_params::<
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
    >(mode)?;

    println!("p434:");
    const CFG_P434: RunForParamsConfig = RunForParamsConfig {
        log_path_length: 10,
        security_bits: 128,
        commitment_size: 142,
    };
    run_for_params::<
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
    >(mode)?;

    {
        println!("503:");
        const CFG_P503: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 10,
            security_bits: 128,
            commitment_size: 158,
        };
        run_for_params::<
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
        >(mode)?;
    }
    {
        println!("p610:");
        const CFG_P610: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 10,
            security_bits: 192,
            commitment_size: 186,
        };
        run_for_params::<
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
        >(mode)?;
    }
    {
        println!("p751:");
        const CFG_P751: RunForParamsConfig = RunForParamsConfig {
            log_path_length: 11,
            security_bits: 256,
            commitment_size: 220,
        };
        run_for_params::<
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
        >(mode)?;
    }
    Ok(())
}
