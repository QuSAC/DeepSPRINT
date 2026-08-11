use std::error::Error;

use polynomial_proving::{
    run_for_params, JInvariantPublicKeyConfig, MaskCheckMode, RadicalPublicKeyConfig,
    RunForParamsConfig, RunForParamsPublicKey,
};
use util::algebra::field::{p434, p503, p610, p751, sqisign};

trait PublicKeyDescriptor {
    fn print();
}

impl PublicKeyDescriptor for RadicalPublicKeyConfig {
    fn print() {
        println!("Radical public key");
    }
}
impl PublicKeyDescriptor for JInvariantPublicKeyConfig {
    fn print() {
        println!("j-invariant public key");
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    fn run_all_for_public_key<PK: RunForParamsPublicKey + PublicKeyDescriptor>(
    ) -> Result<(), Box<dyn Error>> {
        for mode in [MaskCheckMode::Additional, MaskCheckMode::InsidePCS].into_iter() {
            PK::print();
            println!(
                "Mode: {}",
                match mode {
                    MaskCheckMode::Additional => "Masked mask",
                    MaskCheckMode::InsidePCS => "Optimized mask opening",
                }
            );
            const CFG_SQISIGN_8: RunForParamsConfig = RunForParamsConfig {
                log_path_length: 8,
                security_bits: 128,
                commitment_size: 96,
            };
            run_for_params::<
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
                PK,
            >(mode)?;

            const CFG_SQISIGN_9: RunForParamsConfig = RunForParamsConfig {
                log_path_length: 9,
                security_bits: 128,
                commitment_size: 96,
            };
            run_for_params::<
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
                PK,
            >(mode)?;

            const CFG_SQISIGN_10: RunForParamsConfig = RunForParamsConfig {
                log_path_length: 10,
                security_bits: 128,
                commitment_size: 96,
            };
            run_for_params::<
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
                PK,
            >(mode)?;

            const CFG_SQISIGN_11: RunForParamsConfig = RunForParamsConfig {
                log_path_length: 11,
                security_bits: 128,
                commitment_size: 96,
            };
            run_for_params::<
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
                PK,
            >(mode)?;

            {
                const CFG_SQISIGN_12: RunForParamsConfig = RunForParamsConfig {
                    log_path_length: 12,
                    security_bits: 128,
                    commitment_size: 96,
                };
                run_for_params::<
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
                    PK,
                >(mode)?;
            }
        }
        Ok(())
    }

    run_all_for_public_key::<JInvariantPublicKeyConfig>()?;
    run_all_for_public_key::<RadicalPublicKeyConfig>()?;
    Ok(())
}
