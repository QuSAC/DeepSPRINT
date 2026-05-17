mod keys;
pub mod polynomials;

pub use keys::{
    expand_keys, pad_walk, CompressedPrivateKey, JInvariantPublicKey, PrivateKey, PublicKey,
    RadicalPublicKey,
};

use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;
use ark_std::rand::SeedableRng;
use spongefish::codecs::arkworks_algebra::{
    FieldDomainSeparator, FieldToUnitDeserialize, FieldToUnitSerialize, UnitToField,
};
use spongefish::{
    ByteDomainSeparator, BytesToUnitDeserialize, BytesToUnitSerialize, DefaultHash, DomainSeparator,
};
use std::error::Error;
use std::{array, default};
use util::algebra::field::FftField;
use util::CODE_RATE;

use deepfold::{prover::Prover, Commit, SerializableCommit};
use util::{
    algebra::{coset::Coset, polynomial::MultilinearPolynomial},
    random_oracle::RandomOracle,
};

use crate::polynomials::piop_polynomials::{Checker, Mask, QConstructionMode, P, Q};
use crate::polynomials::{HypercubeEvalPoly, PiopPolynomial};

/// Choose how the mask is checked. Either the linear combination is checked
/// inside the polynomial commitment scheme, or we do this seperately at the
/// cost of proof size.
#[derive(PartialEq, Eq, Copy, Clone, Hash)]
pub enum MaskCheckMode {
    InsidePCS,
    Additional,
}

#[derive(default::Default)]
pub struct RunForParamsConfig {
    pub log_path_length: usize,
    pub security_bits: usize,
    /// The easiest way to determine is just to set it to a random value and
    /// check the error.
    pub commitment_size: usize,
}

impl RunForParamsConfig {
    pub const fn variable_count(&self) -> usize {
        2 * self.log_path_length + 1
    }
    pub const fn path_length(&self) -> usize {
        1 << self.log_path_length
    }
    pub const fn path_length_div_64(&self) -> usize {
        self.path_length() / 64
    }
    pub const fn path_length_times_two(&self) -> usize {
        self.path_length() * 2
    }
    pub const fn path_length_times_four(&self) -> usize {
        self.path_length() * 4
    }
    pub const fn log_path_length(&self) -> usize {
        self.log_path_length
    }
    pub const fn log_path_length_plus_one(&self) -> usize {
        self.log_path_length + 1
    }
    pub const fn security_bits(&self) -> usize {
        self.security_bits
    }
    pub const fn q_variable_count(&self) -> usize {
        self.log_path_length() * 2
    }
    pub const fn final_round_evaluations(&self) -> usize {
        4 + 6 * self.log_path_length + 3
    }
}

pub trait RunForParamsPublicKey {
    type PK<F: Field>: PublicKey<F> + From<RadicalPublicKey<F>>;
}
pub struct RadicalPublicKeyConfig;
impl RunForParamsPublicKey for RadicalPublicKeyConfig {
    type PK<F: Field> = RadicalPublicKey<F>;
}
pub struct JInvariantPublicKeyConfig;
impl RunForParamsPublicKey for JInvariantPublicKeyConfig {
    type PK<F: Field> = JInvariantPublicKey<F>;
}

/// A convenience function to prove and then verify.
pub fn run_for_params<
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
    PK: RunForParamsPublicKey,
>(
    mask_check_mode: MaskCheckMode,
) -> Result<(), Box<dyn Error>> {
    // Generate the public and private key
    let compressed_private_key =
        CompressedPrivateKey::<PATH_LENGTH_DIV_64>::rand(StdRng::from_entropy());
    let (public_key, private_key) = expand_keys(&compressed_private_key);
    assert!(private_key.check(&public_key));
    let public_key: PK::PK<F> = public_key.into();

    let random_oracle = RandomOracle::new(
        &mut StdRng::from_entropy(),
        LOG_2_PATH_LENGTH + 2,
        SECURITY_BITS / CODE_RATE,
    );

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
        PK::PK<F>,
    >(&public_key, &private_key, &random_oracle, mask_check_mode)?;

    println!(
        "Proof size: {} (NARG: {}, PCS: {})",
        proof.size(),
        proof.narg.compressed_size(),
        proof.final_evaluation_proof.compressed_size()
    );

    let verified = verify_pok::<
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
        PK::PK<F>,
    >(&public_key, proof, &random_oracle, mask_check_mode)?;

    assert!(verified);

    Ok(())
}

/// The IO pattern determines the pattern of communication for the purposes of
/// Fiat-Shamir.
fn pok_io_pattern<
    const LOG_2_PATH_LENGTH: usize,
    const COMMITMENT_SIZE: usize,
    const FINAL_ROUND_EVALUATIONS: usize,
    F: Field,
    PK: PublicKey<F>,
>(
    mask_check_mode: MaskCheckMode,
) -> DomainSeparator<DefaultHash> {
    let mut ds = DomainSeparator::<DefaultHash>::new("Proof of Knowledge");
    ds = ds.add_bytes(MESSAGE.len(), "message");
    ds = ds.add_bytes(COMMITMENT_SIZE, "commitment");
    ds = FieldDomainSeparator::<F>::challenge_scalars(ds, 2, "e");
    ds = FieldDomainSeparator::<F>::challenge_scalars(ds, LOG_2_PATH_LENGTH, "k");
    if PK::CHECKING_RANDOMNESS_SIZE > 0 {
        for _ in 0..2 {
            ds = FieldDomainSeparator::<F>::challenge_scalars(
                ds,
                PK::CHECKING_RANDOMNESS_SIZE,
                "j invariant checking",
            );
        }
    }

    // The rounds start with a verifier challenge and then the prover message.
    // The first round has no challenge. This is why we do one additional round
    // here.
    for _ in 0..(2 * LOG_2_PATH_LENGTH + 1) {
        ds = FieldDomainSeparator::<F>::add_scalars(ds, 3, "sumcheck_round");
        ds = FieldDomainSeparator::<F>::challenge_scalars(ds, 1, "challenge");
    }
    let final_evaluation_count = match mask_check_mode {
        MaskCheckMode::Additional => FINAL_ROUND_EVALUATIONS,
        MaskCheckMode::InsidePCS => 4 + 1,
    };
    match mask_check_mode {
        MaskCheckMode::Additional => {
            ds = FieldDomainSeparator::<F>::add_scalars(ds, 2, "mask_points");
            ds = FieldDomainSeparator::<F>::challenge_scalars(ds, 1, "mask_challenge");
            ds = FieldDomainSeparator::<F>::add_scalars(
                ds,
                final_evaluation_count,
                "final evaluations",
            );
        }
        MaskCheckMode::InsidePCS => {
            ds = FieldDomainSeparator::<F>::add_scalars(
                ds,
                final_evaluation_count,
                "final evaluations",
            );
        }
    }

    if PK::Checker::OPENING_LEN > 0 {
        ds = FieldDomainSeparator::<F>::add_scalars(
            ds,
            PK::Checker::OPENING_LEN,
            "start statement checking",
        );
        ds = FieldDomainSeparator::<F>::add_scalars(
            ds,
            PK::Checker::OPENING_LEN,
            "end statement checking",
        );
    }

    // The opening proof has a variable size, so we will not include it in the
    // transcript. This is not an issue, since it is the last element to be
    // sent: no verifier challenges are derived from it.
    // ds = ds.absorb(63712, "proof");
    ds
}

pub fn measure_multiplications_for_exclusion<F: Field>(m: usize, mut f: F) -> F {
    for _ in 0..(2 * (m - 1)) {
        f = f * f;
    }
    for _ in 0..(3 * (m - 1)) {
        f = f + f;
    }
    f
}

#[derive(Clone)]
pub struct Proof<F: FftField> {
    pub narg: Vec<u8>,
    pub final_evaluation_proof: deepfold::Proof<F>,
}

impl<F: FftField> Proof<F> {
    #[inline]
    pub fn size(&self) -> usize {
        self.narg.compressed_size() + self.final_evaluation_proof.compressed_size()
    }
}

const MESSAGE: [u8; 32] = [0; 32];

/// Prove knowledge of an isogeny chain. `a` and `b` are vectors with elliptic
/// curve coefficients. Their length must be equal and a power of two.
///
/// # Panic
/// Panics if the length of `a` and `b` is not the same power of 2.
pub fn prove<
    const VARIABLE_COUNT: usize,
    const PATH_LENGTH: usize,
    const PATH_LENGTH_TIMES_TWO: usize,
    const PATH_LENGTH_TIMES_FOUR: usize,
    const LOG_2_PATH_LENGTH: usize,
    const LOG_2_PATH_LENGTH_PLUS_ONE: usize,
    const Q_VARIABLE_COUNT: usize,
    const COMMITMENT_SIZE: usize,
    const FINAL_ROUND_EVALUATIONS: usize,
    F: FftField,
    PK: PublicKey<F>,
>(
    public_key: &PK,
    private_key: &PrivateKey<F>,
    random_oracle: &RandomOracle<F>,
    mask_check_mode: MaskCheckMode,
) -> Result<Proof<F>, Box<dyn Error>> {
    // Sanity checking and parameter generation
    debug_assert_eq!(private_key.len(), PATH_LENGTH);

    // This is the object that enables Fiat-Shamir for us. We feed it data and
    // it spits out challenges.
    let mut merlin =
        pok_io_pattern::<LOG_2_PATH_LENGTH, COMMITMENT_SIZE, FINAL_ROUND_EVALUATIONS, F, PK>(
            mask_check_mode,
        )
        .to_prover_state();

    // Add message
    merlin.add_bytes(&MESSAGE)?;

    // The polynomial commitment has `num_variables` variables, plus
    // - The first two variables select the polynomial to query:
    //   - `00` is A
    //   - `01` is C
    //   - `10` & `11` are the mask
    let mut b_0: Vec<_> = (0..PATH_LENGTH)
        .map(|_| F::rand(merlin.rng()))
        // Then comes the real witness
        .chain(private_key.path_a().iter().copied())
        .collect();
    debug_assert_eq!(b_0.len(), PATH_LENGTH_TIMES_TWO);

    let mut b_1: Vec<_> = (0..PATH_LENGTH)
        .map(|_| F::rand(merlin.rng()))
        // Then comes the real witness
        .chain(private_key.path_c().iter().copied())
        .collect();

    public_key.fix_checker_relation(LOG_2_PATH_LENGTH, &mut b_0, &mut b_1);

    // Fix the constant term of the mask
    Mask::<VARIABLE_COUNT, F>::fix_constant_term(&mut b_0, &mut b_1);

    let b_0 = HypercubeEvalPoly::new(&b_0, 1 + LOG_2_PATH_LENGTH);
    let b_1 = HypercubeEvalPoly::new(&b_1, 1 + LOG_2_PATH_LENGTH);

    // Commit and compress
    let polynomial_evaluations = b_0
        .evals()
        .iter()
        .chain(b_1.evals().iter())
        .copied()
        .collect::<Vec<_>>();
    let polynomial = interpolate_polynomial(&polynomial_evaluations);
    let (prover, commitment) =
        deepfold_commit::<LOG_2_PATH_LENGTH, F>(polynomial, random_oracle, polynomial_evaluations);
    let mut commitment_out = Vec::new();
    let commitment = SerializableCommit::from(commitment);
    commitment.serialize_compressed(&mut commitment_out)?;
    debug_assert_eq!(commitment_out.len(), COMMITMENT_SIZE);
    merlin.add_bytes(&commitment_out)?;

    // Fetch the evaluation for e. We don't need to run the sumcheck protocol
    // over it, since the polynomial will be linear in e and should be zero
    // everywhere.
    let e = merlin.challenge_scalars()?;
    let k = merlin.challenge_scalars()?;
    let statement_checker_randomness = PK::get_checking_randomness(&mut merlin)?;

    let mut p = P::<
        '_,
        VARIABLE_COUNT,
        LOG_2_PATH_LENGTH_PLUS_ONE,
        LOG_2_PATH_LENGTH,
        PATH_LENGTH,
        PATH_LENGTH_TIMES_TWO,
        Q_VARIABLE_COUNT,
        F,
        PK,
    >::new(
        e,
        k,
        &b_0,
        &b_1,
        &statement_checker_randomness,
        QConstructionMode::GrayCodes,
        public_key,
    );

    // Test that the constraints actually hold
    debug_assert!(p.is_system_satisfied());

    // Do sumcheck
    let mut evaluation_points = [F::ZERO; VARIABLE_COUNT];
    //let mut sum = F::ZERO;
    for evaluation_point in evaluation_points.iter_mut() {
        //debug_assert_eq!(p.eval_field_then_sum_hypercube(F::ZERO) + p.eval_field_then_sum_hypercube(F::ONE), sum);
        let evaluations: [F; 3] =
            array::from_fn(|i| p.eval_field_then_sum_hypercube(F::from(1 + i as u64)));
        merlin.add_scalars(&evaluations)?;
        *evaluation_point = merlin.challenge_scalars::<1>()?[0];
        p.fix_variable(*evaluation_point);
        //sum = evaluate_poly_from_evals([sum - evaluations[0], evaluations[0], evaluations[1], evaluations[2]], *evaluation_point);
    }
    // All variables should now be exhausted.
    debug_assert_eq!(p.variable_count(), 0);

    // The final proof evaluations. We only send 3 in the optimized case.
    let evaluations = p.final_evaluations();
    if mask_check_mode == MaskCheckMode::InsidePCS {
        merlin.add_scalars(&evaluations)?;
    } else {
        // In the non-optimized case, we send all the coefficients of the mask
        let mask_at_y_0 = evaluations[4];
        let mask_at_y_1 = p.mask().eval_with_y_eq_1(&evaluation_points);
        merlin.add_scalars(&[mask_at_y_0, mask_at_y_1])?;
        let [r_y]: [F; 1] = merlin.challenge_scalars()?;
        let mut masked_coefficients = p.mask().masked_coefficients(r_y);
        // Note that we don't send the constant term M'(_, 0): this can be
        // derived from the other points. It will need to be checked.

        // We do send the four evaluations needed for checking the final
        // relation.
        masked_coefficients.extend_from_slice(&evaluations[..4]);
        merlin.add_scalars(&masked_coefficients)?;
    }

    if PK::Checker::OPENING_LEN > 0 {
        for checker in p.checkers() {
            merlin.add_scalars(checker.get_openings().as_ref())?;
        }
    }

    let narg = merlin.narg_string().to_vec();

    let mut eval_point = vec![evaluation_points[0], F::ZERO];
    eval_point.extend_from_slice(&evaluation_points[1..1 + LOG_2_PATH_LENGTH]);
    let final_evaluation_proof = prover.generate_proof(eval_point);

    Ok(Proof {
        narg,
        final_evaluation_proof,
    })
}

fn interpolate_cosets<F: FftField>(variable_num: usize) -> Vec<Coset<F>> {
    use util::CODE_RATE;
    let mut interpolate_cosets = vec![Coset::new(1 << (variable_num + CODE_RATE), F::from(1))];
    for i in 1..variable_num + 1 {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    interpolate_cosets
}

pub fn deepfold_commit<const LOG_2_PATH_LENGTH: usize, F: FftField>(
    polynomial: MultilinearPolynomial<F>,
    random_oracle: &RandomOracle<F>,
    hypercube_interpolation: Vec<F>,
) -> (Prover<F>, Commit<F>) {
    use util::STEP;
    debug_assert_eq!(polynomial.variable_num(), 2 + LOG_2_PATH_LENGTH);
    let variable_num = polynomial.variable_num();
    let prover = Prover::new_with_interpolation(
        variable_num,
        &interpolate_cosets(variable_num),
        polynomial,
        random_oracle,
        STEP,
        hypercube_interpolation,
    );
    let commit = prover.commit_polynomial();
    (prover, commit)
}

/// Evaluate a cubic polynomial as defined by four evaluation points, at 0, 1,
/// 2 and 3.
fn evaluate_poly_from_evals<F: Field>(evals: [F; 4], eval_at: F) -> F {
    // First compute the coefficients
    let mut coefficients = [F::ZERO; 4];
    coefficients[0] = evals[0];
    coefficients[1] = evals
        .iter()
        .zip(
            [
                -F::from(11) / F::from(6),
                F::from(3),
                -F::from(3) / F::from(2),
                F::ONE / F::from(3),
            ]
            .iter(),
        )
        .map(|(a, b)| *a * b)
        .sum();
    coefficients[2] = evals
        .iter()
        .zip(
            [
                F::ONE,
                -F::from(5) / F::from(2),
                F::from(2),
                -F::ONE / F::from(2),
            ]
            .iter(),
        )
        .map(|(a, b)| *a * b)
        .sum();
    coefficients[3] = evals
        .iter()
        .zip(
            [
                -F::ONE / F::from(6),
                F::ONE / F::from(2),
                -F::ONE / F::from(2),
                F::ONE / F::from(6),
            ]
            .iter(),
        )
        .map(|(a, b)| *a * b)
        .sum();

    coefficients[0]
        + eval_at * (coefficients[1] + eval_at * (coefficients[2] + eval_at * coefficients[3]))
}

pub fn compute_mask_value<F: Field>(coefficients: &[F], point: &[F], r_y: F) -> F {
    let (coef_0, coef_1) = coefficients.split_at(coefficients.len() / 2);
    point
        .iter()
        .zip(coef_0.chunks_exact(3))
        .map(|(var, values)| {
            let [c, b, a] = values.try_into().unwrap();
            ((a * var + b) * var + c) * var
        })
        .sum::<F>()
        + r_y
            * point
                .iter()
                .zip(coef_1.chunks_exact(3))
                .map(|(var, values)| {
                    let [c, b, a] = values.try_into().unwrap();
                    ((a * var + b) * var + c) * var
                })
                .sum::<F>()
}

pub fn verify_pok<
    const VARIABLE_COUNT: usize,
    const PATH_LENGTH: usize,
    const PATH_LENGTH_TIMES_TWO: usize,
    const PATH_LENGTH_TIMES_FOUR: usize,
    const LOG_2_PATH_LENGTH: usize,
    const LOG_2_PATH_LENGTH_PLUS_ONE: usize,
    const Q_VARIABLE_COUNT: usize,
    const COMMITMENT_SIZE: usize,
    const FINAL_ROUND_EVALUATIONS: usize,
    F: FftField,
    PK: PublicKey<F>,
>(
    public_key: &PK,
    proof: Proof<F>,
    random_oracle: &RandomOracle<F>,
    mask_check_mode: MaskCheckMode,
) -> Result<bool, Box<dyn Error>> {
    assert_eq!(VARIABLE_COUNT, LOG_2_PATH_LENGTH * 2 + 1);

    let mut arthur =
        pok_io_pattern::<LOG_2_PATH_LENGTH, COMMITMENT_SIZE, FINAL_ROUND_EVALUATIONS, F, PK>(
            mask_check_mode,
        )
        .to_verifier_state(&proof.narg);

    // Read message
    let message = arthur.next_bytes()?;
    if message != MESSAGE {
        return Ok(false);
    }

    let commitment: Commit<F> = SerializableCommit::<F>::deserialize_compressed(
        arthur.next_bytes::<COMMITMENT_SIZE>()?.as_slice(),
    )?
    .into();

    let [e_0, e_1]: [F; 2] = arthur.challenge_scalars()?;
    let k: [F; LOG_2_PATH_LENGTH] = arthur.challenge_scalars()?;
    let pk_randomness = PK::challenge_randomness(&mut arthur)?;

    let a_a_i = P::<
        VARIABLE_COUNT,
        LOG_2_PATH_LENGTH_PLUS_ONE,
        LOG_2_PATH_LENGTH,
        PATH_LENGTH,
        PATH_LENGTH_TIMES_TWO,
        Q_VARIABLE_COUNT,
        F,
        PK,
    >::a_a_i(e_1);
    let a_a_j = P::<
        VARIABLE_COUNT,
        LOG_2_PATH_LENGTH_PLUS_ONE,
        LOG_2_PATH_LENGTH,
        PATH_LENGTH,
        PATH_LENGTH_TIMES_TWO,
        Q_VARIABLE_COUNT,
        F,
        PK,
    >::a_a_j(e_1);
    let a_c_j = P::<
        VARIABLE_COUNT,
        LOG_2_PATH_LENGTH_PLUS_ONE,
        LOG_2_PATH_LENGTH,
        PATH_LENGTH,
        PATH_LENGTH_TIMES_TWO,
        Q_VARIABLE_COUNT,
        F,
        PK,
    >::a_c_j(e_1);

    // Do sumcheck
    let mut sum = F::ZERO;
    let mut evaluation_point: [F; VARIABLE_COUNT] = [F::ZERO; VARIABLE_COUNT];
    for evaluation in evaluation_point.iter_mut() {
        let [eval_1, eval_2, eval_3] = arthur.next_scalars()?;
        // Now interpolate
        *evaluation = arthur.challenge_scalars::<1>()?[0];
        sum = evaluate_poly_from_evals([sum - eval_1, eval_1, eval_2, eval_3], *evaluation);
    }
    let x = evaluation_point[0];
    let i = &evaluation_point[1..1 + LOG_2_PATH_LENGTH];
    let j = &evaluation_point[1 + LOG_2_PATH_LENGTH..1 + 2 * LOG_2_PATH_LENGTH];

    // We arrive at our final claim. At the point `evaluation_point`,
    // polynomial `P` must equal `sum`. Now we need to verify this.

    // Obtain the polynomial openings
    let [b_0_i, b_0_j, b_1_i, b_1_j, mask]: [F; 5] = match mask_check_mode {
        MaskCheckMode::InsidePCS => arthur.next_scalars(),
        MaskCheckMode::Additional => {
            let [mask_0, mask_1] = arthur.next_scalars()?;
            let [r_y]: [F; 1] = arthur.challenge_scalars()?;
            let mask_coefficients: [F; FINAL_ROUND_EVALUATIONS] = arthur.next_scalars()?;
            let mask_value = compute_mask_value(
                &mask_coefficients[..FINAL_ROUND_EVALUATIONS - 4],
                &evaluation_point,
                r_y,
            );
            let _mask_constant_term = mask_0 + r_y * mask_1 - mask_value;
            Ok([
                mask_coefficients[FINAL_ROUND_EVALUATIONS - 4],
                mask_coefficients[FINAL_ROUND_EVALUATIONS - 3],
                mask_coefficients[FINAL_ROUND_EVALUATIONS - 2],
                mask_coefficients[FINAL_ROUND_EVALUATIONS - 1],
                mask_0,
            ])
        }
    }?;

    let start_openings = PK::Checker::read_openings(&mut arthur)?;
    let end_openings = PK::Checker::read_openings(&mut arthur)?;

    let mut eval_point = vec![evaluation_point[0], F::ZERO];
    eval_point.extend_from_slice(&evaluation_point[1..1 + LOG_2_PATH_LENGTH]);

    let mut verifier = deepfold::verifier::Verifier::new(
        &mut StdRng::from_entropy(),
        LOG_2_PATH_LENGTH + 2,
        &interpolate_cosets(LOG_2_PATH_LENGTH + 2),
        commitment,
        random_oracle,
        util::STEP,
    );
    verifier.set_open_point(&eval_point);
    let commitment_verifies = verifier.verify(proof.final_evaluation_proof);

    let q_value = Q::<0, F>::direct_eval(i, j, &k);
    let p_value = x
        * e_0
        * q_value
        * ((a_a_j * b_0_j + a_a_i * b_0_i) * (b_0_j - b_0_i) - (b_1_i + a_c_j * b_1_j))
        + public_key.direct_evaluate_constraints(
            x,
            &pk_randomness,
            &[start_openings, end_openings],
        )
        + mask;
    let piop_matches = p_value == sum;

    Ok(piop_matches && commitment_verifies)
}

pub fn interpolate_polynomial<F: Field>(evaluations: &[F]) -> MultilinearPolynomial<F> {
    let mut coefficients: Vec<F> = vec![F::ZERO; evaluations.len()];
    for (index, evaluation) in evaluations.iter().enumerate() {
        // Add eq(i,j) * evaluations[j] to the coefficients.
        //
        // We can interpret index as a sequence of bits, specifying the
        // coordinate. If the ith bit is 0, we need to multiply by (1 - x_i),
        // if it is 1, we need to multiply by x_i.
        //
        // In other words, we should always add 0, -1 or 1 times the
        // evaluation, depending on the i th bits of the evaluation point and
        // coefficient.
        // - Evaluation point bit is 0, and coef bit is 0, 1.
        // - Evaluation point bit is 0, and coef bit is 1, -1.
        // - Evaluation point bit is 1, and coef bit is 0, 0.
        // - Evaluation point bit is 1, and coef bit is 1, 1.
        for (coef_index, coefficient) in coefficients.iter_mut().enumerate() {
            // If there is any bit where evaluation is 1, and coef is 0, set to 0.
            if (!coef_index) & index != 0 {
                continue;
            }

            let mut value = *evaluation;
            // For every bit where evaluation is 0, and coef is 1, invert.
            if ((!index) & coef_index).count_ones() % 2 == 1 {
                value = -value;
            }

            *coefficient += value;
        }
    }
    MultilinearPolynomial::new(coefficients)
}

#[cfg(test)]
mod tests {
    use util::algebra::field::arkfield::Fp2256;

    use crate::evaluate_poly_from_evals;

    #[test]
    fn test_interpolation() {
        let evaluations = [
            Fp2256::from(4),
            Fp2256::from(4),
            Fp2256::from(4),
            Fp2256::from(4),
        ];
        assert_eq!(
            evaluate_poly_from_evals(evaluations, Fp2256::from(100)),
            Fp2256::from(4)
        );
    }
}
