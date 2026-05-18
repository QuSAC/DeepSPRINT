use ark_ff::Field;
use ark_ff::UniformRand;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use bitvec::array::BitArray;
use spongefish::codecs::arkworks_algebra::UnitToField;
use spongefish::VerifierState;
use spongefish::{ProofError, ProverState};

use crate::polynomials::piop_polynomials::j_invariant_checking::JInvariantChecker;
use crate::polynomials::piop_polynomials::j_invariant_checking::LinearPolynomial;
use crate::polynomials::piop_polynomials::ACChecker;
use crate::polynomials::piop_polynomials::Checker;

#[derive(Eq, PartialEq, Clone, Copy)]
pub enum OneOrMinusOne {
    MinusOne,
    One,
}

impl OneOrMinusOne {
    fn into_field<F: Field>(self) -> F {
        match self {
            OneOrMinusOne::One => F::from(1),
            OneOrMinusOne::MinusOne => -F::from(1),
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq)]
pub struct RadicalPublicKey<F: Field> {
    end_a: F,
    end_c: F,
}

impl<F: Field> RadicalPublicKey<F> {
    pub fn from_a_c(end_a: F, end_c: F) -> Self {
        Self { end_a, end_c }
    }
}

/// A version of the public key with just a j-invariant.
#[derive(CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq)]
pub struct JInvariantPublicKey<F: Field>(F);

pub trait PublicKey<F: Field>: Sized {
    type Checker: Checker<F>;
    type CheckingRandomness;
    const CHECKING_RANDOMNESS_SIZE: usize;
    fn get_checking_randomness(
        prover_state: &mut ProverState,
    ) -> Result<Self::CheckingRandomness, ProofError>;
    fn challenge_randomness(
        verifier_state: &mut VerifierState,
    ) -> Result<Self::CheckingRandomness, ProofError>;
    fn fix_checker_relation(&self, log_2_path_size: usize, evals_0: &mut [F], evals_1: &mut [F]);
    fn get_checkers(
        &self,
        log_2_path_size: usize,
        evals_0: &[F],
        evals_1: &[F],
        randomness: &Self::CheckingRandomness,
    ) -> [Self::Checker; 2];
    fn direct_evaluate_constraints(
        &self,
        x: F,
        randomness: &Self::CheckingRandomness,
        openings: &[<Self::Checker as Checker<F>>::Openings; 2],
    ) -> F;
}
impl<F: Field> PublicKey<F> for RadicalPublicKey<F> {
    type Checker = ACChecker<F>;
    type CheckingRandomness = ();
    const CHECKING_RANDOMNESS_SIZE: usize = 0;
    fn challenge_randomness(
        _verifier_state: &mut VerifierState,
    ) -> Result<Self::CheckingRandomness, ProofError> {
        Ok(())
    }
    fn get_checking_randomness(
        _prover_state: &mut ProverState,
    ) -> Result<Self::CheckingRandomness, ProofError> {
        Ok(())
    }
    fn fix_checker_relation(
        &self,
        _log_2_path_size: usize,
        _evals_0: &mut [F],
        _evals_1: &mut [F],
    ) {
    }
    fn get_checkers(
        &self,
        _log_2_path_size: usize,
        _evals_0: &[F],
        _evals_1: &[F],
        _randomness: &(),
    ) -> [Self::Checker; 2] {
        [ACChecker::default(), ACChecker::default()]
    }
    fn direct_evaluate_constraints(
        &self,
        _x: F,
        _randomness: &Self::CheckingRandomness,
        _openings: &[<Self::Checker as Checker<F>>::Openings; 2],
    ) -> F {
        F::ZERO
    }
}
impl<F: Field> PublicKey<F> for JInvariantPublicKey<F> {
    type Checker = JInvariantChecker<F>;
    type CheckingRandomness = [F; 2];
    const CHECKING_RANDOMNESS_SIZE: usize = 1;
    fn challenge_randomness(
        verifier_state: &mut VerifierState,
    ) -> Result<Self::CheckingRandomness, ProofError> {
        verifier_state.challenge_scalars()
    }
    fn get_checking_randomness(
        prover_state: &mut ProverState,
    ) -> Result<Self::CheckingRandomness, ProofError> {
        prover_state.challenge_scalars()
    }
    fn fix_checker_relation(&self, log_2_path_size: usize, evals_0: &mut [F], evals_1: &mut [F]) {
        // 1 maps to A and C, 0 maps to the mask.
        let half_mask_size = evals_0.len() / 2;
        debug_assert_eq!(2_usize.pow(log_2_path_size as u32), half_mask_size);
        debug_assert_eq!(2_usize.pow(log_2_path_size as u32), half_mask_size);
        let start_index = 6 * log_2_path_size + 4;
        assert!(start_index + 10 < half_mask_size);

        let a_start = evals_0[half_mask_size];
        let c_start = evals_1[half_mask_size];
        let a_end = evals_0[2 * half_mask_size - 1];
        let c_end = evals_1[2 * half_mask_size - 1];

        let d_1_start = &mut evals_1[start_index + 1];
        *d_1_start = a_start.pow([2]);
        let d_1_start = *d_1_start;
        let d_1_end = &mut evals_1[start_index + 6];
        *d_1_end = a_end.pow([2]);
        let d_1_end = *d_1_end;

        let d_2_start = &mut evals_1[start_index + 2];
        *d_2_start = c_start.pow([2]);
        let d_2_end = &mut evals_1[start_index + 7];
        *d_2_end = c_end.pow([2]);

        let d_3_start = &mut evals_1[start_index + 3];
        *d_3_start = d_1_start.pow([2]);
        let d_3_end = &mut evals_1[start_index + 8];
        *d_3_end = d_1_end.pow([2]);

        let d_4_start = &mut evals_1[start_index + 4];
        *d_4_start = c_start.inverse().expect("c cannot be 0");
        let d_4_end = &mut evals_1[start_index + 9];
        *d_4_end = c_end.inverse().expect("c cannot be 0");

        let d_5_start = &mut evals_1[start_index + 5];
        *d_5_start = (F::from(4) * c_start - d_1_start)
            .inverse()
            .expect("denominator non-zero");
        let d_5_end = &mut evals_1[start_index + 10];
        *d_5_end = (F::from(4) * c_end - d_1_end)
            .inverse()
            .expect("denominator non-zero");
    }

    fn get_checkers(
        &self,
        log_2_path_size: usize,
        evals_0: &[F],
        evals_1: &[F],
        randomness: &[F; 2],
    ) -> [Self::Checker; 2] {
        // 1 maps to A and C, 0 maps to the mask.
        let half_mask_size = evals_0.len() / 2;
        debug_assert_eq!(2_usize.pow(log_2_path_size as u32), half_mask_size);
        debug_assert_eq!(2_usize.pow(log_2_path_size as u32), half_mask_size);
        let start_index = 6 * log_2_path_size + 4;
        let d_1_start =
            LinearPolynomial::from_evals(evals_0[start_index + 1], evals_1[start_index + 1]);
        let d_2_start =
            LinearPolynomial::from_evals(evals_0[start_index + 2], evals_1[start_index + 2]);
        let d_3_start =
            LinearPolynomial::from_evals(evals_0[start_index + 3], evals_1[start_index + 3]);
        let d_4_start =
            LinearPolynomial::from_evals(evals_0[start_index + 4], evals_1[start_index + 4]);
        let d_5_start =
            LinearPolynomial::from_evals(evals_0[start_index + 5], evals_1[start_index + 5]);
        let d_1_end =
            LinearPolynomial::from_evals(evals_0[start_index + 6], evals_1[start_index + 6]);
        let d_2_end =
            LinearPolynomial::from_evals(evals_0[start_index + 7], evals_1[start_index + 7]);
        let d_3_end =
            LinearPolynomial::from_evals(evals_0[start_index + 8], evals_1[start_index + 8]);
        let d_4_end =
            LinearPolynomial::from_evals(evals_0[start_index + 9], evals_1[start_index + 9]);
        let d_5_end =
            LinearPolynomial::from_evals(evals_0[start_index + 10], evals_1[start_index + 10]);
        let a_start = LinearPolynomial::from_evals(evals_0[0], evals_0[half_mask_size]);
        let c_start = LinearPolynomial::from_evals(evals_1[0], evals_1[half_mask_size]);
        let a_end = LinearPolynomial::from_evals(
            evals_0[half_mask_size - 1],
            evals_0[2 * half_mask_size - 1],
        );
        let c_end = LinearPolynomial::from_evals(
            evals_1[half_mask_size - 1],
            evals_1[2 * half_mask_size - 1],
        );
        let num_variables = 2 * log_2_path_size + 1;
        let checker_1 = JInvariantChecker::new(
            randomness[0],
            starting_j_invariant(),
            a_start,
            c_start,
            d_1_start,
            d_2_start,
            d_3_start,
            d_4_start,
            d_5_start,
            num_variables
        );
        debug_assert_eq!(compute_j_invariant(a_end.eval(F::ONE), c_end.eval(F::ONE)), self.0);
        let checker_2 = JInvariantChecker::new(
            randomness[1],
            self.0,
            a_end,
            c_end,
            d_1_end,
            d_2_end,
            d_3_end,
            d_4_end,
            d_5_end,
            num_variables
        );
        [checker_1, checker_2]
    }

    fn direct_evaluate_constraints(
        &self,
        x: F,
        &[h1, h2]: &Self::CheckingRandomness,
        &[opening1, opening2]: &[<Self::Checker as Checker<F>>::Openings; 2],
    ) -> F {
        JInvariantChecker::eval_with_opening(x, h1, starting_j_invariant(), &opening1)
            + JInvariantChecker::eval_with_opening(x, h2, self.0, &opening2)
    }
}

pub fn compute_j_invariant<F: Field>(a: F, c: F) -> F {
    let a_2 = a.pow([2]);
    F::from(256) * (-a_2 + F::from(3) * c).pow([3]) / (F::from(4) * c - a_2) / c / c
}

impl<F: Field> From<RadicalPublicKey<F>> for JInvariantPublicKey<F> {
    fn from(value: RadicalPublicKey<F>) -> Self {
        let RadicalPublicKey { end_a, end_c } = value;
        JInvariantPublicKey(compute_j_invariant(end_a, end_c))
    }
}

pub struct CompressedPrivateKey<const PATH_LENGTH_DIV_64: usize>(
    BitArray<[u64; PATH_LENGTH_DIV_64]>,
);

impl<const PATH_LENGTH_DIV_64: usize> CompressedPrivateKey<PATH_LENGTH_DIV_64> {
    /// Generate a new private key.
    pub fn rand(mut rng: impl Rng) -> Self {
        let mut bits = BitArray::new([0u64; PATH_LENGTH_DIV_64]);
        for data in bits.data.iter_mut() {
            *data = u64::rand(&mut rng);
        }
        Self(bits)
    }
}

pub struct PrivateKey<F: Field> {
    path_a: Vec<F>,
    path_c: Vec<F>,
}

impl<F: Field> PrivateKey<F> {
    pub fn len(&self) -> usize {
        debug_assert_eq!(self.path_a.len(), self.path_c.len());
        self.path_a.len()
    }

    pub fn path_a(&self) -> &[F] {
        &self.path_a
    }

    pub fn path_c(&self) -> &[F] {
        &self.path_c
    }

    /// Checks if the private key is well-formed
    pub fn check(&self, public_key: &RadicalPublicKey<F>) -> bool {
        if !(self.path_a.last().unwrap() == &public_key.end_a
            && self.path_c.last().unwrap() == &public_key.end_c)
        {
            return false;
        }

        for i in 0..self.path_a.len() - 1 {
            if (-self.path_a[i] / F::from(12)) * (self.path_a[i + 1] - self.path_a[i])
                != self.path_c[i] - self.path_c[i + 1] / F::from(8)
            {
                return false;
            }
            if (self.path_a[i + 1] - self.path_a[i]).square() / F::from(36) != self.path_c[i] {
                return false;
            }
        }

        true
    }
}

/// Expand a private key. This will create the entire path from the directions
/// and from the starting curve.
pub fn expand_keys<const PATH_LENGTH_DIV_64: usize, F: Field>(
    compressed_private_key: &CompressedPrivateKey<PATH_LENGTH_DIV_64>,
) -> (RadicalPublicKey<F>, PrivateKey<F>) {
    let key: Vec<OneOrMinusOne> = compressed_private_key
        .0
        .iter()
        .map(|b| {
            if *b {
                OneOrMinusOne::One
            } else {
                OneOrMinusOne::MinusOne
            }
        })
        // We use 255 steps
        .skip(1)
        .collect();
    let (path_a, path_c) = build_path(&key);
    let end_a = *path_a.last().unwrap();
    let end_c = *path_c.last().unwrap();
    (
        RadicalPublicKey { end_a, end_c },
        PrivateKey { path_a, path_c },
    )
}

#[inline]
pub fn starting_curve<F: Field>() -> (F, F) {
    (F::from(0), F::from(251948161))
}

#[inline]
pub fn starting_j_invariant<F: Field>() -> F {
    let (a_start, c_start) = starting_curve();
    compute_j_invariant(a_start, c_start)
}

#[inline]
pub fn build_path<F: Field>(key: &[OneOrMinusOne]) -> (Vec<F>, Vec<F>) {
    let mut a: Vec<F> = Vec::with_capacity(key.len());
    let mut c: Vec<F> = Vec::with_capacity(key.len());
    let (a_start, c_start) = starting_curve();
    a.push(a_start);
    c.push(c_start);

    let mut a_curr = a_start;
    let mut c_curr = c_start;
    for k in key.iter() {
        let c_sqrt = k.into_field::<F>()
            * (c_curr
                .sqrt()
                .expect("could not follow path - C_i is not a square residue"));
        c_curr = (c_curr + c_curr + c_sqrt * a_curr) * F::from(4);
        a_curr = a_curr + F::from(6) * c_sqrt;
        a.push(a_curr);
        c.push(c_curr);
    }

    (a, c)
}

/// Pad a private key by length `length`.
pub fn pad_walk<F: Field>(length: usize) -> (Vec<F>, Vec<F>) {
    let mut a = Vec::with_capacity(length);
    let mut c = Vec::with_capacity(length);
    a.push(F::from(4));
    c.push(F::from(4));

    let mut a_curr = a[0];
    let mut c_curr = c[0];
    for _ in 0..length {
        let mut c_sqrt = c_curr
            .sqrt()
            .expect("could not follow path - C_i is not a square residue");
        c_sqrt.double_in_place();
        let c_sqrt_6 = c_sqrt + c_sqrt + c_sqrt;
        let a_next = c_sqrt_6 + a_curr;
        c_curr.double_in_place();
        c_curr.double_in_place();
        let mut c_next = c_sqrt * a_curr + c_curr;
        c_next.double_in_place();
        a.push(a_next);
        c.push(c_next);
        a_curr = a_next;
        c_curr = c_next;
    }

    (a, c)
}
