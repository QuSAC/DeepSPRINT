//! These are the polynomials that are custom for our PIOP.
#[cfg(feature = "rayon")]
use rayon::prelude::*;
use std::{marker::PhantomData, ops::MulAssign};

use crate::{
    polynomials::{piop_polynomials::j_invariant_checking::JInvariantChecker, HypercubePoint},
    PublicKey,
};

use super::{HypercubeEvalPoly, PiopPolynomial};
use ark_ff::Field;

pub mod j_invariant_checking;
mod mask;
mod q;

pub use mask::Mask;
pub use q::Q;
use spongefish::{codecs::arkworks_algebra::FieldToUnitDeserialize, ProofResult, VerifierState};

#[derive(Clone)]
pub struct P<
    'a,
    const VARIABLE_COUNT: usize,
    const LOG_2_PATH_LENGTH_PLUS_ONE: usize,
    const LOG_2_PATH_LENGTH: usize,
    const PATH_LENGTH: usize,
    const PATH_LENGTH_TIMES_2: usize,
    const Q_VARIABLE_COUNT: usize,
    F: Field,
    PK: PublicKey<F>,
> {
    e_0: F,
    a_a_j: F,
    a_a_i: F,
    a_c_j: F,
    q: Q<PATH_LENGTH, F>,
    b_0_i: HypercubeEvalPoly<F>,
    b_0_j: HypercubeEvalPoly<F>,
    b_1_i: HypercubeEvalPoly<F>,
    b_1_j: HypercubeEvalPoly<F>,
    variable_count: usize,
    mask: Mask<'a, VARIABLE_COUNT, F>,
    j_invariant_checkers: [PK::Checker; 2],
}

pub trait Checker<F: Field>: Sync {
    fn eval_at_x(&self, x: F) -> F;
    fn fix_variable(&mut self, x: F);
    fn value(&self) -> F;
    fn eval_then_sum_hypercube(&self, val: F) -> F;
    const OPENING_LEN: usize;
    type Openings: AsRef<[F]>;
    fn get_openings(&self) -> Self::Openings;
    fn read_openings(arthur: &mut VerifierState<'_>) -> ProofResult<Self::Openings>;
}
#[derive(Default)]
pub struct ACChecker<F: Field> {
    _a: PhantomData<F>,
}
impl<F: Field> Checker<F> for ACChecker<F> {
    const OPENING_LEN: usize = 0;
    type Openings = [F; 0];
    fn get_openings(&self) -> Self::Openings {
        []
    }
    fn eval_at_x(&self, _x: F) -> F {
        F::ZERO
    }
    fn fix_variable(&mut self, _x: F) {}
    fn value(&self) -> F {
        F::ZERO
    }
    fn eval_then_sum_hypercube(&self, _val: F) -> F {
        F::ZERO
    }
    fn read_openings(_arthur: &mut VerifierState<'_>) -> ProofResult<Self::Openings> {
        Ok([])
    }
}
impl<F: Field> Checker<F> for JInvariantChecker<F> {
    const OPENING_LEN: usize = 7;
    type Openings = [F; 7];
    fn get_openings(&self) -> Self::Openings {
        self.get_openings()
    }
    fn eval_at_x(&self, x: F) -> F {
        self.eval_at_x(x)
    }
    fn fix_variable(&mut self, x: F) {
        self.fix_variable(x);
    }
    fn value(&self) -> F {
        self.value()
    }
    fn eval_then_sum_hypercube(&self, val: F) -> F {
        self.eval_then_sum_hypercube(val)
    }
    fn read_openings(arthur: &mut VerifierState<'_>) -> ProofResult<Self::Openings> {
        arthur.next_scalars()
    }
}

impl<F: Field> MulAssign<F> for HypercubeEvalPoly<F> {
    fn mul_assign(&mut self, rhs: F) {
        for a in self.evals.iter_mut() {
            *a *= rhs;
        }
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum QConstructionMode {
    Naive,
    GrayCodes,
}

impl<
        'a,
        const VARIABLE_COUNT: usize,
        const LOG_2_PATH_LENGTH_PLUS_ONE: usize,
        const LOG_2_PATH_LENGTH: usize,
        const PATH_LENGTH: usize,
        const PATH_LENGTH_TIMES_2: usize,
        const Q_VARIABLE_COUNT: usize,
        F: Field,
        PK: PublicKey<F>,
    >
    P<
        'a,
        VARIABLE_COUNT,
        LOG_2_PATH_LENGTH_PLUS_ONE,
        LOG_2_PATH_LENGTH,
        PATH_LENGTH,
        PATH_LENGTH_TIMES_2,
        Q_VARIABLE_COUNT,
        F,
        PK,
    >
{
    pub fn new(
        e: [F; 2],
        k: [F; LOG_2_PATH_LENGTH],
        // The evaluations of A, including a second half of masking
        b_0: &'a HypercubeEvalPoly<F>,
        // The evaluations of B, including a second half of masking
        b_1: &'a HypercubeEvalPoly<F>,
        statement_checker_randomness: &PK::CheckingRandomness,
        q_construction_mode: QConstructionMode,
        public_key: &PK,
    ) -> Self {
        debug_assert_eq!(PATH_LENGTH_TIMES_2, 2 * (1 << LOG_2_PATH_LENGTH));
        debug_assert_eq!(2 * LOG_2_PATH_LENGTH + 1, VARIABLE_COUNT);
        debug_assert_eq!(2 * LOG_2_PATH_LENGTH, Q_VARIABLE_COUNT);
        debug_assert_eq!(LOG_2_PATH_LENGTH + 1, LOG_2_PATH_LENGTH_PLUS_ONE);
        let e_0 = e[0];
        let a_a_j = Self::a_a_j(e[1]);
        let a_a_i = Self::a_a_i(e[1]);
        let a_c_j = Self::a_c_j(e[1]);

        let q = match q_construction_mode {
            QConstructionMode::GrayCodes => Q::<PATH_LENGTH, F>::new_gray_codes(&k),
            QConstructionMode::Naive => Q::<PATH_LENGTH, F>::new(&k),
        };

        // Construct mask
        let mask = Mask::new_from_evals(&b_0.evals, &b_1.evals);

        let mut b_1_j = b_1.clone();
        // We multiply this value in the polynomial. We do need to store it for
        // later to compensate.
        b_1_j *= a_c_j;

        let j_invariant_checkers = public_key.get_checkers(
            LOG_2_PATH_LENGTH,
            b_0.evals(),
            b_1.evals(),
            statement_checker_randomness,
        );

        Self {
            e_0,
            a_a_i,
            a_a_j,
            a_c_j,
            q,
            b_0_i: b_0.clone(),
            b_0_j: b_0.clone(),
            b_1_i: b_1.clone(),
            b_1_j,
            variable_count: VARIABLE_COUNT,
            mask,
            j_invariant_checkers,
        }
    }

    pub fn checkers(&self) -> &[PK::Checker; 2] {
        &self.j_invariant_checkers
    }

    pub fn mask(&self) -> &Mask<'a, VARIABLE_COUNT, F> {
        &self.mask
    }

    pub fn a_a_j(e_1: F) -> F {
        HypercubeEvalPoly::new(&[F::ZERO, F::ONE / F::from(36)], 1).eval(&[e_1])
    }

    pub fn a_a_i(e_1: F) -> F {
        HypercubeEvalPoly::new(&[-F::ONE / F::from(12), -F::ONE / F::from(36)], 1).eval(&[e_1])
    }

    pub fn a_c_j(e_1: F) -> F {
        HypercubeEvalPoly::new(&[-F::ONE / F::from(8), F::ZERO], 1).eval(&[e_1])
    }

    /// Tests if the constraint system is satisfied. Only works if no variables
    /// have been fixed!
    pub fn is_system_satisfied(&self) -> bool {
        assert_eq!(self.variable_count(), VARIABLE_COUNT);
        (self.eval_field_then_sum_hypercube(F::ZERO) + self.eval_field_then_sum_hypercube(F::ONE))
            .is_zero()
    }

    pub fn eval_field_then_sum_hypercube(&self, f: F) -> F {
        // No variables have been fixed yet
        if self.variable_count() == VARIABLE_COUNT {
            let x = f;
            let i_hypercube_variables_mask = (1 << LOG_2_PATH_LENGTH) - 1;
            #[cfg(feature = "rayon")]
            let iter = (0..(1 << LOG_2_PATH_LENGTH)).into_par_iter();
            #[cfg(not(feature = "rayon"))]
            let iter = 0..(1 << LOG_2_PATH_LENGTH);
            iter.map(HypercubePoint)
                .map(|j| {
                    // We can compute i, since we already know that Q will be zero otherwise
                    let i_hypercube_remainder =
                        HypercubePoint((j.0.wrapping_sub(1)) & i_hypercube_variables_mask);
                    let b_0_i = self
                        .b_0_i
                        .eval_field_then_hypercube(x, i_hypercube_remainder);
                    let b_1_i = self
                        .b_1_i
                        .eval_field_then_hypercube(x, i_hypercube_remainder);
                    let b_0_j = self.b_0_j.eval_field_then_hypercube(x, j);
                    let b_1_j = self.b_1_j.eval_field_then_hypercube(x, j);
                    let ij = HypercubePoint((i_hypercube_remainder.0 << LOG_2_PATH_LENGTH) | j.0);
                    self.q.eval_hypercube(ij)
                        * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i)
                            - (b_1_i + b_1_j))
                })
                .sum::<F>()
                * self.e_0
                * x
                + self.mask.eval_field_then_sum_hypercube(f)
                // Compute the j-invariant checkers
                + self
                    .j_invariant_checkers
                    .as_ref()
                    .iter().map(|c| c.eval_then_sum_hypercube(f)).sum::<F>()
        } else if self.b_0_i.variable_count() > 0 {
            // Subtract 1 due to the variable we will assign now
            let remaining_i_vars = self.b_0_i.variable_count() - 1;
            let i_hypercube_variables_mask = (1 << remaining_i_vars) - 1;
            // Evaluate i
            // We can sum over j, and determine the qualifying i from this
            #[cfg(feature = "rayon")]
            let iter = (0..(1 << LOG_2_PATH_LENGTH)).into_par_iter();
            #[cfg(not(feature = "rayon"))]
            let iter = 0..(1 << LOG_2_PATH_LENGTH);
            iter.map(HypercubePoint)
                .map(|j| {
                    // Now compute j, the first couple of variables are free,
                    // the others are determined by i.
                    let i_hypercube_remainder =
                        HypercubePoint(j.0.wrapping_sub(1) & i_hypercube_variables_mask);
                    let b_0_i = self
                        .b_0_i
                        .eval_field_then_hypercube(f, i_hypercube_remainder);
                    let b_1_i = self
                        .b_1_i
                        .eval_field_then_hypercube(f, i_hypercube_remainder);
                    let b_0_j = self.b_0_j.eval_hypercube(j);
                    let b_1_j = self.b_1_j.eval_hypercube(j);
                    let ij = HypercubePoint(i_hypercube_remainder.0 << LOG_2_PATH_LENGTH | j.0);
                    self.q.eval_field_then_hypercube(f, ij)
                        * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i)
                            - (b_1_i + b_1_j))
                })
                .sum::<F>()
                * self.e_0 // x is included in e_0
                + self.mask.eval_field_then_sum_hypercube(f)
                + self
                    .j_invariant_checkers
                    .as_ref().iter().map(|c| c.eval_then_sum_hypercube(f)).sum::<F>()
        } else {
            // Evaluate j
            debug_assert!(self.variable_count() > 0);
            debug_assert_eq!(self.b_0_i.variable_count(), 0);
            debug_assert_eq!(self.b_1_i.variable_count(), 0);
            debug_assert!(self.b_0_j.variable_count() > 0);
            debug_assert!(self.b_1_j.variable_count() > 0);
            let variables_j = self.b_1_j.variable_count() - 1;
            #[cfg(feature = "rayon")]
            let iter = (0..(1 << variables_j)).into_par_iter();
            #[cfg(not(feature = "rayon"))]
            let iter = 0..(1 << variables_j);
            iter.map(HypercubePoint)
                .map(|j_hypercube_remainder| {
                    let b_0_i = self.b_0_i.get_as_const();
                    let b_1_i = self.b_1_i.get_as_const();
                    let b_0_j = self
                        .b_0_j
                        .eval_field_then_hypercube(f, j_hypercube_remainder);
                    let b_1_j = self
                        .b_1_j
                        .eval_field_then_hypercube(f, j_hypercube_remainder);
                    let ij = j_hypercube_remainder;
                    self.q.eval_field_then_hypercube(f, ij)
                        * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i)
                            - (b_1_i + b_1_j))
                })
                .sum::<F>()
                * self.e_0 // x is included in e_0
                + self.mask.eval_field_then_sum_hypercube(f)
                + self
                    .j_invariant_checkers
                    .as_ref().iter().map(|c| c.eval_then_sum_hypercube(f)).sum::<F>()
        }
    }

    pub fn final_evaluations(&self) -> [F; 5] {
        [
            self.b_0_i.get_as_const(),
            self.b_0_j.get_as_const(),
            self.b_1_i.get_as_const(),
            self.b_1_j.get_as_const() / self.a_c_j,
            self.mask.eval(&[]),
        ]
    }
}

impl<
        'a,
        const VARIABLE_COUNT: usize,
        const LOG_2_PATH_LENGTH_PLUS_ONE: usize,
        const LOG_2_PATH_LENGTH: usize,
        const PATH_LENGTH: usize,
        const PATH_LENGTH_TIMES_2: usize,
        const Q_VARIABLE_COUNT: usize,
        F: Field,
        PK: PublicKey<F>,
    > PiopPolynomial<F>
    for P<
        'a,
        VARIABLE_COUNT,
        LOG_2_PATH_LENGTH_PLUS_ONE,
        LOG_2_PATH_LENGTH,
        PATH_LENGTH,
        PATH_LENGTH_TIMES_2,
        Q_VARIABLE_COUNT,
        F,
        PK,
    >
{
    fn eval(&self, point: &[F]) -> F {
        let x;
        let q_eval_point;
        let b_0_i;
        let b_0_j;
        let b_1_i;
        let b_1_j;
        if self.variable_count() == VARIABLE_COUNT {
            x = point[0];
            let mut i = [x; LOG_2_PATH_LENGTH_PLUS_ONE];
            let mut j = [x; LOG_2_PATH_LENGTH_PLUS_ONE];
            i[1..].copy_from_slice(&point[1..(LOG_2_PATH_LENGTH + 1)]);
            j[1..].copy_from_slice(&point[(LOG_2_PATH_LENGTH + 1)..(2 * LOG_2_PATH_LENGTH + 1)]);
            b_0_i = self.b_0_i.eval(&i);
            b_0_j = self.b_0_j.eval(&j);
            b_1_i = self.b_1_i.eval(&i);
            b_1_j = self.b_1_j.eval(&j);
            q_eval_point = &point[1..];
            self.e_0
                * x
                * self.q.eval(q_eval_point)
                * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i) - (b_1_i + b_1_j))
                + self.mask.eval(point)
                + self
                    .j_invariant_checkers
                    .as_ref()
                    .iter()
                    .map(|c| c.eval_at_x(point[0]))
                    .sum::<F>()
        } else {
            let j_len = LOG_2_PATH_LENGTH.min(point.len());
            let i_end = point.len() - j_len;
            let i = &point[0..i_end];
            let j = &point[i_end..];
            b_0_i = self.b_0_i.eval(i);
            b_0_j = self.b_0_j.eval(j);
            b_1_i = self.b_1_i.eval(i);
            b_1_j = self.b_1_j.eval(j);
            q_eval_point = point;
            self.e_0
                // x has been absorbed in e_0
                * self.q.eval(q_eval_point)
                * (self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i)
                - (b_1_i + b_1_j)
                + self.mask.eval(point)
                + self
                    .j_invariant_checkers
                    .as_ref()
                    .iter()
                    .map(|c| c.value())
                    .sum::<F>()
        }
    }

    fn eval_hypercube(&self, point: super::HypercubePoint) -> F {
        debug_assert!(point.0 < (1 << self.variable_count));
        if self.variable_count() < VARIABLE_COUNT {
            let (i, j) = point.split(LOG_2_PATH_LENGTH, LOG_2_PATH_LENGTH);
            let b_0_i = self.b_0_i.eval_hypercube(i);
            let b_0_j = self.b_0_j.eval_hypercube(j);
            let b_1_i = self.b_1_i.eval_hypercube(i);
            let b_1_j = self.b_1_j.eval_hypercube(j);
            self.e_0
                * self.q.eval_hypercube(point)
                * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i) - (b_1_i + b_1_j))
                + self.mask.eval_hypercube(point)
                + self
                    .j_invariant_checkers
                    .as_ref()
                    .iter()
                    .map(|c| c.value())
                    .sum::<F>()
        } else {
            let (x, ij) = point.split(1, 2 * LOG_2_PATH_LENGTH);
            let (i, j) = ij.split(LOG_2_PATH_LENGTH, LOG_2_PATH_LENGTH);
            let xi = HypercubePoint(x.0 << LOG_2_PATH_LENGTH | i.0);
            let xj = HypercubePoint(x.0 << LOG_2_PATH_LENGTH | j.0);
            let b_0_i = self.b_0_i.eval_hypercube(xi);
            let b_0_j = self.b_0_j.eval_hypercube(xj);
            let b_1_i = self.b_1_i.eval_hypercube(xi);
            let b_1_j = self.b_1_j.eval_hypercube(xj);
            let x = x.value_at_pos(1, 0);
            if x {
                self.e_0
                    * self.q.eval_hypercube(ij)
                    * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i)
                        - (b_1_i + b_1_j))
                    + self.mask.eval_hypercube(point)
            } else {
                F::ZERO
            }
        }
    }

    fn fix_variable(&mut self, var: F) {
        if self.variable_count == VARIABLE_COUNT {
            // Fix x first
            self.b_0_i.fix_variable(var);
            self.b_1_i.fix_variable(var);
            self.b_0_j.fix_variable(var);
            self.b_1_j.fix_variable(var);
            self.mask.fix_variable(var);
            // Multiply x into e_0, to save this multiplication on every later iteration
            self.e_0 *= var;
            debug_assert_eq!(self.b_0_i.variable_count(), LOG_2_PATH_LENGTH);
            debug_assert_eq!(self.b_0_j.variable_count(), LOG_2_PATH_LENGTH);
            debug_assert_eq!(self.b_1_i.variable_count(), LOG_2_PATH_LENGTH);
            debug_assert_eq!(self.b_1_j.variable_count(), LOG_2_PATH_LENGTH);
            debug_assert_eq!(self.mask.variable_count(), 2 * LOG_2_PATH_LENGTH);
        } else {
            self.q.fix_variable(var);
            self.mask.fix_variable(var);
            if self.variable_count > LOG_2_PATH_LENGTH {
                debug_assert_eq!(
                    self.b_0_i.variable_count(),
                    self.variable_count - LOG_2_PATH_LENGTH
                );
                debug_assert_eq!(
                    self.b_1_i.variable_count(),
                    self.variable_count - LOG_2_PATH_LENGTH
                );
                debug_assert_eq!(self.b_0_j.variable_count(), LOG_2_PATH_LENGTH);
                debug_assert_eq!(self.b_1_j.variable_count(), LOG_2_PATH_LENGTH);
                self.b_0_i.fix_variable(var);
                self.b_1_i.fix_variable(var);
            } else {
                debug_assert_eq!(self.b_0_i.variable_count(), 0);
                debug_assert_eq!(self.b_1_i.variable_count(), 0);
                debug_assert_eq!(self.b_0_j.variable_count(), self.variable_count());
                debug_assert_eq!(self.b_1_j.variable_count(), self.variable_count());
                self.b_0_j.fix_variable(var);
                self.b_1_j.fix_variable(var);
            }
        }
        for j_checker in self.j_invariant_checkers.iter_mut() {
            j_checker.fix_variable(var);
        }
        self.variable_count -= 1;
        debug_assert_eq!(
            self.variable_count(),
            self.b_0_i.variable_count() + self.b_0_j.variable_count()
        );
        debug_assert_eq!(
            self.variable_count(),
            self.b_1_i.variable_count() + self.b_1_j.variable_count()
        );
        debug_assert_eq!(self.variable_count(), self.mask.variable_count());
        debug_assert_eq!(self.variable_count(), self.q.variable_count());
    }

    fn variable_count(&self) -> usize {
        self.variable_count
    }

    fn eval_field_then_hypercube(&self, f: F, hypercube_point: HypercubePoint) -> F {
        let x;
        if self.variable_count() == VARIABLE_COUNT {
            // No variables have been fixed yet
            debug_assert_eq!(self.variable_count, VARIABLE_COUNT);
            x = f;
            let (i, j) = hypercube_point.split(LOG_2_PATH_LENGTH, LOG_2_PATH_LENGTH);
            let b_1_i = self.b_1_i.eval_field_then_hypercube(x, i);
            let b_0_i = self.b_0_i.eval_field_then_hypercube(x, i);
            let b_0_j = self.b_0_j.eval_field_then_hypercube(x, j);
            let b_1_j = self.b_1_j.eval_field_then_hypercube(x, j);
            self.e_0
                * x
                * self.q.eval_hypercube(hypercube_point)
                * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i) - (b_1_i + b_1_j))
                + self.mask.eval_field_then_hypercube(x, hypercube_point)
                + self
                    .j_invariant_checkers
                    .as_ref()
                    .iter()
                    .map(|c| c.eval_at_x(f))
                    .sum::<F>()
        } else {
            unimplemented!()
        }
    }
}

#[cfg(test)]
mod tests {
    use std::array;

    use crate::keys::RadicalPublicKey;

    use super::*;
    use ark_ff::{AdditiveGroup, Field};
    use util::algebra::field::arkfield::Fp2256;

    #[test]
    fn test_p_eval_hypercube_vs_field() {
        let e = [Fp2256::from(2), Fp2256::from(3)];
        let k = [Fp2256::from(25235); 8];
        let mut b_0_evals = [Fp2256::from(32423); 2 * (1 << 8)];
        let mut b_1_evals = [Fp2256::from(32423); 2 * (1 << 8)];
        Mask::<17, _>::fix_constant_term(&mut b_0_evals, &mut b_1_evals);
        Mask::<17, _>::new_from_evals(&b_0_evals, &b_1_evals);
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + 8);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + 8);
        let public_key =
            RadicalPublicKey::from_a_c(*b_0_evals.last().unwrap(), *b_1_evals.last().unwrap());
        let mut p: P<17, 9, 8, { 1 << 8 }, { 2 * (1 << 8) }, 16, Fp2256, RadicalPublicKey<_>> =
            P::new(e, k, &b_0, &b_1, &(), QConstructionMode::Naive, &public_key);
        let point_hypercube = HypercubePoint::from_usize(0b1_10111101_10110101);
        let mut iter = point_hypercube.into_point_iter(17);
        let point_field: [Fp2256; 17] = array::from_fn(|_| iter.next().unwrap());

        let eval = p.eval_hypercube(point_hypercube);
        assert_eq!(p.eval(&point_field), eval);

        // Now make the first field element something else
        let f = Fp2256::from(12412);
        let point_hypercube = HypercubePoint::from_usize(0b10111101_10110101);
        let eval = p.eval_field_then_hypercube(f, point_hypercube);
        p.fix_variable(f);
        assert_eq!(p.eval_hypercube(point_hypercube), eval);
    }

    #[test]
    fn test_p_eval_sum_hypercube_vs_field() {
        const LOG_2_PATH_LENGTH: usize = 8;
        const VARIABLE_COUNT: usize = 2 * LOG_2_PATH_LENGTH + 1;
        const PATH_LENGTH: usize = 1 << LOG_2_PATH_LENGTH;
        const PATH_LENGTH_TIMES_2: usize = PATH_LENGTH * 2;
        const Q_VARIABLE_COUNT: usize = 2 * LOG_2_PATH_LENGTH;
        let e = [Fp2256::from(2), Fp2256::from(3)];
        let k = [Fp2256::from(25235); LOG_2_PATH_LENGTH];
        let mut b_0_evals = [Fp2256::from(32423); 2 * (1 << LOG_2_PATH_LENGTH)];
        let mut b_1_evals = [Fp2256::from(32423); 2 * (1 << LOG_2_PATH_LENGTH)];
        Mask::<VARIABLE_COUNT, _>::fix_constant_term(&mut b_0_evals, &mut b_1_evals);
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + LOG_2_PATH_LENGTH);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + LOG_2_PATH_LENGTH);
        let public_key =
            RadicalPublicKey::from_a_c(*b_0_evals.last().unwrap(), *b_1_evals.last().unwrap());
        let p: P<
            VARIABLE_COUNT,
            { LOG_2_PATH_LENGTH + 1 },
            LOG_2_PATH_LENGTH,
            PATH_LENGTH,
            PATH_LENGTH_TIMES_2,
            Q_VARIABLE_COUNT,
            Fp2256,
            RadicalPublicKey<_>,
        > = P::new(e, k, &b_0, &b_1, &(), QConstructionMode::Naive, &public_key);
        let assignment = Fp2256::from(2352);
        let sum = (0..(1 << (VARIABLE_COUNT - 1)))
            .map(HypercubePoint)
            .map(|point| p.eval_field_then_hypercube(assignment, point))
            .sum::<Fp2256>();
        assert_eq!(p.eval_field_then_sum_hypercube(assignment), sum);
    }

    #[test]
    fn test_p_eval_sum_hypercube_vs_field_2() {
        let e = [Fp2256::from(1), Fp2256::from(1)];
        let k = [Fp2256::from(2); 8];
        let mut b_0_evals = [Fp2256::from(1); 2 * (1 << 8)];
        let mut b_1_evals = [Fp2256::from(1); 2 * (1 << 8)];
        Mask::<17, _>::fix_constant_term(&mut b_0_evals, &mut b_1_evals);
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + 8);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + 8);
        let public_key =
            RadicalPublicKey::from_a_c(*b_0_evals.last().unwrap(), *b_1_evals.last().unwrap());
        let p: P<17, 9, 8, { 1 << 8 }, { 2 * (1 << 8) }, 16, Fp2256, RadicalPublicKey<_>> =
            P::new(e, k, &b_0, &b_1, &(), QConstructionMode::Naive, &public_key);
        let assignment = Fp2256::from(1);
        let sum = (0..(1 << 16))
            .map(HypercubePoint)
            .map(|point| p.eval_field_then_hypercube(assignment, point))
            .sum::<Fp2256>();
        assert_eq!(p.eval_field_then_sum_hypercube(assignment), sum);
    }

    #[test]
    fn test_p_eval_sum_hypercube_vs_field_3() {
        const LOG_2_PATH_LENGTH: usize = 8;
        const VARIABLE_COUNT: usize = 2 * LOG_2_PATH_LENGTH + 1;
        const PATH_LENGTH: usize = 1 << LOG_2_PATH_LENGTH;
        const PATH_LENGTH_TIMES_2: usize = PATH_LENGTH * 2;
        const Q_VARIABLE_COUNT: usize = 2 * LOG_2_PATH_LENGTH;
        let e = [Fp2256::from(2), Fp2256::from(3)];
        let k: [_; LOG_2_PATH_LENGTH] = array::from_fn(|i| Fp2256::from(24 * i as u64 + 23));
        let mut b_0_evals: [_; PATH_LENGTH_TIMES_2] =
            array::from_fn(|i| Fp2256::from(3 * i as u64 + 3412));
        let mut b_1_evals: [_; PATH_LENGTH_TIMES_2] =
            array::from_fn(|i| Fp2256::from(2412 * i as u64 + 42));
        Mask::<17, _>::fix_constant_term(&mut b_0_evals, &mut b_1_evals);
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + LOG_2_PATH_LENGTH);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + LOG_2_PATH_LENGTH);
        let public_key =
            RadicalPublicKey::from_a_c(*b_0_evals.last().unwrap(), *b_1_evals.last().unwrap());
        let p: P<
            VARIABLE_COUNT,
            { LOG_2_PATH_LENGTH + 1 },
            LOG_2_PATH_LENGTH,
            PATH_LENGTH,
            PATH_LENGTH_TIMES_2,
            Q_VARIABLE_COUNT,
            Fp2256,
            RadicalPublicKey<_>,
        > = P::new(e, k, &b_0, &b_1, &(), QConstructionMode::Naive, &public_key);
        let assignment = Fp2256::from(2352);
        let sum = (0..(1 << (VARIABLE_COUNT - 1)))
            .map(HypercubePoint)
            .map(|point| p.eval_field_then_hypercube(assignment, point))
            .sum::<Fp2256>();
        assert_eq!(p.eval_field_then_sum_hypercube(assignment), sum);
    }

    #[test]
    fn test_p_eval_sum_hypercube_fix_variable() {
        const LOG_2_PATH_LENGTH: usize = 8;
        const VARIABLE_COUNT: usize = 2 * LOG_2_PATH_LENGTH + 1;
        const PATH_LENGTH: usize = 1 << LOG_2_PATH_LENGTH;
        const PATH_LENGTH_TIMES_2: usize = PATH_LENGTH * 2;
        const Q_VARIABLE_COUNT: usize = 2 * LOG_2_PATH_LENGTH;
        let e = [Fp2256::from(2), Fp2256::from(3)];
        let k: [_; LOG_2_PATH_LENGTH] = array::from_fn(|i| Fp2256::from(8 * i as u64 + 23));
        let mut b_0_evals: [_; PATH_LENGTH_TIMES_2] =
            array::from_fn(|i| Fp2256::from(7 * i as u64 + 3412));
        let mut b_1_evals: [_; PATH_LENGTH_TIMES_2] =
            array::from_fn(|i| Fp2256::from(2342 * i as u64 + 42));
        Mask::<VARIABLE_COUNT, _>::fix_constant_term(&mut b_0_evals, &mut b_1_evals);
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + LOG_2_PATH_LENGTH);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + LOG_2_PATH_LENGTH);
        let public_key =
            RadicalPublicKey::from_a_c(*b_0_evals.last().unwrap(), *b_1_evals.last().unwrap());
        let mut p: P<
            VARIABLE_COUNT,
            { LOG_2_PATH_LENGTH + 1 },
            LOG_2_PATH_LENGTH,
            PATH_LENGTH,
            PATH_LENGTH_TIMES_2,
            Q_VARIABLE_COUNT,
            Fp2256,
            RadicalPublicKey<_>,
        > = P::new(e, k, &b_0, &b_1, &(), QConstructionMode::Naive, &public_key);
        let assignment = Fp2256::from(2352);
        let value = p.eval_field_then_sum_hypercube(assignment);
        p.fix_variable(assignment);
        assert_eq!(
            value,
            p.eval_field_then_sum_hypercube(Fp2256::ZERO)
                + p.eval_field_then_sum_hypercube(Fp2256::ONE)
        );
    }

    #[test]
    fn test_eval() {
        let q = Q::<2, Fp2256>::new(&[Fp2256::ZERO]);
        // Should return eq(i, k) * eq(j, k + 1)
        assert_eq!(q.eval(&[Fp2256::ZERO, Fp2256::ZERO]), Fp2256::ZERO);
        assert_eq!(q.eval(&[Fp2256::ZERO, Fp2256::ONE]), Fp2256::ONE);
        assert_eq!(q.eval(&[Fp2256::ONE, Fp2256::ZERO]), Fp2256::ZERO);
        assert_eq!(q.eval(&[Fp2256::ONE, Fp2256::ONE]), Fp2256::ZERO);

        let q2 = Q::<4, _>::new(&[Fp2256::ZERO, Fp2256::ONE]);
        // x = 1
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ZERO, Fp2256::ZERO, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ZERO, Fp2256::ZERO, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ZERO, Fp2256::ONE, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ZERO, Fp2256::ONE, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ONE, Fp2256::ZERO, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ONE, Fp2256::ZERO, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ONE, Fp2256::ONE, Fp2256::ZERO]),
            Fp2256::ONE
        );
        assert_eq!(
            q2.eval(&[Fp2256::ZERO, Fp2256::ONE, Fp2256::ONE, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ZERO, Fp2256::ZERO, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ZERO, Fp2256::ZERO, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ZERO, Fp2256::ONE, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ZERO, Fp2256::ONE, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ONE, Fp2256::ZERO, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ONE, Fp2256::ZERO, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ONE, Fp2256::ONE, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            q2.eval(&[Fp2256::ONE, Fp2256::ONE, Fp2256::ONE, Fp2256::ONE]),
            Fp2256::ZERO
        );
    }
}
