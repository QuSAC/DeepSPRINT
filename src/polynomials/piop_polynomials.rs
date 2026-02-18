//! These are the polynomials that are custom for our PIOP.
use std::ops::MulAssign;

use crate::polynomials::HypercubePoint;

use super::{HypercubeEvalPoly, PiopPolynomial};
use ark_ff::Field;

mod q;

pub use q::Q;

/// The mask. Internally, the coefficients are represented in lexographical order,
/// so
/// `[var0 deg1, var0 deg2, var0 deg3, var1 deg1]`
/// the last variable will be taken as the constant term.
#[derive(Clone)]
pub struct Mask<'a, const VARIABLE_COUNT: usize, F: Field> {
    b_0_mask: &'a [F],
    b_1_mask: &'a [F],
    /// The current total of all evaluated points. Initially this will just be
    /// the constant term.
    current_total: F,
    variable_count: usize,
}

impl<'a, const VARIABLE_COUNT: usize, F: Field>
    Mask<'a, VARIABLE_COUNT, F>
{
    pub fn new(b_0_mask: &'a [F], b_1_mask: &'a [F]) -> Self {
        // The mask should have 3 coefficients per variable, plus a constant
        // term. This number is then doubled to mask the coefficient opening.
        debug_assert!(b_0_mask.len() == 3 * VARIABLE_COUNT + 1);
        debug_assert!(b_1_mask.len() == 3 * VARIABLE_COUNT + 1);
        let mask = Self {
            b_0_mask,
            b_1_mask,
            current_total: b_0_mask[0],
            variable_count: VARIABLE_COUNT,
        };
        // Check that the constant term was set correctly
        debug_assert!((mask.eval_field_then_sum_hypercube(F::ZERO)
            + mask.eval_field_then_sum_hypercube(F::ONE))
        .is_zero());
        mask
    }

    /// The index of the first free variable. The index of the last free
    /// variable is always `VARIABLE_COUNT - 1`.
    fn current_variable_index(&self) -> usize {
        VARIABLE_COUNT - self.variable_count
    }

    fn eval_field_then_sum_hypercube(&self, point: F) -> F {
        // Evaluate for the current variable
        let starting_index_first = 1 + self.current_variable_index() * 3;
        let starting_index_next = starting_index_first + 3;
        let remaining_variable_count = self.variable_count - 1;
        let [c, b, a] = self.b_0_mask[starting_index_first..starting_index_next]
            .try_into()
            .unwrap();
        (self.current_total + ((a * point + b) * point + c) * point
        // Add the remaining variables. Precisely half the points on the
        // hypercube will give a contribution to any sum.
        + self.b_0_mask[starting_index_next..].iter()
            // Take three values per variable, excluding the variable we
            // already have.
            .take(3 * remaining_variable_count)
            .sum::<F>() / F::from(2))
            * F::from(1 << remaining_variable_count)
    }

    /// Compute the constant term that would fix the coefficients.
    pub fn fix_constant_term(coefficients: &mut [F]) {
        debug_assert!(coefficients.len() == 1 + 3 * VARIABLE_COUNT);
        coefficients[0] = -coefficients[1..].iter().sum::<F>() / F::from(2);
    }

    /// Adds a variable y to mask the mask's coefficients.
    pub fn eval_with_y_eq_1(&self, point: &[F]) -> F {
        self.b_1_mask[0]
        + point
            .iter()
            .zip(self.b_1_mask[1..].chunks_exact(3))
            .map(|(var, values)| {
                let [c, b, a] = values.try_into().unwrap();
                ((a * var + b) * var + c) * var
            })
            .sum::<F>()
    }

    /// Get the masked (non-constant) coefficients.
    pub fn masked_coefficients(&self, r_y: F) -> Vec<F> {
        self.b_0_mask[1..1 + VARIABLE_COUNT * 3]
            .iter()
            .zip(self.b_1_mask[1..1 + VARIABLE_COUNT * 3].iter())
            .map(|(a, b)| *a + r_y * b)
            .collect()
    }
}

impl<'a, const VARIABLE_COUNT: usize, F: Field> PiopPolynomial<F>
    for Mask<'a, VARIABLE_COUNT, F>
{
    fn eval(&self, point: &[F]) -> F {
        let starting_index = 1 + self.current_variable_index() * 3;
        self.current_total
            + point
                .iter()
                .zip(self.b_0_mask[starting_index..].chunks_exact(3))
                .map(|(var, values)| {
                    let [c, b, a] = values.try_into().unwrap();
                    ((a * var + b) * var + c) * var
                })
                .sum::<F>()
    }

    fn fix_variable(&mut self, var: F) {
        let starting_index = 1 + self.current_variable_index() * 3;
        self.variable_count -= 1;
        // Compute a*x^2 + b*x + c
        if let [c, b, a] = self.b_0_mask[starting_index..starting_index + 3] {
            self.current_total += ((a * var + b) * var + c) * var;
        } else {
            panic!("out of bounds!");
        }
    }

    fn eval_hypercube(&self, point: super::HypercubePoint) -> F {
        let starting_index = 1 + self.current_variable_index() * 3;
        self.current_total
            + point
                .into_bool_iter(self.variable_count)
                .zip(self.b_0_mask[starting_index..].chunks_exact(3))
                .flat_map(|(is_one, values)| values.iter().filter(move |_| is_one))
                .sum::<F>()
    }

    fn variable_count(&self) -> usize {
        self.variable_count
    }

    fn eval_field_then_hypercube(&self, f: F, hypercube_point: HypercubePoint) -> F {
        let starting_index = 1 + self.current_variable_index() * 3;
        let [c, b, a] = self.b_0_mask[starting_index..starting_index + 3]
            .try_into()
            .unwrap();

        self.current_total
            // Evaluate first variable as normal
            + ((a * f + b) * f + c) * f
            // Evaluate the remainder as a hypercube point
            + hypercube_point
                .into_bool_iter(self.variable_count - 1)
                .zip(self.b_0_mask[starting_index + 3..].chunks_exact(3))
                .flat_map(|(is_one, values)| values.iter().filter(move |_| is_one) )
                .sum::<F>()
    }
}

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
    >
{
    pub fn new(
        e: [F; 2],
        k: [F; LOG_2_PATH_LENGTH],
        // The evaluations of A, including a second half of masking
        b_0: &'a HypercubeEvalPoly<F>,
        // The evaluations of B, including a second half of masking
        b_1: &'a HypercubeEvalPoly<F>,
        q_construction_mode: QConstructionMode,
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
        let mask_0_borrowed = &b_0.evals()[..1 + VARIABLE_COUNT * 3];
        let mask_1_borrowed = &b_1.evals()[..1 + VARIABLE_COUNT * 3];
        let mask = Mask::new(mask_0_borrowed, mask_1_borrowed);

        let mut b_1_j = b_1.clone();
        // We multiply this value in the polynomial. We do need to store it for
        // later to compensate.
        b_1_j *= a_c_j;

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
        }
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
            (0..(1 << LOG_2_PATH_LENGTH))
                .map(HypercubePoint)
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
        } else if self.b_0_i.variable_count() > 0 {
            // Subtract 1 due to the variable we will assign now
            let remaining_i_vars = self.b_0_i.variable_count() - 1;
            let i_hypercube_variables_mask = (1 << remaining_i_vars) - 1;
            // Evaluate i
            // We can sum over j, and determine the qualifying i from this
            (0..(1 << LOG_2_PATH_LENGTH))
                .map(HypercubePoint)
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
        } else {
            // Evaluate j
            debug_assert!(self.variable_count() > 0);
            debug_assert_eq!(self.b_0_i.variable_count(), 0);
            debug_assert_eq!(self.b_1_i.variable_count(), 0);
            debug_assert!(self.b_0_j.variable_count() > 0);
            debug_assert!(self.b_1_j.variable_count() > 0);
            let variables_j = self.b_1_j.variable_count() - 1;
            (0usize..(1 << variables_j))
                .map(HypercubePoint)
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
                * ((self.a_a_j * b_0_j + self.a_a_i * b_0_i) * (b_0_j - b_0_i) - (b_1_i + b_1_j))
                + self.mask.eval(point)
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
        } else {
            unimplemented!()
        }
    }
}

#[cfg(test)]
mod tests {
    use std::array;

    use super::*;
    use ark_ff::{AdditiveGroup, Field};
    use util::algebra::field::arkfield::Fp2256;

    #[test]
    fn test_p_eval_hypercube_vs_field() {
        let e = [Fp2256::from(2), Fp2256::from(3)];
        let k = [Fp2256::from(25235); 8];
        let mut b_0_evals = [Fp2256::from(32423); 2 * (1 << 8)];
        Mask::<17, _>::fix_constant_term(&mut b_0_evals[..1 + 3 * 17]);
        Mask::<17, _>::new(&b_0_evals[..1 + 3 * 17], &b_0_evals[..1 + 3 * 17]);
        let b_1_evals = [Fp2256::from(32423); 2 * (1 << 8)];
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + 8);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + 8);
        let mut p: P<17, 9, 8, { 1 << 8 }, { 2 * (1 << 8) }, 16, Fp2256> =
            P::new(e, k, &b_0, &b_1, QConstructionMode::Naive);
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
    fn test_mask_eval_hypercube() {
        const ONE: Fp2256 = Fp2256::ONE;
        const MASK_VARIABLE_COUNT: usize = 3 * 17 + 1;
        let mut val = array::from_fn::<_, { MASK_VARIABLE_COUNT }, _>(|i| Fp2256::from(i as u32));
        Mask::<17, _>::fix_constant_term(&mut val);
        let mask: Mask<17, _> = Mask::new(&val, &val);
        assert_eq!(
            mask.eval_hypercube(HypercubePoint(0b1_0110_1101_0010_1010)),
            mask.eval_field_then_hypercube(ONE, HypercubePoint(0b0110_1101_0010_1010))
        );
    }

    #[test]
    fn test_mask_new() {
        const MASK_COEFFICIENT_COUNT: usize = 3 * 17 + 1;
        let mut val: [Fp2256; MASK_COEFFICIENT_COUNT] = [Fp2256::from(32423); MASK_COEFFICIENT_COUNT];
        Mask::<17, _>::fix_constant_term(&mut val);
        Mask::<17, _>::new(&val, &val);
    }

    #[test]
    fn test_mask_eval_hypercube_field() {
        const ONE: Fp2256 = Fp2256::ONE;
        const ZERO: Fp2256 = Fp2256::ZERO;
        const MASK_COEFFICIENT_COUNT: usize = 3 * 17 + 1;
        let mut val: [Fp2256; MASK_COEFFICIENT_COUNT] = array::from_fn(|i| Fp2256::from(i as u32));
        Mask::<17, _>::fix_constant_term(&mut val);
        let mask: Mask<17, _> = Mask::new(&val, &val);
        assert_eq!(
            mask.eval_hypercube(HypercubePoint(0b1_0110_1101_0010_1010)),
            mask.eval(&[
                ONE, ZERO, ONE, ONE, // 1011
                ZERO, ONE, ONE, ZERO, // 0110
                ONE,  // 1
                ZERO, ZERO, ONE, ZERO, // 0010
                ONE, ZERO, ONE, ZERO // 1010
            ])
        );
    }

    #[test]
    fn test_mask_eval_sum_hypercube_vs_field() {
        const MASK_VARIABLE_COUNT: usize = 17;
        let mut val: [_; 3 * MASK_VARIABLE_COUNT + 1] = array::from_fn(|i| Fp2256::from(i as u32));
        Mask::<17, _>::fix_constant_term(&mut val);
        let mask: Mask<17, _> = Mask::new(&val, &val);
        let assignment = Fp2256::from(2352);
        let sum = (0..(1 << 16))
            .map(HypercubePoint)
            .map(|point| mask.eval_field_then_hypercube(assignment, point))
            .sum::<Fp2256>();
        assert_eq!(mask.eval_field_then_sum_hypercube(assignment), sum);
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
        Mask::<VARIABLE_COUNT, _>::fix_constant_term(&mut b_0_evals[..1 + 3 * VARIABLE_COUNT]);
        let b_1_evals = [Fp2256::from(32423); 2 * (1 << LOG_2_PATH_LENGTH)];
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + LOG_2_PATH_LENGTH);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + LOG_2_PATH_LENGTH);
        let p: P<
            VARIABLE_COUNT,
            { LOG_2_PATH_LENGTH + 1 },
            LOG_2_PATH_LENGTH,
            PATH_LENGTH,
            PATH_LENGTH_TIMES_2,
            Q_VARIABLE_COUNT,
            Fp2256,
        > = P::new(e, k, &b_0, &b_1, QConstructionMode::Naive);
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
        Mask::<17, _>::fix_constant_term(&mut b_0_evals[..1 + 3 * 17]);
        let b_1_evals = [Fp2256::from(1); 2 * (1 << 8)];
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + 8);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + 8);
        let p: P<17, 9, 8, { 1 << 8 }, { 2 * (1 << 8) }, 16, Fp2256> =
            P::new(e, k, &b_0, &b_1, QConstructionMode::Naive);
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
        Mask::<17, _>::fix_constant_term(&mut b_0_evals[..1 + 3 * 17]);
        let b_1_evals: [_; PATH_LENGTH_TIMES_2] =
            array::from_fn(|i| Fp2256::from(2412 * i as u64 + 42));
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + LOG_2_PATH_LENGTH);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + LOG_2_PATH_LENGTH);
        let p: P<
            VARIABLE_COUNT,
            { LOG_2_PATH_LENGTH + 1 },
            LOG_2_PATH_LENGTH,
            PATH_LENGTH,
            PATH_LENGTH_TIMES_2,
            Q_VARIABLE_COUNT,
            Fp2256,
        > = P::new(e, k, &b_0, &b_1, QConstructionMode::Naive);
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
        Mask::<VARIABLE_COUNT, _>::fix_constant_term(&mut b_0_evals[..1 + 3 * VARIABLE_COUNT]);
        let b_1_evals: [_; PATH_LENGTH_TIMES_2] =
            array::from_fn(|i| Fp2256::from(2342 * i as u64 + 42));
        let b_0 = HypercubeEvalPoly::new(&b_0_evals, 1 + LOG_2_PATH_LENGTH);
        let b_1 = HypercubeEvalPoly::new(&b_1_evals, 1 + LOG_2_PATH_LENGTH);
        let mut p: P<
            VARIABLE_COUNT,
            { LOG_2_PATH_LENGTH + 1 },
            LOG_2_PATH_LENGTH,
            PATH_LENGTH,
            PATH_LENGTH_TIMES_2,
            Q_VARIABLE_COUNT,
            Fp2256,
        > = P::new(e, k, &b_0, &b_1, QConstructionMode::Naive);
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

    #[test]
    fn test_mask_fix() {
        const VARIABLES: usize = 2 * 8 + 1;
        let mut mask_evals = [Fp2256::ZERO; {1 + 3 * 17}];
        for (a, i) in mask_evals.iter_mut().enumerate() {
            *i = Fp2256::from(a as u64);
        }
        Mask::<VARIABLES, _>::fix_constant_term(&mut mask_evals);
        let mut mask = Mask::<VARIABLES, _>::new(&mask_evals, &mask_evals);

        let assignment = array::from_fn::<_, VARIABLES, _>(|i| Fp2256::from(i as u64));
        let a = mask.eval(&assignment);

        mask.fix_variable(assignment[0]);
        assert_eq!(mask.eval(&assignment[1..]), a);
        mask.fix_variable(assignment[1]);
        assert_eq!(mask.eval(&assignment[2..]), a);
        mask.fix_variable(assignment[2]);
        assert_eq!(mask.eval(&assignment[3..]), a);
    }
}
