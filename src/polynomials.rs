//! This module contains polynomials custom built for the PIOP.
//!
//! As a matter of convention, assignments [F0, F1, F2, F3], are indexed most
//! significant value first.

pub mod piop_polynomials;

use ark_ff::Field;
use std::{fmt::Debug, iter, marker::PhantomData};

#[derive(Eq, PartialEq, Clone, Copy)]
pub struct HypercubePoint(usize);

impl Debug for HypercubePoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#b}", self.0)
    }
}

impl HypercubePoint {
    pub fn from_usize(u: usize) -> Self {
        Self(u)
    }

    /// Iterate through this point, starting with the most significant bit.
    #[cfg(test)]
    fn into_point_iter<F: Field>(self, variable_count: usize) -> impl Iterator<Item = F> {
        debug_assert!(self.0 < (1 << variable_count));
        iter::successors(Some(1 << (variable_count - 1)), |mask| Some(mask >> 1))
            .take(variable_count)
            .map(move |mask| if self.0 & mask == 0 { F::ZERO } else { F::ONE })
    }

    /// Split this point into two. The first element of the tuple will be a
    /// point of size `first_len`, and the second part will be of size
    /// `second_len`.
    ///
    /// # Panics
    /// This function may panic if the point is larger than
    /// `first_len + second_len`.
    fn split(self, first_len: usize, second_len: usize) -> (HypercubePoint, HypercubePoint) {
        debug_assert!(self.0 < (1 << (first_len + second_len)));
        (
            HypercubePoint(self.0 >> second_len),
            HypercubePoint(self.0 & ((1 << second_len) - 1)),
        )
    }

    fn take_last_n(self, n: usize) -> HypercubePoint {
        HypercubePoint(self.0 & ((1 << n) - 1))
    }

    /// Get the value of this hypercube point at position `pos`. Indexing
    /// starts with the most significant bit in the representation.
    fn value_at_pos(self, num_variables: usize, pos: usize) -> bool {
        ((self.0 >> (num_variables - pos - 1)) & 1) != 0
    }

    pub fn into_bool_iter(self, num_variables: usize) -> impl Iterator<Item = bool> {
        (0..num_variables).map(move |pos| self.value_at_pos(num_variables, pos))
    }
}

pub trait PiopPolynomial<F: Field> {
    fn variable_count(&self) -> usize;
    fn eval(&self, point: &[F]) -> F;
    fn eval_hypercube(&self, point: HypercubePoint) -> F;
    fn fix_variable(&mut self, var: F);
    fn fix_all_variables(&mut self, vars: &[F]) -> F {
        debug_assert_eq!(vars.len(), self.variable_count());
        for var in vars {
            self.fix_variable(*var);
        }
        self.eval_hypercube(HypercubePoint(0))
    }
    fn eval_field_then_hypercube(&self, f: F, hypercube_point: HypercubePoint) -> F;
}

/// Implements the multilinear polynomial C * eq(x, p), where p is a fixed
/// point on the boolean hypercube and C is some fixed coefficient.
///
/// Evaluation on the boolean hypercube takes O(1) time, evaluation outside
/// this domain takes time O(`NUM_VARS`).
#[derive(PartialEq, Eq, Clone, Copy)]
pub struct EqFixedPoint<F: Field> {
    point: HypercubePoint,
    phantom_data: PhantomData<F>,
    num_variables: usize,
    coefficient: F,
}

impl<F: Field> EqFixedPoint<F> {
    pub fn new(point: HypercubePoint, num_variables: usize, coefficient: F) -> Self {
        debug_assert!(point.0 < (1 << num_variables));
        Self {
            point,
            phantom_data: PhantomData,
            num_variables,
            coefficient,
        }
    }
}

impl<F: Field> PiopPolynomial<F> for EqFixedPoint<F> {
    fn variable_count(&self) -> usize {
        self.num_variables
    }

    fn eval(&self, point: &[F]) -> F {
        iter::successors(Some(self.point.0), |u| Some(u >> 1))
            .map(|p| p & 1)
            .zip(point.iter().rev())
            .map(|(p, value)| if p == 1 { *value } else { F::ONE - *value })
            .product::<F>()
            * self.coefficient
    }

    fn eval_hypercube(&self, point: HypercubePoint) -> F {
        if point == self.point {
            self.coefficient
        } else {
            F::ZERO
        }
    }

    fn fix_variable(&mut self, var: F) {
        if self.point.value_at_pos(self.num_variables, 0) {
            self.coefficient *= var;
        } else {
            // c * (1 - v) = c - c*v
            self.coefficient -= self.coefficient * var;
        }
        self.num_variables -= 1;
        // We won't update the point, since the higher part will just be ignored
    }

    fn eval_field_then_hypercube(&self, var: F, hypercube_point: HypercubePoint) -> F {
        let first_variable_mask = 1 << (self.num_variables - 1);
        let remainder_mask = first_variable_mask - 1;
        if self.point.0 & remainder_mask == hypercube_point.0 {
            if (self.point.0 & first_variable_mask) != 0 {
                self.coefficient * var
            } else {
                self.coefficient - self.coefficient * var
            }
        } else {
            F::ZERO
        }
    }
}

/// Construct a polynomial from hypercube evaluations.
///
/// Evaluation on a point on the hypercube takes O(1) time. Evaluation on
/// another point takes time O(`NUM_VARS * 2^NUM_VARS`) time.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct HypercubeEvalPoly<F: Field> {
    evals: Vec<F>,
    variable_count: usize,
}

impl<F: Field> HypercubeEvalPoly<F> {
    pub fn new(evals: &[F], variable_count: usize) -> Self {
        debug_assert_eq!(1 << variable_count, evals.len());
        Self {
            evals: evals.to_vec(),
            variable_count,
        }
    }

    pub fn evals(&self) -> &'_ [F] {
        &self.evals
    }

    pub fn evals_mut(&mut self) -> &'_ mut [F] {
        &mut self.evals
    }

    pub fn get_as_const(&self) -> F {
        debug_assert_eq!(self.variable_count(), 0);
        self.evals[0]
    }
}

impl<F: Field> PiopPolynomial<F> for HypercubeEvalPoly<F> {
    fn eval(&self, point: &[F]) -> F {
        debug_assert_eq!(point.len(), self.variable_count);
        self.evals
            .iter()
            .enumerate()
            .take(1 << self.variable_count)
            .map(|(index, &eval)| {
                let eq_point = HypercubePoint(index);
                EqFixedPoint::new(eq_point, self.variable_count, eval).eval(point)
            })
            .sum()
    }

    fn eval_hypercube(&self, point: HypercubePoint) -> F {
        self.evals[point.0]
    }

    fn fix_variable(&mut self, var: F) {
        debug_assert!(self.variable_count > 0);
        self.variable_count -= 1;
        let (first, second) = self.evals.split_at_mut(1 << self.variable_count);
        // This helps performance
        assert!(first.len() <= second.len());
        for (a, b) in first.iter_mut().zip(second.iter()) {
            // a = (1 - var) * a + var * b
            //   = a - a * var + var * b
            //   = a - var * (b - a)
            *a += var * (*b - *a);
        }
    }

    fn variable_count(&self) -> usize {
        self.variable_count
    }

    fn eval_field_then_hypercube(&self, f: F, hypercube_point: HypercubePoint) -> F {
        debug_assert!(hypercube_point.0 < (1 << (self.variable_count - 1)));
        let val_0 = self.evals[hypercube_point.0];
        let val_1 = self.evals[(1 << (self.variable_count - 1)) | hypercube_point.0];
        val_0 + f * (val_1 - val_0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::AdditiveGroup;
    use util::algebra::field::arkfield::Fp2256;

    #[test]
    fn test_hypercube_point_iter() {
        let hypercube_point = HypercubePoint(0b11011);
        assert_eq!(
            hypercube_point
                .into_point_iter::<Fp2256>(5)
                .collect::<Vec<_>>(),
            vec![
                Fp2256::ONE,
                Fp2256::ONE,
                Fp2256::ZERO,
                Fp2256::ONE,
                Fp2256::ONE
            ]
        );
        assert_eq!(
            hypercube_point.into_bool_iter(5).collect::<Vec<_>>(),
            vec![true, true, false, true, true]
        );
    }

    #[test]
    fn test_eq_fixed_point_var_then_hypercube() {
        let eq_fixed_point = EqFixedPoint::new(HypercubePoint(0b1011), 4, Fp2256::from(3));
        assert_eq!(
            eq_fixed_point.eval_hypercube(HypercubePoint(0b1011)),
            eq_fixed_point.eval_field_then_hypercube(Fp2256::ONE, HypercubePoint(0b011))
        );
        assert_eq!(
            eq_fixed_point.eval_hypercube(HypercubePoint(0b0011)),
            eq_fixed_point.eval_field_then_hypercube(Fp2256::ZERO, HypercubePoint(0b011))
        );
    }

    #[test]
    fn test_hypercube_eval_point_var_then_hypercube() {
        const VARS: usize = 4;
        let evals: Vec<_> = (0..(1 << VARS)).map(|i| Fp2256::from(i as u64)).collect();
        let eval_poly = HypercubeEvalPoly::new(&evals, VARS);
        assert_eq!(
            eval_poly.eval_hypercube(HypercubePoint(0b1011)),
            eval_poly.eval_field_then_hypercube(Fp2256::ONE, HypercubePoint(0b011))
        );
        assert_eq!(
            eval_poly.eval_hypercube(HypercubePoint(0b0011)),
            eval_poly.eval_field_then_hypercube(Fp2256::ZERO, HypercubePoint(0b011))
        );
    }

    #[test]
    fn test_hypercube_fix_variable() {
        let mut eval = HypercubeEvalPoly::new(
            [
                Fp2256::from(1),
                Fp2256::from(2),
                Fp2256::from(3),
                Fp2256::from(4),
            ].as_ref(),
            2,
        );
        let point = [Fp2256::from(462), Fp2256::from(34712)];
        let goal = eval.eval(&point);
        eval.fix_variable(point[0]);
        assert_eq!(eval.variable_count(), 1);
        eval.fix_variable(point[1]);
        assert_eq!(eval.variable_count(), 0);
        assert_eq!(eval.eval(&[]), goal);
    }

    #[test]
    fn test_fixed_eq() {
        let fixed_eq_point = EqFixedPoint::<Fp2256>::new(HypercubePoint(0b10), 2, Fp2256::ONE);
        assert_eq!(
            fixed_eq_point.eval_hypercube(HypercubePoint(0b00)),
            Fp2256::ZERO
        );
        assert_eq!(
            fixed_eq_point.eval_hypercube(HypercubePoint(0b01)),
            Fp2256::ZERO
        );
        assert_eq!(
            fixed_eq_point.eval_hypercube(HypercubePoint(0b10)),
            Fp2256::ONE
        );
        assert_eq!(
            fixed_eq_point.eval_hypercube(HypercubePoint(0b11)),
            Fp2256::ZERO
        );
        assert_eq!(
            fixed_eq_point.eval(&[Fp2256::ZERO, Fp2256::ZERO]),
            Fp2256::ZERO
        );
        assert_eq!(
            fixed_eq_point.eval(&[Fp2256::ZERO, Fp2256::ONE]),
            Fp2256::ZERO
        );
        assert_eq!(
            fixed_eq_point.eval(&[Fp2256::ONE, Fp2256::ZERO]),
            Fp2256::ONE
        );
        assert_eq!(
            fixed_eq_point.eval(&[Fp2256::ONE, Fp2256::ONE]),
            Fp2256::ZERO
        );

        assert_eq!(
            fixed_eq_point.eval(&[Fp2256::from(4), Fp2256::from(12)]),
            Fp2256::from(4) * (Fp2256::ONE - Fp2256::from(12))
        );
    }
}
