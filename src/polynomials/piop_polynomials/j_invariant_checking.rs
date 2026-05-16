//! The code to check a j-invariant.

use ark_ff::Field;

use crate::polynomials::{HypercubePoint};

#[derive(Eq, PartialEq, Clone)]
pub struct LinearPolynomial<F: Field>([F; 2]);

impl<F: Field> LinearPolynomial<F> {
    /// `a`: Evaluation at 0
    /// `b`: Evaluation at 1
    pub fn from_evals(a: F, b: F) -> Self {
        LinearPolynomial([a, b])
    }

    pub fn eval(&self, x: F) -> F {
        let [a, b] = self.0;
        if x == F::ZERO {
            a
        } else if x == F::ONE {
            b
        } else {
            a + x * (b - a)
        }
    }
}

/// This polynomial is just in a single variable, which will be fixed first
#[derive(Eq, PartialEq, Clone)]
pub struct JInvariantChecker<F: Field> {
    h: F,
    j: F,
    value: F,
    has_first_variable_been_fixed: bool,
    a: LinearPolynomial<F>,
    c: LinearPolynomial<F>,
    d1: LinearPolynomial<F>,
    d2: LinearPolynomial<F>,
    d3: LinearPolynomial<F>,
    d4: LinearPolynomial<F>,
    d5: LinearPolynomial<F>,
    x: Option<F>
}

impl<F: Field> JInvariantChecker<F> {
    pub fn new(
        h: F,
        j: F,
        a: LinearPolynomial<F>,
        c: LinearPolynomial<F>,
        d1: LinearPolynomial<F>,
        d2: LinearPolynomial<F>,
        d3: LinearPolynomial<F>,
        d4: LinearPolynomial<F>,
        d5: LinearPolynomial<F>,
    ) -> JInvariantChecker<F> {
        debug_assert_eq!(d1.eval(F::ZERO), a.eval(F::ZERO).pow([2]));
        debug_assert_eq!(d2.eval(F::ZERO), c.eval(F::ZERO).pow([2]));
        debug_assert_eq!(d3.0[0], d1.eval(F::ZERO).pow([2]));
        debug_assert_eq!(c.eval(F::ZERO) * d4.eval(F::ZERO), F::ONE);
        debug_assert_eq!(d5.eval(F::ZERO) * (F::from(4) * c.eval(F::ZERO) - d1.eval(F::ZERO)), F::ONE);
        JInvariantChecker {
            h,
            j,
            value: F::ZERO,
            has_first_variable_been_fixed: false,
            a,
            c,
            d1,
            d2,
            d3,
            d4,
            d5,
            x: None
        }
    }

    pub fn eval_at_x(&self, x: F) -> F {
        let a = self.a.eval(x);
        let c = self.c.eval(x);
        let d1 = self.d1.eval(x);
        let d2 = self.d2.eval(x);
        let d3 = self.d3.eval(x);
        let d4 = self.d4.eval(x);
        let d5 = self.d5.eval(x);
        let h = self.h;
        let j = self.j;
        let mut val = d5 * (F::from(3) * c - d1) - F::ONE;
        val = val * h + c * d4 - F::ONE;
        val = val * h + j * d2 * (F::from(4) * c - d1)
            + F::from(256) * (d1*d3 - F::from(3) * c*d2 - F::from(9) * c * d3 + F::from(27) * d2 * d1);
        val = val * h + d3 - d1.pow([2]);
        val = val * h + d2 - c.pow([2]);
        val = val * h + d1 - a.pow([2]);
        val *= h;
        val * x
    }

    pub fn eval_field_then_hypercube(&self, f: F, hypercube_point: HypercubePoint) -> F {
        if self.has_first_variable_been_fixed {
            self.value
        } else {
            // Compute polynomial
            self.eval_at_x(f)
        }
    }

    pub fn value(&self) -> F {
        assert!(self.has_first_variable_been_fixed);
        self.value
    }

    pub fn fix_variable(&mut self, value: F) {
        if !self.has_first_variable_been_fixed {
            let value = self.eval_at_x(value);
            self.value = value;
            self.has_first_variable_been_fixed = true;
            self.x = Some(value);
        }
    }
    pub fn get_openings(&self) -> [F; 7] {
        let x = self.x.unwrap();
        [
            self.a.eval(x),
            self.c.eval(x),
            self.d1.eval(x),
            self.d2.eval(x),
            self.d3.eval(x),
            self.d4.eval(x),
            self.d5.eval(x)
        ]
    }
    pub fn eval_with_opening(h: F, j: F, &[a, c, d1, d2, d3, d4, d5]: &[F; 7]) -> F {
        let mut val = d5 * (F::from(3) * c - d1) - F::ONE;
        val = val * h + c * d4 - F::ONE;
        val = val * h + j * d2 * (F::from(4) * c - d1)
            + F::from(256) * (d1*d3 - F::from(3) * c*d2 - F::from(9) * c * d3 + F::from(27) * d2 * d1);
        val = val * h + d3 - d1.pow([2]);
        val = val * h + d2 - c.pow([2]);
        val = val * h + d1 - a.pow([2]);
        val *= h;
        val
    }
}
