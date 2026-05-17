//! The code to check a j-invariant.

use ark_ff::Field;

use crate::polynomials::HypercubePoint;

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
    num_variables: usize,
    original_num_variables: usize,
    a: LinearPolynomial<F>,
    c: LinearPolynomial<F>,
    d1: LinearPolynomial<F>,
    d2: LinearPolynomial<F>,
    d3: LinearPolynomial<F>,
    d4: LinearPolynomial<F>,
    d5: LinearPolynomial<F>,
    x: Option<F>,
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
        variable_count: usize
    ) -> JInvariantChecker<F> {
        // A^2
        debug_assert_eq!(d1.eval(F::ONE), a.eval(F::ONE).pow([2]));
        // C^2
        debug_assert_eq!(d2.eval(F::ONE), c.eval(F::ONE).pow([2]));
        // A^4
        debug_assert_eq!(d3.eval(F::ONE), d1.eval(F::ONE).pow([2]));
        // 1 / C
        debug_assert_eq!(c.eval(F::ONE) * d4.eval(F::ONE), F::ONE);
        debug_assert_eq!(
            d5.eval(F::ONE) * (F::from(4) * c.eval(F::ONE) - d1.eval(F::ONE)),
            F::ONE
        );
        debug_assert_eq!(
            j * d2.eval(F::ONE) * (F::from(4) * c.eval(F::ONE) - d1.eval(F::ONE))
                - F::from(256)
                    * (-d1.eval(F::ONE) * d3.eval(F::ONE)
                        + F::from(3 * 3) * c.eval(F::ONE) * d3.eval(F::ONE)
                        - F::from(9 * 3) * d2.eval(F::ONE) * d1.eval(F::ONE)
                        + F::from(27) * d2.eval(F::ONE) * c.eval(F::ONE)),
            F::ZERO
        );
        let checker = JInvariantChecker {
            h,
            j,
            value: F::ZERO,
            num_variables: variable_count,
            original_num_variables: variable_count,
            a,
            c,
            d1,
            d2,
            d3,
            d4,
            d5,
            x: None,
        };
        debug_assert_eq!(
            checker.eval_at_x(F::ZERO) + checker.eval_at_x(F::ONE),
            F::ZERO
        );
        return checker;
    }

    pub fn num_variables(&self) -> usize {
        self.num_variables
    }

    pub fn eval_then_sum_hypercube(&self, val: F) -> F {
        assert!(self.num_variables > 0);
        if !self.has_first_variable_been_fixed() {
            self.eval_at_x(val) * F::from(2u32.pow((self.num_variables - 1) as u32))
        } else if self.num_variables == 1 {
            self.value()
        } else {
            self.value() * F::from(2u32.pow((self.num_variables - 1) as u32))
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
        Self::eval_with_opening(x, h, j, &[a, c, d1, d2, d3, d4, d5])
    }

    fn has_first_variable_been_fixed(&self) -> bool {
        self.original_num_variables != self.num_variables
    }

    pub fn eval_field_then_hypercube(&self, f: F, _hypercube_point: HypercubePoint) -> F {
        if self.has_first_variable_been_fixed() {
            self.value
        } else {
            // Compute polynomial
            self.eval_at_x(f)
        }
    }

    pub fn value(&self) -> F {
        assert!(self.has_first_variable_been_fixed());
        self.value
    }

    pub fn fix_variable(&mut self, var: F) {
        assert!(self.num_variables > 0);
        if !self.has_first_variable_been_fixed() {
            let value = self.eval_at_x(var);
            self.value = value;
            self.x = Some(var);
        }
        self.num_variables -= 1;
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
            self.d5.eval(x),
        ]
    }
    pub fn eval_with_opening(x: F, h: F, j: F, &[a, c, d1, d2, d3, d4, d5]: &[F; 7]) -> F {
        let mut val = d5 * (F::from(4) * c - d1) - F::ONE;
        val = val * h + c * d4 - F::ONE;
        val = val * h + j * d2 * (F::from(4) * c - d1)
            - F::from(256)
                * (-d1 * d3 + F::from(3 * 3) * c * d3 - F::from(9 * 3) * d2 * d1
                    + F::from(27) * d2 * c);
        val = val * h + d3 - d1.pow([2]);
        val = val * h + d2 - c.pow([2]);
        val = val * h + d1 - a.pow([2]);
        val *= h;
        val *= x;
        val
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{AdditiveGroup, Field};
use util::algebra::field::sqisign::level_i::Fp2251;

use crate::{keys::compute_j_invariant, polynomials::{HypercubePoint, piop_polynomials::j_invariant_checking::{JInvariantChecker, LinearPolynomial}}};
    #[test]
    fn test_fix_variable() {
        let h = Fp2251::from(14);
        let a = Fp2251::from(10);
        let c = Fp2251::from(5);
        let j = compute_j_invariant(a, c);
        let d1 = a * a;
        let d2 = c * c;
        let d3 = d1 * d1;
        let d4 = c.inverse().unwrap();
        let d5 = (Fp2251::from(4) * c - d1).inverse().unwrap();

        let a = LinearPolynomial::from_evals(141214.into(), a);
        let c = LinearPolynomial::from_evals(35345.into(), c);
        let d1 = LinearPolynomial::from_evals(12749.into(), d1);
        let d2 = LinearPolynomial::from_evals(1124.into(), d2);
        let d3 = LinearPolynomial::from_evals(6459.into(), d3);
        let d4 = LinearPolynomial::from_evals(456.into(), d4);
        let d5 = LinearPolynomial::from_evals(36.into(), d5);

        let mut tester = JInvariantChecker::new(h, j, a, c, d1, d2, d3, d4, d5, 17);
        assert_eq!(tester.eval_field_then_hypercube(0.into(), HypercubePoint(0)), 0.into());
        assert_eq!(tester.eval_field_then_hypercube(1.into(), HypercubePoint(0)), 0.into());

        let x = 124.into();
        let value = tester.eval_at_x(x);
        let mut sum = tester.eval_then_sum_hypercube(x);
        tester.fix_variable(x);
        sum /= Fp2251::from(2);

        for variables in (1..=16).rev() {
            assert_eq!(tester.num_variables(), variables);
            assert_eq!(tester.value(), value);
            assert_eq!(tester.eval_then_sum_hypercube(Fp2251::from(variables as u64)), sum);

            tester.fix_variable(Fp2251::from(variables as u64));
            if tester.num_variables() > 0 {
                assert_eq!(tester.eval_then_sum_hypercube(Fp2251::ZERO) + tester.eval_then_sum_hypercube(Fp2251::ONE), sum);
            }
            sum /= Fp2251::from(2);
        }
    }
}
