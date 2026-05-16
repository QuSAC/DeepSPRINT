use ark_ff::Field;

use crate::polynomials::{HypercubePoint, PiopPolynomial};

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

impl<'a, const VARIABLE_COUNT: usize, F: Field> Mask<'a, VARIABLE_COUNT, F> {
    /// Create a new mask with all the evaluations of the polynomial.
    pub fn new_from_evals(b_0: &'a [F], b_1: &'a [F]) -> Self {
        let b_0_mask = &b_0[1..2 + VARIABLE_COUNT * 3];
        let b_1_mask = &b_1[1..2 + VARIABLE_COUNT * 3];
        Self::new(b_0_mask, b_1_mask)
    }

    /// Create a new mask using just the indices that are accessed.
    fn new(b_0_mask: &'a [F], b_1_mask: &'a [F]) -> Self {
        // The mask should have 3 coefficients per variable, plus a constant
        // term. This number is then doubled to mask the coefficient opening.
        debug_assert_eq!(b_0_mask.len(), 3 * VARIABLE_COUNT + 1);
        debug_assert_eq!(b_1_mask.len(), 3 * VARIABLE_COUNT + 1);
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

    pub fn eval_field_then_sum_hypercube(&self, point: F) -> F {
        // Evaluate for the current variable. The first index contains the constant term.s
        let starting_index_first = 1 + self.current_variable_index() * 3;
        let starting_index_next = starting_index_first + 3;
        let remaining_variable_count = self.variable_count - 1;
        let [c, b, a] = self.b_0_mask[starting_index_first..starting_index_next]
            .try_into()
            .unwrap();
        ((self.current_total + ((a * point + b) * point + c) * point) * F::from(2)
            // Add the remaining variables. Precisely half the points on the
            // hypercube will give a contribution to any sum.
            + self.b_0_mask[starting_index_next..].iter()
            .sum::<F>())
            * F::from((1 << remaining_variable_count) / 2)
    }

    /// Compute the constant term that would fix the coefficients.
    pub fn fix_constant_term(b_0_evals: &'a mut [F], b_1_evals: &'a mut [F]) {
        let n = (VARIABLE_COUNT - 1) / 2;
        let path_length = 2usize.pow(n as u32);
        debug_assert_eq!(b_0_evals.len(), 2 * path_length);
        debug_assert_eq!(b_1_evals.len(), 2 * path_length);

        let coefficients = &mut b_0_evals[1..2 + 3 * VARIABLE_COUNT];
        Self::fix_constant_term_from_mask_slice(coefficients);
    }

    fn fix_constant_term_from_mask_slice(coefficients: &mut [F]) {
        debug_assert_eq!(coefficients.len(), 1 + 3 * VARIABLE_COUNT);
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

impl<'a, const VARIABLE_COUNT: usize, F: Field> PiopPolynomial<F> for Mask<'a, VARIABLE_COUNT, F> {
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

#[cfg(test)]
mod tests {
    use std::array;

    use ark_ff::{AdditiveGroup, Field};
    use util::algebra::field::arkfield::Fp2256;

    use crate::polynomials::{piop_polynomials::mask::Mask, HypercubePoint, PiopPolynomial};

    #[test]
    fn test_mask_fix() {
        const VARIABLES: usize = 2 * 8 + 1;
        let mut mask_evals = [Fp2256::ZERO; { 1 + 3 * 17 }];
        for (a, i) in mask_evals.iter_mut().enumerate() {
            *i = Fp2256::from(a as u64);
        }
        Mask::<VARIABLES, _>::fix_constant_term_from_mask_slice(&mut mask_evals);
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
    #[test]
    fn test_mask_eval_hypercube() {
        const ONE: Fp2256 = Fp2256::ONE;
        const MASK_VARIABLE_COUNT: usize = 3 * 17 + 1;
        let mut val = array::from_fn::<_, { MASK_VARIABLE_COUNT }, _>(|i| Fp2256::from(i as u32));
        Mask::<17, _>::fix_constant_term_from_mask_slice(&mut val);
        let mask: Mask<17, _> = Mask::new(&val, &val);
        assert_eq!(
            mask.eval_hypercube(HypercubePoint(0b1_0110_1101_0010_1010)),
            mask.eval_field_then_hypercube(ONE, HypercubePoint(0b0110_1101_0010_1010))
        );
    }

    #[test]
    fn test_mask_new() {
        const MASK_COEFFICIENT_COUNT: usize = 3 * 17 + 1;
        let mut val: [Fp2256; MASK_COEFFICIENT_COUNT] =
            [Fp2256::from(32423); MASK_COEFFICIENT_COUNT];
        Mask::<17, _>::fix_constant_term_from_mask_slice(&mut val);
        Mask::<17, _>::new(&val, &val);
    }

    #[test]
    fn test_mask_eval_hypercube_field() {
        const ONE: Fp2256 = Fp2256::ONE;
        const ZERO: Fp2256 = Fp2256::ZERO;
        const MASK_COEFFICIENT_COUNT: usize = 3 * 17 + 1;
        let mut val: [Fp2256; MASK_COEFFICIENT_COUNT] = array::from_fn(|i| Fp2256::from(i as u32));
        Mask::<17, _>::fix_constant_term_from_mask_slice(&mut val);
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
        Mask::<17, _>::fix_constant_term_from_mask_slice(&mut val);
        let mask: Mask<17, _> = Mask::new(&val, &val);
        let assignment = Fp2256::from(2352);
        let sum = (0..(1 << 16))
            .map(HypercubePoint)
            .map(|point| mask.eval_field_then_hypercube(assignment, point))
            .sum::<Fp2256>();
        assert_eq!(mask.eval_field_then_sum_hypercube(assignment), sum);
    }
}
