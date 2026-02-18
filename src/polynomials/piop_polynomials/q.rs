use ark_ff::Field;

use crate::polynomials::{EqFixedPoint, HypercubeEvalPoly, HypercubePoint, PiopPolynomial};

/// The polynomial Q, which is defined as
/// `Q(i,j,k)=sum_l eq(i,l)*eq(j,l+1)*eq(k,l)`.
///
/// Since `k` is always supplied in advance, we employ two strategies, one for
/// `i` and one for `j`. At first, we fully evaluate `eq(k,l)`. If `j` is a
/// hypercube point, we know exactly which single entry will get a component.
/// `i` will then be fully determined by `j`. Once we start fixing values of
/// `i`, we enter it in our precomputed values. Then, as soon as we're ready to
/// evaluate for `j`, we evaluate this precomputed polynomial as usual.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Q<const EVALS: usize, F: Field> {
    polyomial_over_j: HypercubeEvalPoly<F>,
    variable_count_i: usize,
}

impl<const EVALS: usize, F: Field> Q<EVALS, F> {
    /// Construct Q using gray codes. This should be more efficient since it
    /// uses fewer multiplications. It requires that none of the entries of k
    /// are 0 or 1, which should be satisfied for randomly chosen k.
    pub fn new_gray_codes(k: &[F]) -> Self {
        // This does not work for values on the boolean hypercube
        debug_assert!(!k.iter().any(|f| f.is_one() || f.is_zero()));
        assert_eq!(1 << k.len(), EVALS);

        let mut evals = vec![F::ZERO; EVALS];
        // The value for the evaluation. At the start, the index will simply be
        // 0. At this point, we evaluate as prod_v (1 - k_v)
        let mut current_value = k.iter().map(|kv| F::ONE - kv).product();
        for i in 0..EVALS {
            // We access as a grey code, to minimize the number of
            // multiplications that we do.
            let gray_code = i ^ (i >> 1);

            // Set the evaluation at this index. Next, we compute the value for
            // the next index. We shift by 1 to account for the difference in
            // indexing between `j` and `k`: the latter is one higher. (So 0
            // has value 0.)
            if gray_code + 1 < EVALS {
                evals[gray_code + 1] = current_value;
            }

            if i + 1 < EVALS {
                // The entry that changes in the next iteration is always the least
                // significant zero entry in the index.
                let entry = i.trailing_ones();
                // Check if the entry is currently on.
                let currently_on = ((1 << entry) & gray_code) != 0;
                if currently_on {
                    // Switch off if it is currently on.
                    current_value *= (F::ONE - k[k.len() - 1 - entry as usize])
                        / k[k.len() - 1 - entry as usize];
                } else {
                    // Switch on if it is currently off.
                    current_value *= k[k.len() - 1 - entry as usize]
                        / (F::ONE - k[k.len() - 1 - entry as usize]);
                }
            }
        }

        Self {
            variable_count_i: k.len(),
            polyomial_over_j: HypercubeEvalPoly::new(&evals, k.len()),
        }
    }

    pub fn new(k: &[F]) -> Self {
        assert_eq!(1 << k.len(), EVALS);
        let evals: Vec<_> = (0..EVALS)
            .map(|i| {
                if i == 0 {
                    F::ZERO
                } else {
                    HypercubePoint(i - 1)
                        .into_bool_iter(k.len())
                        .zip(k.iter())
                        .map(|(b, v)| if b { *v } else { F::ONE - v })
                        .product()
                }
            })
            .collect();
        Self {
            variable_count_i: k.len(),
            polyomial_over_j: HypercubeEvalPoly::new(&evals, k.len()),
        }
    }

    /// Compute Q directly in case all variables are fixed
    pub fn direct_eval(i: &[F], j: &[F], k: &[F]) -> F {
        (0..(1 << i.len()) - 1)
            .map(HypercubePoint)
            .map(|l| {
                EqFixedPoint::new(l, i.len(), F::ONE).eval(i)
                    * EqFixedPoint::new(HypercubePoint(l.0 + 1), j.len(), F::ONE).eval(j)
                    * EqFixedPoint::new(l, i.len(), F::ONE).eval(k)
            })
            .sum()
    }
}

impl<const EVALS: usize, F: Field> PiopPolynomial<F> for Q<EVALS, F> {
    fn variable_count(&self) -> usize {
        self.variable_count_i + self.polyomial_over_j.variable_count()
    }

    fn eval(&self, point: &[F]) -> F {
        self.polyomial_over_j
            .evals()
            .iter()
            .enumerate()
            .skip(1)
            .map(|(j, val)| {
                // Compute the value of l that yields this value of eq(l + 1, j) * eq(l, k)
                let l = j - 1;
                // Compute the value of ij that yields this value
                let free_variables_mask = (1 << self.variable_count()) - 1;
                let ij = HypercubePoint(
                    (l << self.polyomial_over_j.variable_count() | j) & free_variables_mask,
                );
                // Compute eq(l, i) * eq(l+1, j)
                ij.into_bool_iter(self.variable_count())
                    .zip(point.iter())
                    .map(|(b, p)| if b { *p } else { F::ONE - p })
                    .product::<F>()
                    * val
            })
            .sum()
    }

    fn eval_hypercube(&self, point: crate::polynomials::HypercubePoint) -> F {
        if self.variable_count_i == self.polyomial_over_j.variable_count() {
            let (i, j) = point.split(self.variable_count_i, self.variable_count_i);
            if i.0 + 1 == j.0 || self.variable_count() == 0 {
                self.polyomial_over_j.eval_hypercube(j)
            } else {
                F::ZERO
            }
        } else {
            // TODO: This is not actually used in the implementation, but it could be more efficient
            let p = point.into_bool_iter(self.variable_count()).next().unwrap();
            let new_point = point.take_last_n(self.variable_count() - 1);
            if p {
                self.eval_field_then_hypercube(F::ONE, new_point)
            } else {
                self.eval_field_then_hypercube(F::ZERO, new_point)
            }
        }
    }

    fn eval_field_then_hypercube(
        &self,
        f: F,
        hypercube_point: crate::polynomials::HypercubePoint,
    ) -> F {
        if self.variable_count_i > 0 {
            // Only one term will be non-zero, that corresponding to l+1 = j. We
            // can therefore compute l from j.

            // The number of variables remaining for i. The first will be a
            // field element.
            let hypercube_variables_i = self.variable_count_i - 1;
            let (i_remainder, j) = hypercube_point.split(
                hypercube_variables_i,
                self.polyomial_over_j.variable_count(),
            );
            // If j = 0, this will fail, but then Q should be zero.
            if let Some(l) = j.0.checked_sub(1).map(HypercubePoint) {
                // We evaluate the part of the hypercube, if i=l this should be
                // 1, and 0 otherwise.
                if i_remainder.take_last_n(hypercube_variables_i)
                    == l.take_last_n(hypercube_variables_i)
                {
                    // The part for eq(l, k)*eq(l+1, j) is straightforward as
                    // well, since this was precomputed already.
                    let eq_k_l_eq_lplus1_j = self.polyomial_over_j.eval_hypercube(j);
                    // Finally, compute eq(l_v, i_v), wjhere v is the index of
                    // the free variable.
                    let iv_eq_lv = if l.value_at_pos(self.variable_count_i, 0) {
                        f
                    } else {
                        F::ONE - f
                    };
                    iv_eq_lv * eq_k_l_eq_lplus1_j
                } else {
                    F::ZERO
                }
            } else {
                F::ZERO
            }
        } else {
            self.polyomial_over_j
                .eval_field_then_hypercube(f, hypercube_point)
        }
    }

    fn fix_variable(&mut self, var: F) {
        if self.variable_count_i > 0 {
            // Fix a variable i, by fixing all evaluations
            for (j, val) in self
                .polyomial_over_j
                .evals_mut()
                .iter_mut()
                .enumerate()
                .skip(1)
            {
                let i = HypercubePoint(j - 1);
                if i.value_at_pos(self.variable_count_i, 0) {
                    *val *= var;
                } else {
                    *val -= var * *val;
                }
            }
            self.variable_count_i -= 1;
        } else {
            self.polyomial_over_j.fix_variable(var);
        }
    }
}

#[cfg(test)]
mod test {
    use crate::polynomials::{piop_polynomials::Q, HypercubePoint, PiopPolynomial};
    use ark_ff::Zero;
    use ark_ff::{AdditiveGroup, Field};
    use std::array;
    use util::algebra::field::arkfield::Fp2256;

    #[test]
    fn naive_method_matches_gray_coding() {
        const EVALS: usize = 1 << 8;
        let k: [Fp2256; 8] = array::from_fn(|i| Fp2256::from(i as u64 + 2));
        assert_eq!(
            Q::<EVALS, Fp2256>::new_gray_codes(&k),
            Q::<EVALS, _>::new(&k)
        );
    }

    #[test]
    fn eval_hypercube_matches_field_then_hypercube() {
        const EVALS: usize = 1 << 8;
        let k: [Fp2256; 8] = array::from_fn(|i| Fp2256::from(i as u64 + 2));
        let q = Q::<EVALS, Fp2256>::new_gray_codes(&k);
        assert_eq!(
            q.eval_hypercube(HypercubePoint(0b1001_1011_1011_0111)),
            q.eval_field_then_hypercube(Fp2256::ONE, HypercubePoint(0b001_1011_1011_0111))
        );
        assert_eq!(
            q.eval_hypercube(HypercubePoint(0b0101_0101_0101_0110)),
            q.eval_field_then_hypercube(Fp2256::ZERO, HypercubePoint(0b101_0101_0101_0110))
        );
    }

    #[test]
    fn test_q_zero_for_zero_j() {
        let k: [Fp2256; 8] = array::from_fn(|i| Fp2256::from(i as u64 + 2));
        let q = Q::<{ 1 << 8 }, _>::new_gray_codes(&k);
        for i in 0..(1 << 8) {
            let i = HypercubePoint(i);
            let j = HypercubePoint(0);
            let ij = HypercubePoint(i.0 << 8 | j.0);
            assert_eq!(q.eval_hypercube(ij), Fp2256::ZERO);
        }
    }

    #[test]
    fn test_zeros() {
        let k = [Fp2256::from(2); 8];
        let q = Q::<{ 1 << 8 }, _>::new(&k);
        for j in 0..(1 << 8) {
            for i in 0..(1 << 8) {
                if j != i + 1 {
                    let ij = HypercubePoint((i << 8) | j);
                    assert!(q.eval_hypercube(ij).is_zero());
                }
            }
        }
    }

    #[test]
    fn test_fix_variable() {
        let k: [Fp2256; 8] = array::from_fn(|i| Fp2256::from(10 * i as u64 + 24));
        let mut q = Q::<{ 1 << 8 }, _>::new_gray_codes(&k);

        let f = Fp2256::from(2342);
        let value = q.eval_field_then_hypercube(f, HypercubePoint(0b001_1011_1011_0111));
        q.fix_variable(f);
        let value_2 =
            q.eval_field_then_hypercube(Fp2256::ZERO, HypercubePoint(0b01_1011_1011_0111));
        assert_eq!(value, value_2);
        let value_3 = q.eval_field_then_hypercube(Fp2256::ZERO, HypercubePoint(0b1_1011_1011_0111));
        assert_eq!(value, value_3);
        let value_4 = q.eval_field_then_hypercube(Fp2256::ONE, HypercubePoint(0b1011_1011_0111));
        assert_eq!(value, value_4);
        let value_5 = q.eval_field_then_hypercube(Fp2256::ONE, HypercubePoint(0b011_1011_0111));
        assert_eq!(value, value_5);
        let value_6 = q.eval_field_then_hypercube(Fp2256::ZERO, HypercubePoint(0b11_1011_0111));
        assert_eq!(value, value_6);
        let value_7 = q.eval_field_then_hypercube(Fp2256::ONE, HypercubePoint(0b1_1011_0111));
        assert_eq!(value, value_7);
        let value_8 = q.eval_field_then_hypercube(Fp2256::ONE, HypercubePoint(0b1011_0111));
        assert_eq!(value, value_8);
    }

    #[test]
    fn test_fix_variable_everywhere() {
        const LOG_2_PATH_LENGTH: usize = 2;
        const EVALS: usize = 1 << LOG_2_PATH_LENGTH;
        let k: [Fp2256; LOG_2_PATH_LENGTH] = array::from_fn(|i| Fp2256::from(2 << i));
        let q = Q::<EVALS, _>::new_gray_codes(&k);
        let f = Fp2256::from(2);
        for eval_point in 0..(1 << (q.variable_count() - 1)) {
            let mut q = q.clone();
            let value = q.eval_field_then_hypercube(f, HypercubePoint(eval_point));
            q.fix_variable(f);
            let value_2 = q.eval_hypercube(HypercubePoint(eval_point));
            assert_eq!(value, value_2);
        }
    }

    #[test]
    fn test_fix_variable_everywhere_last() {
        const LOG_2_PATH_LENGTH: usize = 2;
        const EVALS: usize = 1 << LOG_2_PATH_LENGTH;
        let k: [Fp2256; LOG_2_PATH_LENGTH] = array::from_fn(|i| Fp2256::from(2 << i));
        // Construct q
        let mut q = Q::<EVALS, _>::new_gray_codes(&k);
        // Fix variables until we have two variables left
        while q.variable_count() > 2 {
            q.fix_variable(Fp2256::from(146));
        }

        let f = Fp2256::from(2);
        for eval_point in 0..(1 << (q.variable_count() - 1)) {
            let mut q = q.clone();
            let value = q.eval_field_then_hypercube(f, HypercubePoint(eval_point));
            q.fix_variable(f);
            let value_2 = q.eval_hypercube(HypercubePoint(eval_point));
            assert_eq!(value, value_2);
        }
    }

    #[test]
    fn test_direct_eval() {
        const LOG_2_PATH_LENGTH: usize = 2;
        const EVALS: usize = 1 << LOG_2_PATH_LENGTH;
        let k: [Fp2256; LOG_2_PATH_LENGTH] = array::from_fn(|i| Fp2256::from(i as u64 + 1));
        let i: [Fp2256; LOG_2_PATH_LENGTH] = array::from_fn(|i| Fp2256::from(i as u64 + 20));
        let j: [Fp2256; LOG_2_PATH_LENGTH] = array::from_fn(|i| Fp2256::from(i as u64 + 30));
        let mut q = Q::<EVALS, _>::new(&k);
        let direct_eval_value = Q::<EVALS, _>::direct_eval(&i, &j, &k);
        for val in i.iter().chain(j.iter()) {
            q.fix_variable(*val);
        }
        assert_eq!(q.eval_hypercube(HypercubePoint(0)), direct_eval_value);
    }
}
