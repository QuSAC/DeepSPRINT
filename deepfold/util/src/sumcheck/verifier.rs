use ark_ff::Field;

use crate::{random_oracle::RandomOracle};

pub struct SumcheckVerifier<T: Field> {
    total_round: usize,
    oracle: RandomOracle<T>,
    sumcheck_values: Vec<(T, T, T)>,
}

impl<T: Field> SumcheckVerifier<T> {
    pub fn new(total_round: usize, oracle: &RandomOracle<T>) -> Self {
        SumcheckVerifier {
            total_round,
            oracle: oracle.clone(),
            sumcheck_values: vec![],
        }
    }

    pub fn receive_sumcheck_value(&mut self, value: (T, T, T)) {
        self.sumcheck_values.push(value);
    }

    // The final output of sumcheck is the random point and a random evaluation of polynomial
    // Return them
    pub fn verify(&self) -> (Vec<T>, T) {
        let mut sum = self.sumcheck_values[0].0 + self.sumcheck_values[0].1;
        for k in 0..self.total_round {
            let challenge = self.oracle.folding_challenges[k];
            sum = self.process_sumcheck(challenge, sum, self.sumcheck_values[k]);
        }
        (
            self.oracle.folding_challenges[0..self.total_round].to_vec(),
            sum,
        )
    }

    fn process_sumcheck(&self, challenge: T, last_sum: T, sumcheck_values: (T, T, T)) -> T {
        let x_0 = sumcheck_values.0;
        let x_1 = sumcheck_values.1;
        let x_2 = sumcheck_values.2;
        let two_inverse = T::from(2).inverse().unwrap();
        assert_eq!(last_sum, x_0 + x_1);
        let sum =
            x_0 * (T::from(1) - challenge) * (T::from(2) - challenge) * two_inverse
                + x_1 * challenge * (T::from(2) - challenge)
                + x_2 * challenge * (challenge - T::from(1)) * two_inverse;
        sum
    }
}
