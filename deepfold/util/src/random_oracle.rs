use ark_ff::Field;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct RandomOracle<T: Field> {
    pub beta: T,
    pub rlc: T,
    pub folding_challenges: Vec<T>,
    pub deep: Vec<T>,
    pub alpha: Vec<T>,
    pub query_list: Vec<usize>,
}

impl<T: Field> RandomOracle<T> {
    pub fn new(rng: &mut impl Rng, total_round: usize, query_num: usize) -> Self {
        RandomOracle {
            beta: T::rand(rng),
            rlc: T::rand(rng),
            folding_challenges: (0..total_round)
                .into_iter()
                .map(|_| T::rand(rng))
                .collect(),
            deep: (0..total_round)
                .into_iter()
                .map(|_| T::rand(rng))
                .collect(),
            alpha: (0..total_round)
                .into_iter()
                .map(|_| T::rand(rng))
                .collect(),
            query_list: (0..query_num)
                .into_iter()
                .map(|_| rand::thread_rng().r#gen())
                .collect(),
        }
    }
}
