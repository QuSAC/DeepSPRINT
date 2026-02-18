use std::mem::size_of;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use util::{
    algebra::{field::FftField, polynomial::Polynomial},
    merkle_tree::{MERKLE_ROOT_SIZE, MerkleRoot},
    query_result::QueryResult,
};

pub mod prover;
pub mod verifier;

#[derive(Clone)]
pub struct DeepEval<T: FftField> {
    point: Vec<T>,
    first_eval: T,
    else_evals: Vec<T>,
}

impl<T: FftField> DeepEval<T> {
    pub fn new(point: Vec<T>, poly_hypercube: Vec<T>) -> Self {
        DeepEval {
            point: point.clone(),
            first_eval: Self::evaluatioin_at(point, poly_hypercube),
            else_evals: vec![],
        }
    }

    fn evaluatioin_at(point: Vec<T>, mut poly_hypercube: Vec<T>) -> T {
        let mut len = poly_hypercube.len();
        assert_eq!(len, 1 << point.len());
        for v in point.into_iter() {
            len >>= 1;
            for i in 0..len {
                poly_hypercube[i] *= T::ONE - v;
                let tmp = poly_hypercube[i + len] * v;
                poly_hypercube[i] += tmp;
            }
        }
        poly_hypercube[0]
    }

    pub fn append_else_eval(&mut self, poly_hypercube: Vec<T>) {
        let mut point = self.point[self.else_evals.len()..].to_vec();
        point[0] += T::ONE;
        self.else_evals
            .push(Self::evaluatioin_at(point, poly_hypercube));
    }

    pub fn verify(&self, challenges: &Vec<T>) -> T {
        let (_, challenges) = challenges.split_at(challenges.len() - self.point.len());
        let mut y_0 = self.first_eval;
        assert_eq!(self.point.len(), self.else_evals.len());
        for ((x, eval), challenge) in self
            .point
            .iter()
            .zip(self.else_evals.iter())
            .zip(challenges.into_iter())
        {
            let y_1 = eval.clone();
            y_0 += (y_1 - y_0) * (challenge.clone() - x.clone());
        }
        y_0
    }
}

pub struct Commit<T: FftField> {
    merkle_root: MerkleRoot,
    deep: T,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SerializableCommit<T: FftField + CanonicalSerialize + CanonicalDeserialize> {
    merkle_root: MerkleRoot,
    deep: T,
}

impl<T: FftField + CanonicalSerialize + CanonicalDeserialize> From<Commit<T>>
    for SerializableCommit<T>
{
    fn from(value: Commit<T>) -> Self {
        SerializableCommit {
            merkle_root: value.merkle_root,
            deep: value.deep,
        }
    }
}
impl<T: FftField + CanonicalSerialize + CanonicalDeserialize> From<SerializableCommit<T>>
    for Commit<T>
{
    fn from(value: SerializableCommit<T>) -> Self {
        Commit {
            merkle_root: value.merkle_root,
            deep: value.deep,
        }
    }
}

#[derive(Clone)]
pub struct Proof<T: FftField> {
    merkle_root: Vec<MerkleRoot>,
    query_result: Vec<QueryResult<T>>,
    deep_evals: Vec<(T, Vec<T>)>,
    shuffle_evals: Vec<T>,
    evaluation: T,
    final_value: T,
    final_poly: Polynomial<T>,
}

impl<T: FftField> Proof<T> {
    #[deprecated(note = "use `CanonicalSerialize`")]
    #[allow(deprecated)]
    pub fn size(&self) -> usize {
        self.merkle_root.len() * MERKLE_ROOT_SIZE
            + self
                .query_result
                .iter()
                .fold(0, |acc, x| acc + x.proof_size())
            + (self.deep_evals.iter().fold(0, |acc, x| acc + x.1.len())
                + self.shuffle_evals.len()
                + 2)
                * size_of::<T>()
    }

    pub fn evaluation(&self) -> T {
        self.evaluation
    }
}

// Does not include `final_poly`, as with the proof size
impl<T> CanonicalSerialize for Proof<T>
where
    T: FftField + CanonicalSerialize,
{
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        for root in &self.merkle_root {
            root.serialize_with_mode(&mut writer, compress)?;
        }
        for q in &self.query_result {
            q.serialize_with_mode(&mut writer, compress)?;
        }
        for (_, deep_eval) in &self.deep_evals {
            for eval in deep_eval {
                eval.serialize_with_mode(&mut writer, compress)?;
            }
        }
        for s in &self.shuffle_evals {
            s.serialize_with_mode(&mut writer, compress)?;
        }
        self.evaluation.serialize_with_mode(&mut writer, compress)?;
        self.final_value
            .serialize_with_mode(&mut writer, compress)?;
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.merkle_root
            .iter()
            .map(|root| root.serialized_size(compress))
            .sum::<usize>()
            + self
                .query_result
                .iter()
                .map(|q| q.serialized_size(compress))
                .sum::<usize>()
            + self
                .deep_evals
                .iter()
                .flat_map(|d| d.1.iter().map(|a| a.serialized_size(compress)))
                .sum::<usize>()
            + self
                .shuffle_evals
                .iter()
                .map(|s| s.serialized_size(compress))
                .sum::<usize>()
            + self.evaluation.serialized_size(compress)
            + self.final_value.serialized_size(compress)
    }
}

#[cfg(test)]
mod tests {
    use crate::{prover::Prover, verifier::Verifier};
    use ark_ff::FftField;
    use ark_serialize::CanonicalSerialize;
    use csv::Writer;
    use util::algebra::field::FftField;
    use util::{CODE_RATE, SECURITY_BITS, SIZE, STEP};
    use util::{
        algebra::{coset::Coset, field::arkfield::Fp2256, polynomial::MultilinearPolynomial},
        random_oracle::RandomOracle,
    };

    #[allow(deprecated)]
    fn output_proof_size<T: FftField + CanonicalSerialize>(variable_num: usize) -> usize {
        let polynomial = MultilinearPolynomial::random_polynomial(variable_num);
        let mut interpolate_cosets =
            vec![Coset::new(1 << (variable_num + CODE_RATE), T::from_int(1))];
        for i in 1..variable_num + 1 {
            interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
        }
        let oracle = RandomOracle::new(variable_num, SECURITY_BITS / CODE_RATE);
        let prover = Prover::new(variable_num, &interpolate_cosets, polynomial, &oracle, STEP);
        let commit = prover.commit_polynomial();
        let verifier = Verifier::new(variable_num, &interpolate_cosets, commit, &oracle, STEP);
        let point = verifier.get_open_point();
        let proof = prover.generate_proof(point);
        let size = proof.compressed_size();
        assert_eq!(size, proof.size());
        assert!(verifier.verify(proof));
        size
    }

    #[test]
    fn test_proof_size() {
        let mut wtr = Writer::from_path("deepfold.csv").unwrap();
        let range = 10..SIZE;
        for i in range.clone() {
            let proof_size = output_proof_size::<Fp2256>(i);
            wtr.write_record(&[i.to_string(), proof_size.to_string()])
                .unwrap();
        }
    }
}
