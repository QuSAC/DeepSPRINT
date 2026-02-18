use std::mem::size_of;

use ark_serialize::CanonicalSerialize;
use util::{
    algebra::{field::MyField, polynomial::MultilinearPolynomial},
    merkle_tree::{MerkleRoot, MERKLE_ROOT_SIZE},
    query_result::QueryResult,
};

pub mod prover;
pub mod verifier;

#[derive(Clone)]
pub struct DeepEval<T: MyField> {
    point: Vec<T>,
    first_eval: T,
    else_evals: Vec<T>,
}

impl<T: MyField> DeepEval<T> {
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
                poly_hypercube[i] *= T::from_int(1) - v;
                let tmp = poly_hypercube[i + len] * v;
                poly_hypercube[i] += tmp;
            }
        }
        poly_hypercube[0]
    }

    pub fn append_else_eval(&mut self, poly: &MultilinearPolynomial<T>) {
        let mut point = self.point[self.else_evals.len()..].to_vec();
        point[0] += T::from_int(1);
        self.else_evals.push(poly.evaluate(&point));
    }

    pub fn verify(&self, challenges: &Vec<T>, out_evals: &Vec<T>) -> T {
        let (_, challenges) = challenges.split_at(challenges.len() - self.point.len());
        let mut y_0 = self.first_eval;
        assert_eq!(self.point.len(), self.else_evals.len());
        assert_eq!(challenges.len(), out_evals.len() + 1);
        for ((x, eval), (i, challenge)) in self
            .point
            .iter()
            .zip(self.else_evals.iter())
            .zip(challenges.into_iter().enumerate())
        {
            let y_1 = eval.clone();
            y_0 += (y_1 - y_0) * (challenge.clone() - x.clone());
            if i.clone() < out_evals.len() {
                y_0 += out_evals[i] * challenge.clone() * challenge.clone();
            }
        }
        y_0
    }
}

pub struct Commit<T: MyField> {
    merkle_root: MerkleRoot,
    deep: T,
}

#[derive(Clone)]
pub struct Proof<T: MyField> {
    merkle_root: Vec<MerkleRoot>,
    query_result: (Vec<QueryResult<T>>, Vec<QueryResult<T>>),
    deep_evals: Vec<(T, Vec<T>)>,
    shuffle_evals: Vec<T>,
    out_evals: Vec<Vec<T>>,
    evaluation: T,
    final_value: T,
}

impl<T: MyField> Proof<T> {
    #[deprecated(note = "use `CanonicalSerialize`")]
    #[allow(deprecated)]
    pub fn size(&self) -> usize {
        self.merkle_root.len() * MERKLE_ROOT_SIZE
            + self
                .query_result
                .0
                .iter()
                .fold(0, |acc, x| acc + x.proof_size())
            + self
                .query_result
                .1
                .iter()
                .fold(0, |acc, x| acc + x.proof_size())
            + (self.deep_evals.iter().fold(0, |acc, x| acc + x.1.len())
                + self.shuffle_evals.len()
                + self.out_evals.iter().fold(0, |acc, x| acc + x.len())
                + 2)
                * size_of::<T>()
    }
}

// Does not include `final_poly`, as with the proof size
impl<T> CanonicalSerialize for Proof<T>
where
    T: MyField + CanonicalSerialize,
{
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        for root in &self.merkle_root {
            root.serialize_with_mode(&mut writer, compress)?;
        }
        for q in &self.query_result.0 {
            q.serialize_with_mode(&mut writer, compress)?;
        }
        for q in &self.query_result.1 {
            q.serialize_with_mode(&mut writer, compress)?;
        }
        for (_, deep_eval) in &self.deep_evals {
            for d in deep_eval {
                d.serialize_with_mode(&mut writer, compress)?;
            }
        }
        for s in &self.shuffle_evals {
            s.serialize_with_mode(&mut writer, compress)?;
        }
        for o in &self.out_evals {
            for t in o {
                t.serialize_with_mode(&mut writer, compress)?;
            }
        }
        self.evaluation.serialize_with_mode(&mut writer, compress)?;
        self.final_value
            .serialize_with_mode(&mut writer, compress)?;
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.merkle_root
            .iter()
            .map(|r| r.serialized_size(compress))
            .sum::<usize>()
            + self
                .query_result
                .0
                .iter()
                .map(|q| q.serialized_size(compress))
                .sum::<usize>()
            + self
                .query_result
                .1
                .iter()
                .map(|q| q.serialized_size(compress))
                .sum::<usize>()
            + self
                .deep_evals
                .iter()
                .map(|d| {
                    d.1.iter()
                        .map(|d| d.serialized_size(compress))
                        .sum::<usize>()
                })
                .sum::<usize>()
            + self
                .shuffle_evals
                .iter()
                .map(|s| s.serialized_size(compress))
                .sum::<usize>()
            + self
                .out_evals
                .iter()
                .map(|s| s.iter().map(|t| t.serialized_size(compress)).sum::<usize>())
                .sum::<usize>()
            + self.evaluation.serialized_size(compress)
            + self.final_value.serialized_size(compress)
    }
}

#[cfg(test)]
mod tests {
    use ark_serialize::CanonicalSerialize;

    use crate::{prover::Prover, verifier::Verifier};
    use csv::Writer;
    use util::{
        algebra::{
            coset::Coset,
            field::{arkfield::Fp2256, MyField},
            polynomial::MultilinearPolynomial,
        },
        random_oracle::RandomOracle,
    };
    use util::{CODE_RATE, SECURITY_BITS, SIZE};

    #[allow(deprecated)]
    fn output_proof_size<T: MyField + CanonicalSerialize>(variable_num: usize) -> usize {
        let mut interpolate_cosets =
            vec![Coset::new(1 << (variable_num + CODE_RATE), T::from_int(1))];
        for i in 1..variable_num {
            interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
        }
        let oracle = RandomOracle::new(variable_num, SECURITY_BITS / CODE_RATE);
        let polynomials = (0..variable_num)
            .rev()
            .map(|x| MultilinearPolynomial::random_polynomial(x + 1))
            .collect::<Vec<MultilinearPolynomial<T>>>();
        let prover = Prover::new(variable_num, &interpolate_cosets, polynomials, &oracle);
        let commit = prover.commit_polynomial();
        let verifier = Verifier::new(variable_num, &interpolate_cosets, commit, &oracle);
        let point = verifier.get_open_point();
        let proof = prover.generate_proof(point);
        let size = proof.compressed_size();
        assert_eq!(size, proof.size());
        assert!(verifier.verify(proof));
        size
    }

    #[test]
    fn test_proof_size() {
        let mut wtr = Writer::from_path("batch.csv").unwrap();
        let range = 10..SIZE;
        for i in range.clone() {
            let proof_size = output_proof_size::<Fp2256>(i);
            wtr.write_record(&[i.to_string(), proof_size.to_string()])
                .unwrap();
        }
    }
}
