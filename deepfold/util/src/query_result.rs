use ark_serialize::CanonicalSerialize;

use crate::algebra::field::{as_bytes_vec, FftField};
use crate::merkle_tree::MerkleTreeVerifier;
use std::collections::HashMap;
use std::mem::size_of;

#[derive(Clone, Debug)]
pub struct QueryResult<T: FftField> {
    pub proof_bytes: Vec<u8>,
    // Cauchy: Why use hashmap rather than Vec here?
    pub proof_values: HashMap<usize, T>,
}

impl<T: FftField> QueryResult<T> {
    pub fn verify_merkle_tree(
        &self,
        leaf_indices: &Vec<usize>,
        leaf_size: usize,
        merkle_verifier: &MerkleTreeVerifier,
    ) -> bool {
        let len = merkle_verifier.leave_number;

        let leaves: Vec<Vec<u8>> = leaf_indices
            .iter()
            .map(|x| {
                as_bytes_vec(
                    &(0..leaf_size)
                        .map(|j| {
                            self.proof_values
                                .get(&(x.clone() + j * len))
                                .unwrap()
                                .clone()
                        })
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        let res = merkle_verifier.verify(self.proof_bytes.clone(), leaf_indices, &leaves);
        assert!(res);
        res
    }

    #[deprecated(note = "use `CanonicalSerialize`")]
    pub fn proof_size(&self) -> usize {
        self.proof_bytes.len() + self.proof_values.len() * size_of::<T>()
    }
}

impl<T> CanonicalSerialize for QueryResult<T>
where
    T: FftField + CanonicalSerialize,
{
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        for byte in &self.proof_bytes {
            byte.serialize_with_mode(&mut writer, compress)?;
        }
        // TODO: This is incorrect. We have evaluations at random points. The
        // evaluation points are Fiat-Shamir derived, so we don't need to
        // serialize them. We do however need to serialize them in some
        // canonical order. For now, we can't deserialize anyway since the
        // evaluation points are included in the struct.
        for val in self.proof_values.values() {
            val.serialize_with_mode(&mut writer, compress)?;
        }
        Ok(())
    }
    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.proof_bytes
            .iter()
            .map(|b| b.serialized_size(compress))
            .sum::<usize>()
            + self
                .proof_values
                .values()
                .map(|val| val.serialized_size(compress))
                .sum::<usize>()
    }
}
