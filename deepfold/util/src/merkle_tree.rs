use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Valid};
use rs_merkle::{Hasher, MerkleProof, MerkleTree};

#[derive(Debug, Clone)]
pub struct Blake3Algorithm {}

impl Hasher for Blake3Algorithm {
    type Hash = [u8; MERKLE_ROOT_SIZE];

    fn hash(data: &[u8]) -> [u8; MERKLE_ROOT_SIZE] {
        blake3::hash(data).into()
    }
}

pub const MERKLE_ROOT_SIZE: usize = 32;
#[derive(Clone)]
pub struct MerkleTreeProver {
    pub merkle_tree: MerkleTree<Blake3Algorithm>,
    leave_num: usize,
}

#[derive(Debug, Clone)]
pub struct MerkleTreeVerifier {
    pub merkle_root: MerkleRoot,
    pub leave_number: usize,
}

impl MerkleTreeProver {
    pub fn new(leaf_values: Vec<Vec<u8>>) -> Self {
        let leaves = leaf_values
            .iter()
            .map(|x| Blake3Algorithm::hash(x))
            .collect::<Vec<_>>();
        let merkle_tree = MerkleTree::<Blake3Algorithm>::from_leaves(&leaves);
        Self {
            merkle_tree,
            leave_num: leaf_values.len(),
        }
    }

    pub fn leave_num(&self) -> usize {
        self.leave_num
    }

    pub fn commit(&self) -> MerkleRoot {
        MerkleRoot(self.merkle_tree.root().unwrap())
    }

    pub fn open(&self, leaf_indices: &Vec<usize>) -> Vec<u8> {
        self.merkle_tree.proof(leaf_indices).to_bytes()
    }
}

impl MerkleTreeVerifier {
    pub fn new(leave_number: usize, merkle_root: &MerkleRoot) -> Self {
        Self {
            leave_number,
            merkle_root: merkle_root.clone(),
        }
    }

    pub fn verify(
        &self,
        proof_bytes: Vec<u8>,
        indices: &Vec<usize>,
        leaves: &Vec<Vec<u8>>,
    ) -> bool {
        let proof = MerkleProof::<Blake3Algorithm>::try_from(proof_bytes).unwrap();
        let leaves_to_prove: Vec<[u8; MERKLE_ROOT_SIZE]> =
            leaves.iter().map(|x| Blake3Algorithm::hash(x)).collect();
        proof.verify(
            self.merkle_root.0,
            indices,
            &leaves_to_prove,
            self.leave_number,
        )
    }
}

#[derive(Debug, Clone)]
pub struct MerkleRoot([u8; MERKLE_ROOT_SIZE]);
impl MerkleRoot {
    pub fn get_root(
        proof_bytes: Vec<u8>,
        index: usize,
        leaf: Vec<u8>,
        total_leaves_count: usize,
    ) -> Self {
        let proof = MerkleProof::<Blake3Algorithm>::try_from(proof_bytes).unwrap();
        let leaf_hashes = vec![Blake3Algorithm::hash(&leaf)];
        Self(
            proof
                .root(&vec![index], &leaf_hashes, total_leaves_count)
                .unwrap(),
        )
    }
}

impl CanonicalSerialize for MerkleRoot {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        for b in self.0 {
            // It is probably overkill to serialize each byte individually, but we want to prevent the size to be included.
            b.serialize_with_mode(&mut writer, compress)?;
        }
        Ok(())
    }
    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.0.iter().map(|b| b.serialized_size(compress)).sum()
    }
}

impl Valid for MerkleRoot {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

impl CanonicalDeserialize for MerkleRoot {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let mut root = MerkleRoot([0; MERKLE_ROOT_SIZE]);
        for b in root.0.iter_mut() {
            *b = u8::deserialize_with_mode(&mut reader, compress, validate)?;
        }
        Ok(root)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::as_bytes_vec;

    #[test]
    fn blake3() {
        let hash_res = Blake3Algorithm::hash("data".as_bytes());
        let hex_string = hex::encode(hash_res);
        assert_eq!(
            "28a249c2e4d3a92bc0a16ed8f1b5cf83ca20415ee12e502b096624902bbc97bd",
            hex_string
        );
    }
}
