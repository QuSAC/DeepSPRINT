use std::hash::Hash;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::algebra::field::FftField;

use super::MyField;
use ark_ff::{Fp2ConfigWrapper};
use ark_ff::{
    fp2, AdditiveGroup, Field, Fp256, Fp2Config, MontBackend, MontConfig, MontFp, UniformRand,
};
use ark_serialize::{CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize, CanonicalSerializeWithFlags, Valid};
use num_traits::Zero;
use num_traits::One;

#[derive(MontConfig)]
#[modulus = "2261564242916331941866620800950935700259179388000792266395655937654553313279"]
#[generator = "3"]
pub struct FqConfig251;

pub type Field251 = Fp256<MontBackend<FqConfig251, 4>>;

pub type Fp2256 = fp2::Fp2<Fp2Config251>;
pub struct Fp2Config251;
impl Fp2Config for Fp2Config251 {
    type Fp = Field251;
    const NONRESIDUE: Self::Fp = MontFp!("2261564242916331941866620800950935700259179388000792266395655937654553313278");
    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        MontFp!("1"),
        MontFp!("2261564242916331941866620800950935700259179388000792266395655937654553313278")
    ];
}

impl FftField for Fp2256 {
    const TWO_ADICITY: u32 = 249;
    const TWO_ADIC_ROOT_OF_UNITY: Self = Fp2256 {
        c0: MontFp!("1205025138543176579400216426469194982307642063028969157762157265575227047065"),
        c1: MontFp!("1016267646470308916753499664310755376375900607926030331828540496222603030995")
    };
}

#[cfg(test)]
mod tests {
    use crate::algebra::field::test_fft_field_impl;
    #[test]
    fn test_two_adicity() {
        test_fft_field_impl::<super::Fp2256>();
    }
}
