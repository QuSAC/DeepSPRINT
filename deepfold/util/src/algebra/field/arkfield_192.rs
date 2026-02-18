use std::cmp;
use std::hash::Hash;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::algebra::field::FftField;

use super::MyField;
use ark_ff::{Fp2ConfigWrapper, Fp384};
use ark_ff::{
    fp2, AdditiveGroup, Field, Fp2Config, MontBackend, MontConfig, MontFp, UniformRand,
};
use ark_serialize::{CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize, CanonicalSerializeWithFlags, Valid};
use num_traits::Zero;
use num_traits::One;

#[derive(MontConfig)]
#[modulus = "21240143965243898950369170053983666816800796950485279844440876913226209392447715095215502903023899837622322901024767"]
#[generator = "3"]
pub struct FqConfig384;

pub type Field384 = Fp384<MontBackend<FqConfig384, 6>>;

pub type Fp2384 = fp2::Fp2<Fp2Config384>;
pub struct Fp2Config384;
impl Fp2Config for Fp2Config384 {
    type Fp = Field384;
    const NONRESIDUE: Self::Fp =
        MontFp!("21240143965243898950369170053983666816800796950485279844440876913226209392447715095215502903023899837622322901024766");
    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        MontFp!("1"),
        MontFp!("21240143965243898950369170053983666816800796950485279844440876913226209392447715095215502903023899837622322901024766"),
    ];
}

impl FftField for Fp2384 {
    const TWO_ADICITY: u32 = 378;
    const TWO_ADIC_ROOT_OF_UNITY: Self = Self {
        c0: MontFp!("16398322775289304959099330544527976706163512667538289916097814970894100616623012276626434150986151271587366138159230"),
        c1: MontFp!("2751546644444135832395674673784949863711811716729994201225531213918079940388963099472589400274300703244729887648983")
    };
}


#[cfg(test)]
mod tests {
    use crate::algebra::field::test_fft_field_impl;
    #[test]
    fn test_two_adicity() {
        test_fft_field_impl::<super::Fp2384>();
    }
}
