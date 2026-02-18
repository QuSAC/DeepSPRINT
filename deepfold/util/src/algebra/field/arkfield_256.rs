use std::hash::Hash;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::algebra::field::FftField;

use super::MyField;
use ark_ff::{Fp2ConfigWrapper, Fp512};
use ark_ff::{
    fp2, AdditiveGroup, Field, Fp256, Fp2Config, MontBackend, MontConfig, MontFp, UniformRand,
};
use ark_serialize::{CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize, CanonicalSerializeWithFlags, Valid};
use halo2curves::ff_ext::inverse;
use num_traits::ops::mul_add;
use num_traits::Zero;
use num_traits::One;

#[derive(MontConfig)]
#[modulus = "864175120484581453683482079962486176185193500155369104423588921177379322250834082489183304374038697487834084609675858746433355728113743766078731283595263"]
#[generator = "3"]
pub struct FqConfig512;

pub type Field512 = Fp512<MontBackend<FqConfig512, 8>>;

pub struct Fp2Config512;
impl Fp2Config for Fp2Config512 {
    type Fp = Field512;
    const NONRESIDUE: Self::Fp =
        MontFp!("864175120484581453683482079962486176185193500155369104423588921177379322250834082489183304374038697487834084609675858746433355728113743766078731283595262");
    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        MontFp!("1"),
        MontFp!("864175120484581453683482079962486176185193500155369104423588921177379322250834082489183304374038697487834084609675858746433355728113743766078731283595262"),
    ];
}
pub type Fp2512 = fp2::Fp2<Fp2Config512>;

impl FftField for Fp2512 {
    const TWO_ADICITY: u32 = 504;
    const TWO_ADIC_ROOT_OF_UNITY: Self = Self {
        c0: MontFp!("230305430540098244123025424958320986723918644073698572553990584973970303622483635385703949484096938866916264392117681055820472987130474547896203739654023"),
        c1: MontFp!("648274971811642665659778556073232260884853363249672153380717486659113964826441401594664833035933488481267410003613514484575050416792521578809087988162073")
    };
}

#[cfg(test)]
mod tests {
    use crate::algebra::field::test_fft_field_impl;
    #[test]
    fn test_two_adicity() {
        test_fft_field_impl::<super::Fp2512>();
    }
}
