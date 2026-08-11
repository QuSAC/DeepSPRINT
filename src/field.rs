use ark_ff::{BigInt, BigInteger, Fp2, Fp2Config, Fp256, MontBackend, MontFp, PrimeField};
use curve25519_dalek::Scalar;

#[derive(ark_ff::MontConfig)]
#[modulus = "7237005577332262213973186563042994240857116359379907606001950938285454250989"]
#[generator = "2"]
pub struct FpConfig253;

pub type Fp253 = Fp256<MontBackend<FpConfig253, 4>>;
pub type Fp2253 = Fp2<Fp2Config253>;

pub struct Fp2Config253;
impl Fp2Config for Fp2Config253 {
    type Fp = Fp253;
    const NONRESIDUE: Self::Fp = MontFp!("2");
    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        MontFp!("1"),
        MontFp!("7237005577332262213973186563042994240857116359379907606001950938285454250988"),
    ];
}

pub fn fp_to_scalar(fp: Fp253) -> Scalar {
    let bytes = fp.into_bigint().to_bytes_le();
    let byte_array: [u8; 32] = bytes.try_into().expect("Field element should be 32 bytes");
    Scalar::from_bytes_mod_order(byte_array)
}

pub fn fp2_to_scalar(f: Fp2253) -> (Scalar, Scalar) {
    (fp_to_scalar(f.c0), fp_to_scalar(f.c1))
}

pub fn scalar_to_fp(s: Scalar) -> Fp253 {
    let limbs: [u64; 4] = s
        .to_bytes()
        .chunks(64 / 8)
        .map(|chunk| {
            let mut limb_bytes = [0u8; 8];
            limb_bytes.copy_from_slice(chunk);
            u64::from_le_bytes(limb_bytes)
        })
        .collect::<Vec<u64>>()
        .try_into()
        .expect("Scalar should have 4 limbs");
    Fp253::from_bigint(BigInt(limbs)).unwrap()
}

pub fn scalar_to_fp2(s1: Scalar, s2: Scalar) -> Fp2253 {
    Fp2253::new(scalar_to_fp(s1), scalar_to_fp(s2))
}

#[cfg(test)]
mod tests {
    use ark_ff::Field;

    use super::*;

    #[test]
    fn test_fp_to_scalar_and_back() {
        let original_fp = Fp253::from(123456789u64);
        let scalar = fp_to_scalar(original_fp);
        let converted_fp = scalar_to_fp(scalar);
        assert_eq!(original_fp, converted_fp);
    }

    #[test]
    fn test_fp2_to_scalar_and_back() {
        let original_fp2 = Fp2253::new(Fp253::from(123456789u64), Fp253::from(987654321u64));
        let (scalar1, scalar2) = fp2_to_scalar(original_fp2);
        let converted_fp2 = scalar_to_fp2(scalar1, scalar2);
        assert_eq!(original_fp2, converted_fp2);
    }

    #[test]
    fn test_multiplication_conversion() {
        let a = Fp253::from(123456789u64);
        let b = Fp253::from(987654321u64);
        let product_fp = a * b;

        let scalar_a = fp_to_scalar(a);
        let scalar_b = fp_to_scalar(b);
        let product_scalar = scalar_a * scalar_b;

        let converted_product_fp = scalar_to_fp(product_scalar);
        assert_eq!(product_fp, converted_product_fp);
    }

    #[test]
    fn test_constant_multiplication_and_conversion_fp2() {
        let a = Fp2253::new(Fp253::from(123456789u64), Fp253::from(987654321u64));
        let b = Fp2253::new(Fp253::from(111111111u64), Fp253::from(0));
        let product_fp2 = a * b;

        let (scalar_a1, scalar_a2) = fp2_to_scalar(a);
        let (scalar_b1, scalar_b2) = fp2_to_scalar(b);
        assert_eq!(scalar_b2, Scalar::ZERO); // Ensure the second component of b is zero
        let product_scalar1 = scalar_a1 * scalar_b1;
        let product_scalar2 = scalar_a2 * scalar_b1;

        let converted_product_fp2 = scalar_to_fp2(product_scalar1, product_scalar2);
        assert_eq!(product_fp2, converted_product_fp2);
    }

    #[test]
    fn test_nonresidue() {
        assert!(Fp2Config253::NONRESIDUE.sqrt().is_none());
    }
}
