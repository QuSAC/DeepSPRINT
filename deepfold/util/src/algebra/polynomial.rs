use std::ops::{Add, Mul};

use ark_ff::{Field};
use rand::Rng;

use crate::{algebra::field::FftField, batch_bit_reverse};

use super::coset::Coset;

#[derive(Debug, Clone)]
pub struct Polynomial<T: FftField> {
    pub coefficients: Vec<T>,
}

impl<T: FftField> Polynomial<T> {
    pub fn new(mut coefficients: Vec<T>) -> Polynomial<T> {
        let zero = T::ZERO;
        while *coefficients.last().unwrap() == zero {
            coefficients.pop();
        }
        Polynomial { coefficients }
    }

    pub fn coefficients(&self) -> &Vec<T> {
        &self.coefficients
    }

    pub fn random_polynomial(rng: &mut impl Rng, degree: usize) -> Polynomial<T> {
        Polynomial {
            coefficients: (0..degree).map(|_| T::rand(rng)).collect(),
        }
    }

    pub fn degree(&self) -> usize {
        let n = self.coefficients.len();
        if n == 0 {
            0
        } else {
            n - 1
        }
    }

    pub fn evaluation_at(&self, x: T) -> T {
        let mut res = T::ZERO;
        for i in self.coefficients.iter().rev() {
            res *= x;
            res += *i;
        }
        res
    }

    pub fn evaluation_over_coset(&self, coset: &Coset<T>) -> Vec<T> {
        coset.fft(self.coefficients.clone())
    }

    pub fn over_vanish_polynomial(
        &self,
        vanishing_polynomial: &VanishingPolynomial<T>,
    ) -> Polynomial<T> {
        let degree = vanishing_polynomial.degree;
        let low_term = vanishing_polynomial.shift;
        let mut coeff = vec![];
        let mut remnant = self.coefficients.clone();
        for i in (degree..self.coefficients.len()).rev() {
            let tmp = remnant[i] * low_term;
            coeff.push(remnant[i]);
            remnant[i - degree] += tmp;
        }
        coeff.reverse();
        Polynomial::new(coeff)
    }
}

#[derive(Debug, Clone)]
pub struct VanishingPolynomial<T: FftField> {
    degree: usize,
    shift: T,
}

impl<T: FftField> VanishingPolynomial<T> {
    pub fn new(coset: &Coset<T>) -> VanishingPolynomial<T> {
        let degree = coset.size();
        VanishingPolynomial {
            degree,
            shift: coset.shift().pow(&[degree as u64]),
        }
    }

    pub fn evaluation_at(&self, x: T) -> T {
        x.pow(&[self.degree as u64]) - self.shift
    }
}

#[derive(Debug, Clone)]
pub struct MultilinearPolynomial<T: Field> {
    coefficients: Vec<T>,
}

impl<T: Field> Add for MultilinearPolynomial<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let mut l = self.coefficients().clone();
        let mut r = rhs.coefficients.clone();
        if r.len() < l.len() {
            (l, r) = (r, l);
        }

        while l.len() < r.len() {
            l.extend(vec![T::ZERO; l.len()]);
        }

        assert_eq!(l.len(), r.len());
        Self::new(
            l
                .iter()
                .zip(r.iter())
                .map(|(&f, &g)| f + g)
                .collect::<Vec<T>>(),
        )
    }
}

impl<T: Field> Mul<T> for MultilinearPolynomial<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(
            self.coefficients()
                .iter()
                .map(|&f| f * rhs)
                .collect::<Vec<T>>(),
        )
    }
}

impl<T: Field> MultilinearPolynomial<T> {
    pub fn coefficients(&self) -> &Vec<T> {
        &self.coefficients
    }

    pub fn evaluate_from_hypercube(point: Vec<T>, mut poly_hypercube: Vec<T>) -> T {
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

    pub fn evaluate_hypercube(&self) -> Vec<T> {
        let log_n = self.variable_num();
        let n = self.coefficients.len();
        let rank = batch_bit_reverse(log_n);
        let mut res = self.coefficients.clone();
        for i in 0..n {
            if i < rank[i] {
                (res[i], res[rank[i]]) = (res[rank[i]], res[i]);
            }
        }
        for i in 0..log_n {
            let m = 1 << i;
            for j in (0..n).step_by(m * 2) {
                for k in 0..m {
                    let tmp = res[j + k];
                    res[j + k + m] += tmp;
                }
            }
        }
        res
    }

    pub const fn new(coefficients: Vec<T>) -> Self {
        let len = coefficients.len();
        if !len.is_power_of_two() {
            panic!("not a power of 2")
        }
        MultilinearPolynomial { coefficients }
    }

    pub fn folding(&self, parameter: T) -> Self {
        let coefficients = Self::folding_vector(&self.coefficients, parameter);
        MultilinearPolynomial { coefficients }
    }

    pub fn fold_self(&mut self, parameter: T) {
        self.coefficients = Self::folding_vector(&self.coefficients, parameter);
    }

    pub fn add_mult(&mut self, poly: &MultilinearPolynomial<T>, k: T) {
        assert_eq!(self.coefficients.len(), poly.coefficients.len());
        for (i, j) in self.coefficients.iter_mut().zip(poly.coefficients.iter()) {
            *i += k * j.clone();
        }
    }

    fn folding_vector(v: &Vec<T>, parameter: T) -> Vec<T> {
        let len = v.len();
        assert_eq!(len & (len - 1), 0);
        let mut res = vec![];
        for i in (0..v.len()).step_by(2) {
            res.push(v[i] + parameter * v[i + 1]);
        }
        res
    }

    pub fn random_polynomial(rng: &mut impl Rng, variable_num: usize) -> Self {
        MultilinearPolynomial {
            coefficients: (0..(1 << variable_num))
                .map(|_| T::rand(rng))
                .collect(),
        }
    }

    pub fn evaluate(&self, point: &Vec<T>) -> T {
        let len = self.coefficients.len();
        assert_eq!(1 << point.len(), self.coefficients.len());
        let mut res = self.coefficients.clone();
        for (index, coeff) in point.iter().enumerate() {
            for i in (0..len).step_by(2 << index) {
                let x = *coeff * res[i + (1 << index)];
                res[i] += x;
            }
        }
        res[0]
    }

    pub fn evaluate_as_polynomial(&self, point: T) -> T {
        let mut res = T::ZERO;
        for i in self.coefficients.iter().rev() {
            res *= point;
            res += *i;
        }
        res
    }

    pub fn variable_num(&self) -> usize {
        self.coefficients.len().ilog2() as usize
    }
}

pub struct EqMultilinear<T: Field> {
    b: Vec<T>,
}

impl<T: Field> EqMultilinear<T> {
    pub fn evaluate_hypercube(&self) -> Vec<T> {
        let mut stack = vec![T::ONE];
        for b in self.b.iter() {
            let new_stack = stack
                .iter()
                .flat_map(|prod| {
                    [
                        prod.clone() * (T::ONE - b.clone()),
                        prod.clone() * b.clone(),
                    ]
                })
                .collect();
            stack = new_stack;
        }
        stack
    }

    pub fn new(b: Vec<T>) -> Self {
        EqMultilinear { b }
    }

    pub fn evaluate(&self, point: &Vec<T>) -> T {
        let mut res = T::ONE;
        for (x, b) in point.iter().zip(self.b.iter()) {
            res *=
                b.clone() * x.clone() + (T::ONE - b.clone()) * (T::ONE - x.clone());
        }
        res
    }
}
