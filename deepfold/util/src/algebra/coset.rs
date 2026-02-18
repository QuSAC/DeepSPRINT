use super::{polynomial::Polynomial};
use crate::{algebra::field::FftField, batch_bit_reverse};
use ark_ff::{Field};

#[derive(Debug, Clone, Copy)]
struct Radix2Domain<T: Field> {
    order: usize,
    omega: T,
}

impl<T: Field> Radix2Domain<T> {
    #[inline]
    pub fn new(order: usize, omega: T) -> Self {
        Radix2Domain { order, omega }
    }

    #[inline]
    pub fn order(&self) -> usize {
        self.order
    }

    #[inline]
    pub fn omega(&self) -> T {
        self.omega
    }

    #[inline]
    pub fn fft(&self, a: &mut Vec<T>) {
        _fft(a, self.omega);
    }

    #[inline]
    pub fn ifft(&self, a: &mut Vec<T>) {
        _fft(a, self.omega.inverse().unwrap());
        let t = T::from(self.order as u64).inverse().unwrap();
        for i in a {
            *i *= t;
        }
    }

    #[inline]
    pub fn coset_fft(&self, a: &mut Vec<T>, shift: T) {
        multiply_by_coset(a, shift);
        self.fft(a);
    }

    #[inline]
    pub fn coset_ifft(&self, a: &mut Vec<T>, shift: T) {
        self.ifft(a);
        multiply_by_coset(a, shift.inverse().unwrap());
    }
}

fn multiply_by_coset<T: Field>(a: &mut Vec<T>, shift: T) {
    let mut t = shift;
    for i in 1..a.len() {
        a[i] *= t;
        t *= shift;
    }
}

fn _fft<T: Field>(a: &mut Vec<T>, omega: T) {
    let n = a.len();
    let log_n = n.ilog2() as usize;
    let rank = batch_bit_reverse(log_n);
    for i in 0..n {
        if i < rank[i] {
            (a[i], a[rank[i]]) = (a[rank[i]], a[i]);
        }
    }
    let mut log_m = 0usize;
    for _i in 0..log_n {
        let w_m = omega.pow(&[(n >> (log_m + 1)) as u64]);
        let m = 1 << log_m;
        for j in (0..n).step_by(m * 2) {
            let mut w = T::ONE;
            for k in 0..m {
                let t = w * a[j + k + m];
                a[j + k + m] = a[j + k] - t;
                a[j + k] += t;
                w *= w_m;
            }
        }
        log_m += 1;
    }
}

use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct Coset<T: FftField> {
    elements: Arc<Vec<T>>,
    elements_inv: Arc<Vec<T>>,
    fft_eval_domain: Radix2Domain<T>,
    shift: T,
}

impl<T: FftField> Coset<T> {
    pub fn mult(poly1: &Polynomial<T>, poly2: &Polynomial<T>) -> Polynomial<T> {
        let degree = {
            let max_d = std::cmp::max(poly1.degree(), poly2.degree()) + 1;
            let mut d = 1;
            while d < max_d {
                d <<= 1;
            }
            d << 1
        };
        let domain = Radix2Domain::new(degree, T::get_root_of_unity(degree as u32));
        let mut coeff1 = poly1.coefficients().clone();
        let len = coeff1.len();
        coeff1.append(&mut (len..degree).into_iter().map(|_| T::ZERO).collect());
        let mut coeff2 = poly2.coefficients().clone();
        let len = coeff2.len();
        coeff2.append(&mut (len..degree).into_iter().map(|_| T::ZERO).collect());
        domain.fft(&mut coeff1);
        domain.fft(&mut coeff2);
        for i in 0..degree {
            coeff1[i] *= coeff2[i];
        }
        domain.ifft(&mut coeff1);
        let poly = Polynomial::new(coeff1);
        poly
    }

    pub fn new(order: usize, shift: T) -> Self {
        assert!(!shift.is_zero());
        let omega = T::get_root_of_unity(order as u32);
        let elements = std::iter::successors(Some(shift), |&last| Some(last * omega))
            .take(order)
            .collect();
        let omega_inv = omega.pow(&[order as u64 - 1]);
        let elements_inv = std::iter::successors(Some(shift), |&last| Some(last * omega_inv))
            .take(order)
            .collect();
        Coset {
            elements: Arc::new(elements),
            elements_inv: Arc::new(elements_inv),
            fft_eval_domain: Radix2Domain::new(order, omega),
            shift,
        }
    }

    pub fn order(&self) -> usize {
        self.fft_eval_domain.order
    }

    pub fn pow(&self, index: usize) -> Coset<T> {
        assert_eq!(index & (index - 1), 0);
        let lowbit = (index as i64 & (-(index as i64))) as usize;
        Coset::new(self.order() / lowbit, self.shift.pow(&[index as u64]))
    }

    pub fn generator(&self) -> T {
        self.fft_eval_domain.omega
    }

    pub fn element_at(&self, index: usize) -> T {
        self.elements[index]
    }

    pub fn element_inv_at(&self, index: usize) -> T {
        self.elements_inv[index]
    }

    pub fn all_elements_inv(&self) -> Vec<T> {
        (*self.elements_inv).clone()
    }

    pub fn all_elements(&self) -> Vec<T> {
        (*self.elements).clone()
    }

    pub fn size(&self) -> usize {
        self.fft_eval_domain.order()
    }

    pub fn fft(&self, mut coeff: Vec<T>) -> Vec<T> {
        let n = self.size() - coeff.len();
        for _i in 0..n {
            coeff.push(T::ZERO);
        }
        self.fft_eval_domain.coset_fft(&mut coeff, self.shift);
        coeff
    }

    pub fn ifft(&self, mut evals: Vec<T>) -> Vec<T> {
        if evals.len() == 1 {
            return vec![evals[0]];
        };
        assert_eq!(self.size(), evals.len());
        self.fft_eval_domain.coset_ifft(&mut evals, self.shift);
        evals
    }

    pub fn shift(&self) -> T {
        self.shift
    }
}
