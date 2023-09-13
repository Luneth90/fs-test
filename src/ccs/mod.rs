use ark_ec::CurveGroup;
use ark_std::{log2, One, Zero};
use std::ops::Neg;
use thiserror::Error;

pub mod r1cs;
use r1cs::*;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Relation not satisfied")]
    NotSatisfied,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CCS<C: CurveGroup> {
    // number of rows in matrix
    pub m: usize,
    // number of cols in matrix, or n = |z|
    pub n: usize,
    // number of io
    pub l: usize,
    // number of matrices
    pub t: usize,
    // number of multisets
    pub q: usize,
    // max cardinality
    pub d: usize,
    // s = log(m), dimension of x
    pub s: usize,
    // s_prime = log(n), dimension of y
    pub s_prime: usize,

    // vector of matrices
    pub m_vec: Vec<SparseMatrix<C::ScalarField>>,
    // vector of multisets
    pub s_vec: Vec<Vec<usize>>,
    // vector of constants
    pub v: Vec<C::ScalarField>,
}

impl<C> CCS<C>
where
    C: CurveGroup,
{
    pub fn from_r1cs(r1cs: R1CS<C::ScalarField>) -> Self {
        CCS {
            m: r1cs.a.n_rows,
            n: r1cs.a.n_cols,
            l: r1cs.l,
            t: 3,
            q: 2,
            d: 2,
            s: log2(r1cs.a.n_rows) as usize,
            s_prime: log2(r1cs.a.n_cols) as usize,
            s_vec: vec![vec![0, 1], vec![2]],
            v: vec![C::ScalarField::one(), C::ScalarField::one().neg()],
            m_vec: vec![r1cs.a, r1cs.b, r1cs.c],
        }
    }

    pub fn is_satisfied(&self, z: &[C::ScalarField]) -> Result<(), Error> {
        let mut r = vec![C::ScalarField::zero(); self.m];
        for q_i in 0..self.q {
            let mut s_set = Vec::new();
            for s in &self.s_vec[q_i] {
                // s is multiset
                s_set.push(&self.m_vec[*s]);
            }
            // first each s * z, then hadamard each other in s
            let hadamard_vec = vec![C::ScalarField::one(); self.m];
            let s_z_vec: Vec<_> = s_set.iter().map(|s_m| vec_mul_matrix(z, s_m)).collect();
            let res = s_z_vec
                .iter()
                .fold(hadamard_vec, |acc, x| hadamard(&acc, x));
            // second multiply c
            let c_s = scalar_mul_vec(self.v[q_i], &res);
            // third add each other in r
            r = vec_add_vec(&r, &c_s);
        }
        for e in r {
            if !e.is_zero() {
                return Err(Error::NotSatisfied);
            }
        }
        //last check sum value = 0
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ccs::r1cs::tests::{get_test_r1cs, get_test_z};
    use ark_pallas::Projective;

    pub fn get_test_ccs<C: CurveGroup>() -> CCS<C> {
        let r1cs = get_test_r1cs::<C::ScalarField>();
        CCS::from_r1cs(r1cs)
    }
    #[test]
    fn test_ccs() {
        let ccs = get_test_ccs::<Projective>();
        let z = get_test_z(2);
        ccs.is_satisfied(&z).unwrap();
    }
}
