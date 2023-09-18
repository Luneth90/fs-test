use ark_ff::PrimeField;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseMatrix<F: PrimeField> {
    //number of rows
    pub n_rows: usize,
    //number of cols
    pub n_cols: usize,
    // position and not zero value, others are zero
    pub vals: Vec<(usize, usize, F)>,
}

pub fn dense_matrix_to_sparse<F: PrimeField>(m: Vec<Vec<F>>) -> SparseMatrix<F> {
    let mut sm = SparseMatrix::<F> {
        n_rows: m.len(),
        n_cols: m[0].len(),
        vals: Vec::new(),
    };

    for (i, m_i) in m.iter().enumerate() {
        for (j, v) in m_i.iter().enumerate() {
            if !v.is_zero() {
                sm.vals.push((i, j, *v));
            }
        }
    }
    sm
}

pub fn to_f_vec<F: PrimeField>(v: Vec<usize>) -> Vec<F> {
    let mut f_v = Vec::new();
    for val in v.iter() {
        f_v.push(F::from(*val as u64));
    }
    f_v
}

pub fn to_f_matrix<F: PrimeField>(m: Vec<Vec<usize>>) -> Vec<Vec<F>> {
    let mut f_m = Vec::new();
    for m_i in m.iter() {
        let mut f_m_i = Vec::new();
        for m_ij in m_i.iter() {
            f_m_i.push(F::from(*m_ij as u64));
        }
        f_m.push(f_m_i);
    }
    f_m
}

pub fn vec_mul_matrix<F: PrimeField>(z: &[F], m: &SparseMatrix<F>) -> Vec<F> {
    let mut v = vec![F::zero(); m.n_rows];
    for (i, j, val) in &m.vals {
        v[*i] += *val * z[*j];
    }
    v
}

pub fn hadamard<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    a.iter().zip(b).map(|(v1, v2)| *v1 * v2).collect()
}

pub fn scalar_mul_vec<F: PrimeField>(c: F, v: &[F]) -> Vec<F> {
    v.iter().map(|a| c * a).collect()
}

pub fn vec_add_vec<F: PrimeField>(v1: &[F], v2: &[F]) -> Vec<F> {
    v1.iter().zip(v2).map(|(v1, v2)| *v1 + v2).collect()
}

pub fn vec_sub_vec<F: PrimeField>(v1: &[F], v2: &[F]) -> Vec<F> {
    v1.iter().zip(v2.iter()).map(|(v1, v2)| *v1 - v2).collect()
}

pub struct R1CS<F: PrimeField> {
    //io length
    pub l: usize,
    pub a: SparseMatrix<F>,
    pub b: SparseMatrix<F>,
    pub c: SparseMatrix<F>,
}

impl<F: PrimeField> R1CS<F> {
    // z = (1, x, w)
    // split z to (w, x) pair
    pub fn split_z(&self, z: &[F]) -> (Vec<F>, Vec<F>) {
        (z[self.l + 1..].to_vec(), z[1..self.l + 1].to_vec())
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    pub fn get_test_r1cs<F: PrimeField>() -> R1CS<F> {
        //R1CS: x^3 + x + 5 = y)
        // x * x = z1
        // z1 * x = z2
        // (z2 + x) * 1 = z3
        // (z3 + 5) * 1 = y
        // z = (1,x,y,z1,z2,z3)
        // when x =3, z = (1,3,35,9,27,30)
        let a = dense_matrix_to_sparse(to_f_matrix::<F>(vec![
            vec![0, 1, 0, 0, 0, 0],
            vec![0, 0, 0, 1, 0, 0],
            vec![0, 1, 0, 0, 1, 0],
            vec![5, 0, 0, 0, 0, 1],
        ]));
        let b = dense_matrix_to_sparse(to_f_matrix::<F>(vec![
            vec![0, 1, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0],
            vec![1, 0, 0, 0, 0, 0],
            vec![1, 0, 0, 0, 0, 0],
        ]));
        let c = dense_matrix_to_sparse(to_f_matrix::<F>(vec![
            vec![0, 0, 0, 1, 0, 0],
            vec![0, 0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 0, 1],
            vec![0, 0, 1, 0, 0, 0],
        ]));

        R1CS::<F> { l: 1, a, b, c }
    }

    //z = (1,x,w)
    pub fn get_test_z<F: PrimeField>(input: usize) -> Vec<F> {
        to_f_vec(vec![
            1,
            input,
            input * input * input + input + 5,
            input * input,
            input * input * input,
            input * input * input + input,
        ])
    }
}
