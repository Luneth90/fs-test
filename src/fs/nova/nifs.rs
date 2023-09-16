use std::marker::PhantomData;

use ark_ec::CurveGroup;
use ark_std::One;

use crate::{
    ccs::r1cs::{hadamard, scalar_mul_vec, vec_add_vec, vec_mul_matrix, vec_sub_vec, R1CS},
    pedersen::{Params as PedersenParams, Pedersen, Proof as PedersenProof},
    transcript::Transcript,
};

use super::{CommittedInstance, Witness};

pub struct NIFS<C: CurveGroup> {
    _phantom: PhantomData<C>,
}

impl<C: CurveGroup> NIFS<C> {
    pub fn compute_t(
        r1cs: &R1CS<C::ScalarField>,
        u1: C::ScalarField,
        u2: C::ScalarField,
        z1: &[C::ScalarField],
        z2: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let (a, b, c) = (&r1cs.a, &r1cs.b, &r1cs.c);

        let az1 = vec_mul_matrix(z1, a);
        let az2 = vec_mul_matrix(z2, a);
        let bz1 = vec_mul_matrix(z1, b);
        let bz2 = vec_mul_matrix(z2, b);
        let cz1 = vec_mul_matrix(z1, c);
        let cz2 = vec_mul_matrix(z2, c);

        let az1_bz2 = hadamard(&az1, &bz2);
        let az2_bz1 = hadamard(&az2, &bz1);
        let u1cz2 = scalar_mul_vec(u1, &cz2);
        let u2cz1 = scalar_mul_vec(u2, &cz1);

        vec_sub_vec(
            &vec_sub_vec(&vec_add_vec(&az1_bz2, &az2_bz1), &u1cz2),
            &u2cz1,
        )
    }

    pub fn fold_witness(
        w1: &Witness<C>,
        w2: &Witness<C>,
        t: &[C::ScalarField],
        r: C::ScalarField,
        r_t: C::ScalarField,
    ) -> Witness<C> {
        let r2 = r * r;
        let e = vec_add_vec(
            &vec_add_vec(&w1.e, &scalar_mul_vec(r, t)),
            &scalar_mul_vec(r2, &w2.e),
        );
        let r_e = w1.r_e + r * r_t + r2 * w2.r_e;
        let w = vec_add_vec(&w1.w, &scalar_mul_vec(r, &w2.w));
        let r_w = w1.r_w + r * w2.r_w;
        Witness { e, r_e, w, r_w }
    }

    pub fn fold_committed_instance(
        r: C::ScalarField,
        cm_t: &C,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
    ) -> CommittedInstance<C> {
        let r2 = r * r;
        let cm_e = ci1.cm_e + cm_t.mul(r) + ci2.cm_e.mul(r2);
        let u = ci1.u + r * ci2.u;
        let cm_w = ci1.cm_w + ci2.cm_w.mul(r);
        let x = vec_add_vec(&ci1.x, &scalar_mul_vec(r, &ci2.x));
        CommittedInstance { cm_e, u, cm_w, x }
    }

    ///Call fold method to generate new (w,ci,t,cm_t)
    pub fn prove(
        params: &PedersenParams<C>,
        r: C::ScalarField,
        r1cs: &R1CS<C::ScalarField>,
        w1: &Witness<C>,
        ci1: &CommittedInstance<C>,
        w2: &Witness<C>,
        ci2: &CommittedInstance<C>,
    ) -> (Witness<C>, CommittedInstance<C>, Vec<C::ScalarField>, C) {
        let u1 = ci1.u;
        let u2 = ci2.u;
        let z1 = [vec![ci1.u], ci1.x.to_vec(), w1.w.to_vec()].concat();
        let z2 = [vec![ci2.u], ci2.x.to_vec(), w2.w.to_vec()].concat();
        let t = Self::compute_t(r1cs, u1, u2, &z1, &z2);
        //r_t = 1, because cm_t do not need hiding property
        let r_t = C::ScalarField::one();
        let cm_t = Pedersen::commit(&r_t, params, &t);

        let w = Self::fold_witness(w1, w2, &t, r, r_t);
        let ci = Self::fold_committed_instance(r, &cm_t, ci1, ci2);
        (w, ci, t, cm_t)
    }

    ///Just generate ci
    pub fn verify(
        r: C::ScalarField,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        cm_t: &C,
    ) -> CommittedInstance<C> {
        Self::fold_committed_instance(r, cm_t, &ci1, &ci2)
    }

    ///Just verify fold method
    pub fn verify_fold_instance(
        r: C::ScalarField,
        ci: &CommittedInstance<C>,
        ci1: &CommittedInstance<C>,
        ci2: &CommittedInstance<C>,
        cm_t: &C,
    ) -> bool {
        let r2 = r * r;
        if ci.cm_e != ci1.cm_e + cm_t.mul(r) + ci2.cm_e.mul(r2) {
            return false;
        }

        if ci.u != ci1.u + r * ci2.u {
            return false;
        }

        if ci.cm_w != ci1.cm_w + ci2.cm_w.mul(r) {
            return false;
        }

        if ci.x != vec_add_vec(&ci1.x, &scalar_mul_vec(r, &ci2.x)) {
            return false;
        }

        return true;
    }

    /// use pedersen commitment to getnerate proof
    pub fn prove_commitments(
        ts: &mut impl Transcript<C>,
        params: &PedersenParams<C>,
        w: &Witness<C>,
        ci: &CommittedInstance<C>,
        t: &Vec<C::ScalarField>,
        cm_t: &C,
    ) -> (PedersenProof<C>, PedersenProof<C>, PedersenProof<C>) {
        let cm_w_proof = Pedersen::prove(&ci.cm_w, &w.w, &w.r_w, params, ts);
        let cm_e_proof = Pedersen::prove(&ci.cm_e, &w.e, &w.r_e, params, ts);
        let cm_t_proof = Pedersen::prove(cm_t, t, &C::ScalarField::one(), params, ts);
        (cm_t_proof, cm_w_proof, cm_e_proof)
    }

    /// finaly verify pedersen proof
    pub fn verify_commitments(
        ts: &mut impl Transcript<C>,
        params: &PedersenParams<C>,
        ci: &CommittedInstance<C>,
        cm_t: C,
        cm_t_proof: PedersenProof<C>,
        cm_w_proof: PedersenProof<C>,
        cm_e_proof: PedersenProof<C>,
    ) -> bool {
        if !Pedersen::verify(ci.cm_w, cm_w_proof, params, ts) {
            return false;
        }
        if !Pedersen::verify(ci.cm_e, cm_e_proof, params, ts) {
            return false;
        }
        if !Pedersen::verify(cm_t, cm_t_proof, params, ts) {
            return false;
        }
        return true;
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use ark_pallas::{Fr, Projective};

    use crate::{
        ccs::r1cs::tests::{get_test_r1cs, get_test_z},
        transcript::poseidon::{tests::poseidon_test_config, PoseidonTranscript},
    };
    use ark_std::UniformRand;

    use super::*;

    pub fn check_relaxed_r1cs<F: PrimeField>(r1cs: &R1CS<F>, z: Vec<F>, u: F, e: &[F]) {
        let az = vec_mul_matrix(&z, &r1cs.a);
        let bz = vec_mul_matrix(&z, &r1cs.b);
        let cz = vec_mul_matrix(&z, &r1cs.c);
        assert!(hadamard(&az, &bz) == vec_add_vec(&e, &scalar_mul_vec(u, &cz)));
    }

    #[test]
    fn test_nifs_fold_one() {
        let r1cs = get_test_r1cs();
        let z1 = get_test_z(3);
        let z2 = get_test_z(4);
        let (w1, x1) = r1cs.split_z(&z1);
        let (w2, x2) = r1cs.split_z(&z2);
        let w1 = Witness::<Projective>::new(w1.clone(), r1cs.a.n_rows);
        let w2 = Witness::new(w2.clone(), r1cs.a.n_rows);

        let mut rng = ark_std::test_rng();
        let params = Pedersen::new_params(&mut rng, r1cs.a.n_cols);

        let r = Fr::rand(&mut rng);
        let ci1 = w1.commit(&params, x1);
        let ci2 = w2.commit(&params, x2);

        let (w, _, t, cm_t) = NIFS::prove(&params, r, &r1cs, &w1, &ci1, &w2, &ci2);
        //nifs verify
        let ci = NIFS::verify(r, &ci1, &ci2, &cm_t);

        //check relaxed r1cs relation
        let z = [vec![ci.u], ci.x.to_vec(), w.w.to_vec()].concat();
        let z_aux = vec_add_vec(&z1, &scalar_mul_vec(r, &z2));
        assert_eq!(z, z_aux);

        check_relaxed_r1cs(&r1cs, z1, ci1.u, &w1.e);
        check_relaxed_r1cs(&r1cs, z2, ci2.u, &w2.e);
        check_relaxed_r1cs(&r1cs, z, ci.u, &w.e);

        let ci_expected = w.commit(&params, ci.x.clone());
        assert_eq!(ci_expected.cm_e, ci.cm_e);
        assert!(NIFS::verify_fold_instance(r, &ci, &ci1, &ci2, &cm_t));

        //generate pedersen commitment
        let config = poseidon_test_config();
        let mut ts_prove = PoseidonTranscript::new(&config);
        let mut ts_verify = PoseidonTranscript::new(&config);
        let (cm_t_proof, cm_w_proof, cm_e_proof) =
            NIFS::prove_commitments(&mut ts_prove, &params, &w, &ci, &t, &cm_t);
        let v = NIFS::verify_commitments(
            &mut ts_verify,
            &params,
            &ci,
            cm_t,
            cm_t_proof,
            cm_w_proof,
            cm_e_proof,
        );
        assert!(v);
    }

    #[test]
    fn test_nifs_fold_loop() {
        let r1cs = get_test_r1cs();
        let mut z1 = get_test_z(3);
        let (w1, x1) = r1cs.split_z(&z1);

        let mut rng = ark_std::test_rng();
        let params = Pedersen::new_params(&mut rng, r1cs.a.n_cols);

        let mut w1 = Witness::<Projective>::new(w1.clone(), r1cs.a.n_rows);
        let mut ci1 = w1.commit(&params, x1);
        let mut t1 = Vec::new();
        let mut cm_t1 = ci1.cm_w.clone();
        check_relaxed_r1cs(&r1cs, z1.clone(), ci1.u, &w1.e);

        let n = 10;
        for i in 0..n {
            let z2 = get_test_z(i + 4);
            let (w2, x2) = r1cs.split_z(&z2);
            let w2 = Witness::<Projective>::new(w2.clone(), r1cs.a.n_rows);
            let ci2 = w2.commit(&params, x2);
            check_relaxed_r1cs(&r1cs, z2.clone(), ci2.u, &w2.e);

            let r = Fr::rand(&mut rng);
            let (w3, _, t, cm_t) = NIFS::prove(&params, r, &r1cs, &w1, &ci1, &w2, &ci2);
            //nifs verify
            let ci3 = NIFS::verify(r, &ci1, &ci2, &cm_t);
            //
            //check relaxed r1cs relation
            let z3 = [vec![ci3.u], ci3.x.to_vec(), w3.w.to_vec()].concat();
            let z_aux = vec_add_vec(&z1, &scalar_mul_vec(r, &z2));
            assert_eq!(z3, z_aux);
            check_relaxed_r1cs(&r1cs, z3.clone(), ci3.u, &w3.e);

            let ci_expected = w3.commit(&params, ci3.x.clone());
            assert_eq!(ci_expected.cm_e, ci3.cm_e);
            assert!(NIFS::verify_fold_instance(r, &ci3, &ci1, &ci2, &cm_t));

            z1 = z3;
            w1 = w3;
            ci1 = ci3;
            cm_t1 = cm_t;
            t1 = t;
        }
        //generate pedersen commitment
        let config = poseidon_test_config();
        let mut ts_prove = PoseidonTranscript::new(&config);
        let mut ts_verify = PoseidonTranscript::new(&config);
        let (cm_t_proof, cm_w_proof, cm_e_proof) =
            NIFS::prove_commitments(&mut ts_prove, &params, &w1, &ci1, &t1, &cm_t1);
        let v = NIFS::verify_commitments(
            &mut ts_verify,
            &params,
            &ci1,
            cm_t1,
            cm_t_proof,
            cm_w_proof,
            cm_e_proof,
        );
        assert!(v);
    }
}
