use ark_ec::CurveGroup;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use std::marker::PhantomData;

use crate::{
    ccs::r1cs::{scalar_mul_vec, vec_add_vec},
    transcript::Transcript,
};

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Params<C: CurveGroup> {
    pub h: C,
    pub generators: Vec<C::Affine>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Proof<C: CurveGroup> {
    pub r_commit: C,
    pub u: Vec<C::ScalarField>,
    pub r_u: C::ScalarField,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Pedersen<C: CurveGroup> {
    _c: PhantomData<C>,
}

impl<C: CurveGroup> Pedersen<C> {
    pub fn new_params<R: Rng>(rng: &mut R, max: usize) -> Params<C> {
        let g = std::iter::repeat_with(|| C::Affine::rand(rng))
            .take(max.next_power_of_two())
            .collect();

        Params {
            h: C::rand(rng),
            generators: g,
        }
    }

    pub fn commit(r: &C::ScalarField, params: &Params<C>, v: &Vec<C::ScalarField>) -> C {
        //h*r + <g, v>
        params.h.mul(r) + C::msm(&params.generators[..v.len()], v).unwrap()
    }

    pub fn prove(
        cm: &C,
        v: &Vec<C::ScalarField>,
        r: &C::ScalarField,
        params: &Params<C>,
        transcript: &mut impl Transcript<C>,
    ) -> Proof<C> {
        transcript.absorb_point(cm);
        let r1 = transcript.get_challenge();
        let d = transcript.get_challenges(v.len());
        // r_commit = h*r1 + <g,d>
        let r_commit = params.h.mul(r1) + C::msm(&params.generators[..d.len()], &d).unwrap();
        transcript.absorb_point(&r_commit);
        let e = transcript.get_challenge();
        // u = d + v*e
        let u = vec_add_vec(&d, &scalar_mul_vec(e, &v));
        //r_u = r1 + e*r
        let r_u = r1 + e * r;
        Proof { r_commit, u, r_u }
    }

    pub fn verify(
        cm: C,
        proof: Proof<C>,
        params: &Params<C>,
        transcript: &mut impl Transcript<C>,
    ) -> bool {
        transcript.absorb_point(&cm);
        transcript.get_challenge(); //r1
        transcript.get_challenges(proof.u.len()); //d
        transcript.absorb_point(&proof.r_commit);
        let e = transcript.get_challenge();
        // r_commit + cm*e = h*r_u + <g,u>
        let lhs = proof.r_commit + cm.mul(e);
        let rhs = params.h.mul(proof.r_u)
            + C::msm(&params.generators[..proof.u.len()], &proof.u).unwrap();
        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use crate::transcript::poseidon::tests::poseidon_test_config;
    use crate::transcript::poseidon::PoseidonTranscript;
    use ark_pallas::{Fr, Projective};

    use super::*;

    #[test]
    fn test_pedersen_vec() {
        let mut rng = ark_std::test_rng();
        const MAX: usize = 10;
        let params = Pedersen::<Projective>::new_params(&mut rng, MAX);
        let poseidon_config = poseidon_test_config::<Fr>();

        let mut ts_prove = PoseidonTranscript::<Projective>::new(&poseidon_config);
        let mut ts_verify = PoseidonTranscript::<Projective>::new(&poseidon_config);
        let v = vec![Fr::rand(&mut rng); MAX];
        let r = Fr::rand(&mut rng);
        let cm = Pedersen::<Projective>::commit(&r, &params, &v);
        let proof = Pedersen::<Projective>::prove(&cm, &v, &r, &params, &mut ts_prove);
        let verify = Pedersen::<Projective>::verify(cm, proof, &params, &mut ts_verify);
        assert!(verify);
    }
} /* tests */
