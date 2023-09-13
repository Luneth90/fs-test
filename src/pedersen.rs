use ark_ec::CurveGroup;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use std::marker::PhantomData;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Params<C: CurveGroup> {
    pub h: C,
    pub generators: Vec<C::Affine>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Prove<C: CurveGroup> {
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

    pub fn commit(r: C::ScalarField, params: &Params<C>, v: &Vec<C::ScalarField>) -> C {
        //h*r + <g, v>
        params.h.mul(r) + C::msm(&params.generators[..v.len()], v).unwrap()
    }

    pub fn prove(cm: &C) -> Prove<C> {

        unimplemented!()
    }
}

