use ark_ec::CurveGroup;
use ark_std::{One, Zero};
use crate::pedersen::Params as PedersenParams;

pub mod circuits;
pub mod nifs;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CommittedInstance<C: CurveGroup> {
    pub cm_e: C,
    pub u: C::ScalarField,
    pub cm_w: C,
    pub x: Vec<C::ScalarField>,
}

impl<C: CurveGroup> CommittedInstance<C> {
    pub fn empty() -> Self {
        CommittedInstance {
            cm_e: C::zero(),
            u: C::ScalarField::one(),
            cm_w: C::zero(),
            x: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Witness<C: CurveGroup> {
    pub e: Vec<C::ScalarField>,
    pub r_e: C::ScalarField,
    pub w: Vec<C::ScalarField>,
    pub r_w: C::ScalarField,
}

impl<C: CurveGroup> Witness<C> {
    pub fn new(witness: Vec<C::ScalarField>, e_len: usize) -> Self {
        Self {
            e: vec![C::ScalarField::zero(); e_len],
            r_e: C::ScalarField::one(),
            w: witness,
            r_w: C::ScalarField::one(),
        }
    }

    pub fn commit(
        &self,
        params: &PedersenParams<C>,
        x: Vec<C::ScalarField>,
    ) -> CommittedInstance<C> {
        unimplemented!();
    }
}
