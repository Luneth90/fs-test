use std::marker::PhantomData;

use ark_ec::{AffineRepr, CurveGroup};
use ark_r1cs_std::prelude::{Boolean, CurveVar};
use ark_relations::r1cs::SynthesisError;

pub type CF<C> = <<C as CurveGroup>::Affine as AffineRepr>::BaseField;

pub struct ECRLC<C: CurveGroup, GC: CurveVar<C, CF<C>>> {
    _c: PhantomData<C>,
    _gc: PhantomData<GC>,
}

impl<C: CurveGroup, GC: CurveVar<C, CF<C>>> ECRLC<C, GC> {
    pub fn check(
        r_bits: Vec<Boolean<CF<C>>>,
        p1: GC,
        p2: GC,
        p3: GC,
    ) -> Result<(), SynthesisError> {
        p3.enforce_equal(&(p1 + p2.scalar_mul_le(r_bits.iter())?))?;

        Ok(())
    }
}
