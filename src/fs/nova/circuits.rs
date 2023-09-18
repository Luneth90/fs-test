use std::marker::PhantomData;

use ark_ec::CurveGroup;
use ark_r1cs_std::{
    fields::fp::FpVar,
    prelude::{AllocVar, AllocationMode, Boolean, CurveVar, EqGadget},
};
use ark_relations::r1cs::{Namespace, SynthesisError};

use crate::fs::circuits::cyclefold::ECRLC;

use super::CommittedInstance;

//pub type CF1<C> = <<C as CurveGroup>::Affine as AffineRepr>::ScalarField;
//pub type CF2<C> = <<C as CurveGroup>::Affine as AffineRepr>::BaseField;

#[derive(Debug, Clone)]
pub struct CommittedInstanceE1Var<C: CurveGroup> {
    u: FpVar<C::ScalarField>,
    x: Vec<FpVar<C::ScalarField>>,
}

impl<C: CurveGroup> AllocVar<CommittedInstance<C>, C::ScalarField> for CommittedInstanceE1Var<C> {
    fn new_variable<T: std::borrow::Borrow<CommittedInstance<C>>>(
        cs: impl Into<Namespace<C::ScalarField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        f().and_then(|ci| {
            let cs = cs.into();
            let u = FpVar::<C::ScalarField>::new_variable(cs.clone(), || Ok(ci.borrow().u), mode)?;
            let x = Vec::new_variable(cs, || Ok(ci.borrow().x.clone()), mode)?;
            Ok(Self { u, x })
        })
    }
}

#[derive(Debug, Clone)]
pub struct CommittedInstanceE2Var<C: CurveGroup, GC: CurveVar<C, C::BaseField>> {
    _c: PhantomData<C>,
    cm_e: GC,
    cm_w: GC,
}

impl<C, GC> AllocVar<CommittedInstance<C>, C::BaseField> for CommittedInstanceE2Var<C, GC>
where
    C: CurveGroup,
    GC: CurveVar<C, C::BaseField>,
{
    fn new_variable<T: std::borrow::Borrow<CommittedInstance<C>>>(
        cs: impl Into<Namespace<C::BaseField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        f().and_then(|ci| {
            let cs = cs.into();
            let cm_e = GC::new_variable(cs.clone(), || Ok(ci.borrow().cm_e), mode)?;
            let cm_w = GC::new_variable(cs, || Ok(ci.borrow().cm_w), mode)?;
            Ok(Self {
                _c: PhantomData,
                cm_e,
                cm_w,
            })
        })
    }
}

pub struct NIFSGadget<C: CurveGroup> {
    _c: PhantomData<C>,
}

impl<C: CurveGroup> NIFSGadget<C> {
    pub fn verify(
        r: FpVar<C::ScalarField>,
        ci1: CommittedInstanceE1Var<C>,
        ci2: CommittedInstanceE1Var<C>,
        ci3: CommittedInstanceE1Var<C>,
    ) -> Result<(), SynthesisError> {
        ci3.u.enforce_equal(&(ci1.u + &r * ci2.u))?;
        let ci3_x = ci1
            .x
            .iter()
            .zip(ci2.x)
            .map(|(v1, v2)| v1 + &r * &v2)
            .collect::<Vec<FpVar<C::ScalarField>>>();
        ci3.x.enforce_equal(&ci3_x)?;
        Ok(())
    }
}

pub struct NIFSCycleGadget<C: CurveGroup, GC: CurveVar<C, C::BaseField>> {
    _c: PhantomData<C>,
    _gc: PhantomData<GC>,
}

impl<C: CurveGroup, GC: CurveVar<C, C::BaseField>> NIFSCycleGadget<C, GC> {
    pub fn verify(
        r_bits: Vec<Boolean<C::BaseField>>,
        cm_t: GC,
        ci1: CommittedInstanceE2Var<C, GC>,
        ci2: CommittedInstanceE2Var<C, GC>,
        ci3: CommittedInstanceE2Var<C, GC>,
    ) -> Result<(), SynthesisError> {
        // cm(E) check: ci3.cmE == ci1.cmE + r * cmT + r^2 * ci2.cmE
        ci3.cm_e.enforce_equal(
            &(ci1.cm_e
                + cm_t.scalar_mul_le(r_bits.iter())?
                + ci2
                    .cm_e
                    .scalar_mul_le(r_bits.iter())?
                    .scalar_mul_le(r_bits.iter())?),
        )?;
        ECRLC::check(r_bits, ci1.cm_w, ci2.cm_w, ci3.cm_w)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::BigInteger;
    use ark_ff::PrimeField;
    use ark_pallas::constraints::GVar;
    use ark_pallas::{Fq, Fr, Projective};
    use ark_relations::r1cs::ConstraintSystem;

    use super::*;
    use crate::{
        ccs::r1cs::tests::{get_test_r1cs, get_test_z},
        fs::nova::{nifs::NIFS, Witness},
        pedersen::Pedersen,
    };
    use ark_std::UniformRand;

    #[test]
    fn test_nifs_gadget() {
        let r1cs = get_test_r1cs();
        let z1 = get_test_z(3);
        let z2 = get_test_z(4);
        let (w1, x1) = r1cs.split_z(&z1);
        let (w2, x2) = r1cs.split_z(&z2);
        let w1 = Witness::<Projective>::new(w1, r1cs.a.n_rows);
        let w2 = Witness::new(w2, r1cs.a.n_rows);

        let mut rng = ark_std::test_rng();
        let max = r1cs.a.n_cols;
        let params = Pedersen::new_params(&mut rng, max);

        let ci1 = w1.commit(&params, x1);
        let ci2 = w2.commit(&params, x2);

        let r = Fr::rand(&mut rng);
        let (_w3, ci3, _t, cm_t) = NIFS::prove(&params, r, &r1cs, &w1, &ci1, &w2, &ci2);

        let cs = ConstraintSystem::<Fr>::new_ref();
        let r_var = FpVar::<Fr>::new_witness(cs.clone(), || Ok(r)).unwrap();
        let ci1_var = CommittedInstanceE1Var::new_witness(cs.clone(), || Ok(ci1.clone())).unwrap();
        let ci2_var = CommittedInstanceE1Var::new_witness(cs.clone(), || Ok(ci2.clone())).unwrap();
        let ci3_var = CommittedInstanceE1Var::new_witness(cs.clone(), || Ok(ci3.clone())).unwrap();
        NIFSGadget::verify(r_var, ci1_var, ci2_var, ci3_var).unwrap();
        assert!(cs.is_satisfied().unwrap());

        let cs = ConstraintSystem::<Fq>::new_ref();
        let r_bits = BigInteger::to_bits_le(&Fr::into_bigint(r));
        let r_bits_var = Vec::<Boolean<Fq>>::new_witness(cs.clone(), || Ok(r_bits)).unwrap();
        let cm_t_var = GVar::new_witness(cs.clone(), || Ok(cm_t)).unwrap();

        let ci1_var =
            CommittedInstanceE2Var::<Projective, GVar>::new_witness(cs.clone(), || Ok(ci1.clone()))
                .unwrap();
        let ci2_var =
            CommittedInstanceE2Var::<Projective, GVar>::new_witness(cs.clone(), || Ok(ci2.clone()))
                .unwrap();
        let ci3_var =
            CommittedInstanceE2Var::<Projective, GVar>::new_witness(cs.clone(), || Ok(ci3.clone()))
                .unwrap();

        NIFSCycleGadget::verify(r_bits_var, cm_t_var, ci1_var, ci2_var, ci3_var).unwrap();
        assert!(cs.is_satisfied().unwrap());
    }
} /* test */
