use ark_crypto_primitives::sponge::poseidon::{PoseidonConfig, PoseidonSponge};
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ec::AffineRepr;
use ark_ec::{CurveGroup, Group};
use ark_ff::BigInteger;
use ark_ff::Field;
use ark_ff::PrimeField;

use super::Transcript;

pub struct PoseidonTranscript<C: CurveGroup> {
    sponge: PoseidonSponge<C::ScalarField>,
}

impl<C: CurveGroup> Transcript<C> for PoseidonTranscript<C>
where
    <C as Group>::ScalarField: Absorb,
{
    type TranscriptConfig = PoseidonConfig<C::ScalarField>;

    fn new(poseidon_config: &Self::TranscriptConfig) -> Self {
        Self {
            sponge: PoseidonSponge::<C::ScalarField>::new(poseidon_config),
        }
    }

    fn absorb(&mut self, v: &C::ScalarField) {
        self.sponge.absorb(&v);
    }

    fn absorb_vec(&mut self, v: &[C::ScalarField]) {
        self.sponge.absorb(&v);
    }

    fn absorb_point(&mut self, p: &C) {
        self.sponge.absorb(&prepare_point(p));
    }

    fn get_challenge(&mut self) -> C::ScalarField {
        let c = self.sponge.squeeze_field_elements(1);
        self.sponge.absorb(&c[0]);
        c[0]
    }

    fn get_challenges(&mut self, n: usize) -> Vec<C::ScalarField> {
        let c = self.sponge.squeeze_field_elements(n);
        self.sponge.absorb(&c);
        c
    }
}

fn prepare_point<C: CurveGroup>(p: &C) -> Vec<C::ScalarField> {
    let p_affine = p.into_affine();
    let p_xy = &p_affine.xy().unwrap();
    let x = p_xy
        .0
        .to_base_prime_field_elements()
        .next()
        .expect("a")
        .into_bigint();
    let y = p_xy
        .1
        .to_base_prime_field_elements()
        .next()
        .expect("a")
        .into_bigint();
    vec![
        C::ScalarField::from_le_bytes_mod_order(&x.to_bytes_le()),
        C::ScalarField::from_le_bytes_mod_order(&y.to_bytes_le()),
    ]
}

#[cfg(test)]
pub mod tests {
    use ark_crypto_primitives::sponge::poseidon::{find_poseidon_ark_and_mds, PoseidonConfig};
    use ark_ff::PrimeField;
    use ark_pallas::{Fr, Projective};
    use super::*;

    pub fn poseidon_test_config<F: PrimeField>() -> PoseidonConfig<F> {
        let full_rounds = 8;
        let partial_rounds = 31;
        let alpha = 5;
        let rate = 2;

        let (ark, mds) = find_poseidon_ark_and_mds(
            F::MODULUS_BIT_SIZE as u64,
            rate,
            full_rounds,
            partial_rounds,
            0,
        );
        PoseidonConfig::new(
            full_rounds as usize,
            partial_rounds as usize,
            alpha,
            mds,
            ark,
            rate,
            1,
        )
    }

    #[test]
    fn test_transcript_challenge() {
        let config = poseidon_test_config::<Fr>(); 
        let mut tr = PoseidonTranscript::<Projective>::new(&config);
        tr.absorb(&Fr::from(42u32));
        let _c = tr.get_challenge();
    }
}
