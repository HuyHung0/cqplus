use ark_ec::{self, AffineCurve, PairingEngine};
use ark_ff::{self, Field, One};

use crate::{
    data_structures::{Proof, VerificationKey},
    error::Error,
};

pub fn verify<E: PairingEngine>(
    vk: &VerificationKey<E>,
    proof: &Proof<E>,
    commit_poly_t2: E::G2Affine,
    value_vartheta: E::Fr,
) -> Result<(), Error> {
    // Check first equation
    let pair_1 = E::pairing(proof.commit_poly_a, commit_poly_t2);
    let pair_2_1 = proof
        .commit_poly_a
        .mul(proof.value_beta + proof.value_eta)
        .into()
        + proof.commit_poly_m.mul(-E::Fr::one()).into()
        + proof
            .commit_poly_s
            .mul(proof.value_eta * proof.value_eta)
            .into();
    let pair_2_2 = vk.commit_poly_u2;
    let pair_2 = E::pairing(pair_2_1, pair_2_2);

    let pair_3_1 = proof
        .commit_poly_b
        .mul(proof.value_eta / value_vartheta)
        .into();
    let pair_3_2 = vk.commit_poly_zu2;
    let pair_3 = E::pairing(pair_3_1, pair_3_2);

    let pair_4_1 = proof.commit_poly_q;
    let pair_4_2 = vk.commit_poly_vu2;
    let pair_4 = E::pairing(pair_4_1, pair_4_2);

    let rhs_1 = proof.commit_poly_r_c.mul(proof.value_eta).into();
    let rhs_2 = vk.commit_1;
    let rhs = E::pairing(rhs_1, rhs_2);

    let lhs = pair_1 * pair_2 * pair_3.inverse().unwrap() * pair_4.inverse().unwrap();

    if lhs != rhs {
        // dbg!("lhr is:{rhs}",lhs);
        // dbg!("rhs is:{rhs}",rhs);
        return Err(Error::VerificationEquation1Error);
    }

    // Check second equation
    let pair_1_1 = proof.commit_poly_b
        + proof.commit_poly_d.mul(proof.value_eta).into()
        + proof.commit_value_b_gamma.mul(-E::Fr::one()).into();
    let pair_1_2 = vk.commit_1;
    let pair_1 = E::pairing(pair_1_1, pair_1_2);

    let pair_2 = E::pairing(proof.commit_poly_p, proof.commit_s_minus_gamma);

    if pair_1 != pair_2 {
        return Err(Error::VerificationEquation2Error);
    }
    Ok(())
}
