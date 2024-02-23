//! Create Structured Reference String (SRS) from a random number or a specific secret s


use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::{One, PrimeField};


use ark_std::rand::RngCore;
// TODO: Check if this is secure

use ark_std::UniformRand;
use std::{cmp::max, iter};

/// Create srs from rng.
/// Why it is called unsafe?
/// Possible answer (not sure): Because it uses the rng to generate the srs.
/// Input is N1, N2, and a random number generator
pub fn unsafe_setup_from_rng<E: PairingEngine, R: RngCore>(
    big_n1: usize,
    big_n2: usize,
    rng: &mut R,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let tau = E::Fr::rand(rng);

    let size = max(big_n1 + 1, big_n2 + 1);

    let powers_of_tau: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(*p * tau))
        .take(size)
        .collect();
    // Start with the value E::Fr::one(), which is the identity element of the finite field `Fr`
    // and then keep multiplying by `tau` to get the powers of `tau`
    // The `take` method returns an iterator that yields the first `n` elements of the iterator

    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_tau
        .iter()
        .take(big_n1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_tau
        .iter()
        .take(big_n2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();

        (srs_g1, srs_g2)
}

/// Create srs from specific tau
/// Input is N1, N2, and tau
pub fn unsafe_setup_from_tau<E: PairingEngine>(
    big_n1: usize,
    big_n2: usize,
    tau: E::Fr,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let size = max(big_n1 + 1, big_n2 + 1);
    let powers_of_tau: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(*p * tau))
        .take(size)
        .collect();

    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_tau
        .iter()
        .take(big_n1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_tau
        .iter()
        .take(big_n2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();
    (srs_g1, srs_g2)
}
