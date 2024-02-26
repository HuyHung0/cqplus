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
    big_n1: usize, // Max power of srs1. Length of srs1 = big_n1 + 1
    big_n2: usize, // Max power of srs2. Length of srs2 = big_n2 + 1
    rng: &mut R,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    // We use the notation s as in the paper. In the repo cq, they use the notation tau.
    let s = E::Fr::rand(rng);

    let max_size = max(big_n1 + 1, big_n2 + 1);

    // Prepare the power of s: from 1=s^0 to s^max_size
    // Start with the value E::Fr::one(), which is the identity element of the finite field `E::Fr`
    // and then keep multiplying by `s` to get the powers of `s`
    // The `take` method returns an iterator that yields the first `max_size` elements of the iterator
    let powers_of_s: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(*p * s))
        .take(max_size)
        .collect();
    
    // Generate the generator of the group G1 and G2
    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_s
        .iter()
        .take(big_n1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_s
        .iter()
        .take(big_n2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();

        (srs_g1, srs_g2)
}

/// Create srs from specific s
/// Input is N1, N2, and s
pub fn unsafe_setup_from_s<E: PairingEngine>(
    big_n1: usize,
    big_n2: usize,
    s: E::Fr,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let size = max(big_n1 + 1, big_n2 + 1);
    let powers_of_s: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(*p * s))
        .take(size)
        .collect();

    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_s
        .iter()
        .take(big_n1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_s
        .iter()
        .take(big_n2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();
    (srs_g1, srs_g2)
}
