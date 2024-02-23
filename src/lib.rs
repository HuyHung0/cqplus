//! Implementation of cq+
//!
//! The name of functions follow the notation in Figure 1 on page 19 of the paper.

// Turn of warning for dead_code
// #![allow(dead_code)]
// #![allow(unused_imports)]
// #![allow(unused_variables)]

pub mod data_structures;
pub mod derive;
pub mod error;
pub mod kzg;
pub mod polynomials;
pub mod prover;
pub mod srs;
pub mod table;
pub mod verifier;

#[cfg(test)]
mod tests {
    use crate::{data_structures, derive, prover, srs, table::Table, verifier};
    use ark_bn254::Bn254;
    use ark_ec::PairingEngine;
    use ark_std::{
        rand::{rngs::StdRng, Rng, RngCore},
        test_rng, UniformRand,
    };

    fn prepare<E: PairingEngine, R: RngCore>(
        big_n1: usize,
        big_n2: usize,
        big_n: usize,
        small_n: usize,
        rng: &mut R,
    ) -> (
        Table<E::Fr>,
        Vec<E::Fr>,
        Vec<E::G1Affine>,
        Vec<E::G2Affine>,
        E::Fr,
        data_structures::VerificationKey<E>,
        E::G2Affine,
        E::Fr,
    ) {
        // Generate random and srs
        let secret_s = E::Fr::rand(rng);
        let (srs1, srs2) = srs::unsafe_setup_from_tau::<E>(big_n1, big_n2, secret_s);

        let table_values: Vec<_> = (0..big_n).map(|_| E::Fr::rand(rng)).collect();
        let table_t = Table::new(&table_values).unwrap();

        // Create random index from the table
        let subvector_indices: Vec<usize> = (0..small_n)
            .map(|_| rng.gen_range(0..small_n - 1))
            .collect();

        // Create a subvector from the table
        let vector_f: Vec<E::Fr> = subvector_indices
            .iter()
            .map(|&i| table_t.values[i])
            .collect();

        let (ek, vk, commit_poly_t2, value_vartheta) =
            derive::derive::<E>(&srs1, &srs2, &table_t, big_n, small_n);

        (
            table_t,
            vector_f,
            srs1,
            srs2,
            secret_s,
            vk,
            commit_poly_t2,
            value_vartheta,
        )
    }
    #[test]
    fn test_correct() {
        let two: usize = 2;
        let big_n = two.pow(4); // 64
        let small_n = two.pow(2); // 8

        let big_n1 = big_n + 4;
        let big_n2 = big_n + 6;

        let mut rng = test_rng();

        let (table_t, vector_f, srs1, srs2, secret_s, vk, commit_poly_t2, value_vartheta) =
            prepare::<Bn254, StdRng>(big_n1, big_n2, big_n, small_n, &mut rng);

        let proof = prover::prove::<Bn254>(secret_s, big_n1, big_n2, &table_t, &vector_f).unwrap();

        let result = verifier::verify::<Bn254>(&vk, &proof, commit_poly_t2, value_vartheta);
        assert!(result.is_ok());
    }
}
