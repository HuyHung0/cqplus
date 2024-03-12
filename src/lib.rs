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
    use crate::{derive, prover, srs, table::Table, verifier};
    use ark_bn254::Bn254;
    use ark_ec::PairingEngine;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::{
        rand::{rngs::StdRng, Rng, RngCore},
        test_rng, UniformRand,
    };

    fn prepare1<E: PairingEngine, R: RngCore>(
        big_n1: usize,
        big_n2: usize,
        big_n: usize,
        small_n: usize,
        rng: &mut R,
    ) -> (
        Vec<E::G1Affine>, // srs1
        Vec<E::G2Affine>, // srs2
        GeneralEvaluationDomain<E::Fr>, //set_k
        GeneralEvaluationDomain<E::Fr>, //set_h
    ) {
        // Create the group K and H of size N and n respectively
        let set_k = GeneralEvaluationDomain::<E::Fr>::new(big_n).unwrap();
        let set_h = GeneralEvaluationDomain::<E::Fr>::new(small_n).unwrap();
        
        // Generate srs1 and srs2 from rng
        let (srs1, srs2) = srs::unsafe_setup_from_rng::<E,R>(big_n1, big_n2, rng);
        
        (
            srs1,
            srs2,
            set_k,
            set_h,
        )
    }

    fn prepare2<E: PairingEngine, R: RngCore>(
        big_n: usize,
        small_n: usize,
        rng: &mut R,
    ) -> (
        Table<E::Fr>, // Table t
        Vec<E::Fr>, // Vector f
    ) {

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

        (
            table_t,
            vector_f,
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

        // Prepare the srs's, groups
        let (srs1, srs2, set_k, set_h) = prepare1::<Bn254, StdRng>(big_n1, big_n2, big_n, small_n, &mut rng);

        // Prepare the table and vector
        let (table_t, vector_f) = prepare2::<Bn254, StdRng>(big_n, small_n, &mut rng);
        
        let (ek, vk, commit_poly_t2, vartheta) =
            derive::derive::<Bn254>(&srs1, &srs2, &table_t, &set_k, &set_h);
        

        let proof = prover::prove::<Bn254>(srs1,srs2,&table_t,&vector_f,&set_k,&set_h,&ek, vartheta).unwrap();

        let result = verifier::verify::<Bn254>(&vk, &proof, commit_poly_t2, vartheta);
        assert!(result.is_ok());
    }
}
