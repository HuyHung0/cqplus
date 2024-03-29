//! Compute Kzg commitments of a polynomial in G1, G2. Taken from the cq crate.

use std::marker::PhantomData;

use ark_ec::{msm::VariableBaseMSM, PairingEngine};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Polynomial};


/// Provide function `commit_g1` and `commit_g2` to commit to a polynomial in G1, G2 using KZG. If $$P(X)=p_0 + p_1 X +\ldots + p_n X^n$$ is a polynomial, and $$srs1 =\\{ g^{s^i} \\}\_{i=0}^n = \\{[s^i]\\}_{i=0}^{n}$$ is a srs in G1, then
/// $$[P(s)]_1= g^{P(s)}= g^{p_0 + p_1 s +\ldots + p_n s^n} = \sum\_{i=1}^n p_i[s^i].$$
/// We use the struct `sms::VariableBaseMSM` from crate `ark-ec` to compute the multiplication of group elements by a scalar more efficiently.
pub struct Kzg<E: PairingEngine> {
    _e: PhantomData<E>,
}

impl<E: PairingEngine> Kzg<E> {
    pub fn commit_g1(srs: &[E::G1Affine], poly: &DensePolynomial<E::Fr>) -> E::G1Projective {
        if srs.len() - 1 < poly.degree() {
            // Max power of srs should be equal or greater than the degree of the polynomial.
            // Index of power of srs runs from 0 to Max power (len = max power +1)
            panic!(
                "SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                srs.len()
            );
        }
        let coeff_scalars: Vec<_> = poly.coeffs.iter().map(|c| c.into_repr()).collect();

        // The crate ark_ec provides module `msm` and struct `VariableBaseMSM` to compute the multiplication of group elements by a scalar more efficiently.
        VariableBaseMSM::multi_scalar_mul(srs, &coeff_scalars)
    }

    pub fn commit_g2(srs: &[E::G2Affine], poly: &DensePolynomial<E::Fr>) -> E::G2Projective {
        if srs.len() - 1 < poly.degree() {
            // Max power of srs should be equal or greater than the degree of the polynomial.
            // Index of power of srs runs from 0 to Max power (len = max power +1)
            panic!(
                "SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                srs.len()
            );
        }
        let coeff_scalars: Vec<_> = poly.coeffs.iter().map(|c| c.into_repr()).collect();
        VariableBaseMSM::multi_scalar_mul(srs, &coeff_scalars)
    }

    // pub fn open_g1(
    //     srs: &[E::G1Affine],
    //     poly: &DensePolynomial<E::Fr>,
    //     challenge: E::Fr,
    // ) -> (E::Fr, E::G1Affine) {
    //     let q = poly / &DensePolynomial::from_coefficients_slice(&[-challenge, E::Fr::one()]);
    //     if srs.len() - 1 < q.degree() {
    //         panic!(
    //             "Open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
    //             q.degree(),
    //             srs.len()
    //         );
    //     }
    //     let proof = Self::commit_g1(srs, &q);
    //     (poly.evaluate(&challenge), proof.into())
    // }

    // pub fn batch_open_g1(
    //     srs: &[E::G1Affine],
    //     polys: &[DensePolynomial<E::Fr>],
    //     opening_challenge: E::Fr,
    //     separation_challenge: E::Fr,
    // ) -> E::G1Affine {
    //     let powers_of_gamma = iter::successors(Some(separation_challenge), |p| {
    //         Some(*p * separation_challenge)
    //     });

    //     let mut batched = polys[0].clone();
    //     for (p_i, gamma_pow_i) in polys.iter().skip(1).zip(powers_of_gamma) {
    //         batched += (gamma_pow_i, p_i);
    //     }

    //     let q = &batched
    //         / &DensePolynomial::from_coefficients_slice(&[-opening_challenge, E::Fr::one()]);

    //     if srs.len() - 1 < q.degree() {
    //         panic!(
    //             "Batch open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
    //             q.degree(),
    //             srs.len()
    //         );
    //     }

    //     Self::commit_g1(srs, &q).into()
    // }
}
