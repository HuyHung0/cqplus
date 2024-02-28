//! Return the Evaluation key and Verification key

use ark_ec::PairingEngine;
use ark_ff::{FftField, One, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
//use ark_ff::{One,Zero};

// use ark_bn254::Fr;

use crate::{
    data_structures::{EvaluationKey, VerificationKey},kzg::Kzg, polynomials, table::Table
};

/// Compute poly_vanish_k_minus_h, $\nu_{K\H}(X)$
pub fn poly_vanish_k_minus_h<F: FftField>(set_k: &[F], set_h: &[F]) -> DensePolynomial<F> {
    let poly_vanish_k: DensePolynomial<F> = polynomials::poly_u(set_k.len());
    let poly_vanish_h: DensePolynomial<F> = polynomials::poly_u(set_h.len());
    let (quotient, reminder) = DenseOrSparsePolynomial::divide_with_q_and_r(
        &poly_vanish_k.into(),
        &poly_vanish_h.into(),
    ).unwrap();
    if reminder != DensePolynomial::zero() {
        panic!("The remainder of the division is not zero");
    }
    quotient
}

/// Compute $r_j^K(X)$, where $j = 1,..,N$ and $N=|K|$.
pub fn poly_rk<F: FftField>(set: &[F], big_n1: usize, big_n: usize) -> Vec<DensePolynomial<F>> {
    // Declare some basic variables and polynomials
    let mu = big_n1 - big_n + 2;
    let u: DensePolynomial<F> = polynomials::poly_u(mu);
    let list_lagrange_polys = polynomials::poly_lagrange_basis_all(set);

    // Create a vector to store the result
    let mut r = Vec::with_capacity(set.len());

    for (i, _) in set.iter().enumerate() {
        // l_i is the i-th Lagrange basis polynomial over the domain K
        let l_i = &list_lagrange_polys[i];

        let r_i = &DensePolynomial::from_coefficients_slice(&l_i.coeffs[1..]) * &u;

        r.push(r_i);
    }
    r
}

/// Compute $r_j^H(X)$, where $j=1,...,n$ and $n=|H|$
/// set_k = K, set_h = H
pub fn poly_rh<F: FftField>(
    set_k: &[F],
    set_h: &[F],
    big_n1: usize,
    big_n: usize,
    _small_n: usize, //not used
) -> Vec<DensePolynomial<F>> {
    let mu = big_n1 - big_n + 2;
    let mut r = Vec::with_capacity(set_h.len());

    let list_lagrange_polys = polynomials::poly_lagrange_basis_all(set_h);

    let u: DensePolynomial<F> = polynomials::poly_u(mu);

    // Find the difference between the two sets and compute the vanishing polynomial on the difference set

    let vanish_k_minus_h = poly_vanish_k_minus_h(set_k, set_h);

    for (i, _) in set_h.iter().enumerate() {
        let l_i = &list_lagrange_polys[i];
        // l_i is the i-th Lagrange basis polynomial over the domain H

        let poly1 = l_i * &vanish_k_minus_h;
        let r_i = &DensePolynomial::from_coefficients_slice(&poly1.coeffs[1..]) * &u;

        r.push(r_i);
    }
    r
}

/// Compute $Q_j(X)$, where $j=1,...,N$ and $n=|K|$
pub fn poly_q_j<F: FftField>(
    table: &Table<F>,
    set_k: &[F],
    _set_h: &[F], //not used
    _big_n1: usize, //not used
    _big_n: usize, //not used
    _small_n: usize, //not used
) -> Vec<DensePolynomial<F>> {

    let mut r = Vec::with_capacity(set_k.len());

    let list_lagrange_polys = polynomials::poly_lagrange_basis_all(set_k);
    let t_x = polynomials::poly_t(table, set_k);

    let vanish_k: DensePolynomial<F> = polynomials::poly_u(set_k.len());
    


    for (i, _) in set_k.iter().enumerate() {
        let t_i = table.values[i];
        let poly_t_i = DensePolynomial::from_coefficients_vec(vec![-t_i]);
        let l_i = &list_lagrange_polys[i];
        // l_i is the i-th Lagrange basis polynomial over the domain K

        let numerator = &(t_x.clone() + poly_t_i) * l_i;
        let (quotient, reminder) = DenseOrSparsePolynomial::divide_with_q_and_r(
            &numerator.into(),
            &vanish_k.clone().into(),
        )
        .unwrap();
    if reminder != DensePolynomial::zero() {
        panic!("The remainder of the division is not zero");
    }
        r.push(quotient);
    }
    r
}

/// Return the Evaluation key and Verification key
pub fn derive<E: PairingEngine>(
    srs1: &Vec<E::G1Affine>,
    srs2: &Vec<E::G2Affine>,
    table: &Table<E::Fr>,
    big_n: usize,
    small_n: usize,
) -> (EvaluationKey<E>, VerificationKey<E>, E::G2Affine, E::Fr) {

    let big_n1 = srs1.len()-1;
    // let big_n2 = srs2.len()-1;

    //create two groups of elements
    let set_k = GeneralEvaluationDomain::<E::Fr>::new(big_n).unwrap();
    let set_h = GeneralEvaluationDomain::<E::Fr>::new(small_n).unwrap();

    //convert the elements to vectors
    let set_k: Vec<E::Fr> = set_k.elements().collect();
    let set_h: Vec<E::Fr> = set_h.elements().collect();

    let poly_rk = poly_rk(&set_k, big_n1, big_n);
    let poly_rh = poly_rh(&set_k, &set_h, big_n1, big_n, small_n);
    let poly_q_j = poly_q_j(&table, &set_k, &set_h, big_n1, big_n, small_n);

    let commit_poly_rk: Vec<E::G1Affine> = poly_rk
        .iter()
        .map(|poly| Kzg::<E>::commit_g1(&srs1, poly).into())
        .collect();
    let commit_poly_rh: Vec<E::G1Affine> = poly_rh
        .iter()
        .map(|poly| Kzg::<E>::commit_g1(&srs1, poly).into())
        .collect();
    let commit_poly_q_j: Vec<E::G1Affine> = poly_q_j
        .iter()
        .map(|poly| Kzg::<E>::commit_g1(&srs1, poly).into())
        .collect();

    let poly_u = polynomials::poly_u(big_n1 - big_n + 2);
    let commit_poly_u = Kzg::<E>::commit_g1(&srs1, &poly_u).into();
    let commit_poly_u2 = Kzg::<E>::commit_g2(&srs2, &poly_u).into();


    let poly_v = polynomials::poly_u(set_k.len());
    let commit_poly_v = Kzg::<E>::commit_g1(&srs1, &poly_v).into();

    let s_poly_v = polynomials::poly_x_vanish(&set_k);
    let commit_s_poly_v = Kzg::<E>::commit_g1(&srs1, &s_poly_v).into();

    let poly_t = polynomials::poly_t(&table, &set_k);
    let commit_poly_t = Kzg::<E>::commit_g1(&srs1, &poly_t).into();
    let commit_poly_t2 = Kzg::<E>::commit_g2(&srs2, &poly_t).into();


    // Verification Key
    let poly_1 = polynomials::poly_1();
    let commit_1 = Kzg::<E>::commit_g2(&srs2, &poly_1).into();



    let poly_vanish_k_minus_h = poly_vanish_k_minus_h(&set_k, &set_h);
    let poly_zu = &poly_vanish_k_minus_h * &poly_u;
    let commit_poly_zu2 = Kzg::<E>::commit_g2(&srs2, &poly_zu).into();

   

    let poly_vu = &poly_v * &poly_u;
    let commit_poly_vu2 = Kzg::<E>::commit_g2(&srs2, &poly_vu).into();

    
    let poly_tu = &poly_t * &poly_u;
    let commit_poly_tu2 = Kzg::<E>::commit_g2(&srs2, &poly_tu).into();



    let v1 = big_n / small_n;
    let mut value_vartheta = E::Fr::zero();
    for _ in 0..v1 {
        value_vartheta = value_vartheta + E::Fr::one();
    }

    (
        EvaluationKey {
            commit_poly_rk,
            commit_poly_rh,
            commit_poly_u,
            commit_poly_v,
            commit_s_poly_v,
            commit_poly_q_j,
            commit_poly_t,
        },
        VerificationKey {
            commit_1,
            commit_poly_u2,
            commit_poly_zu2,
            commit_poly_vu2,
            commit_poly_tu2,
        },
        commit_poly_t2,
        value_vartheta,
    )
}

#[cfg(test)]
mod test_derive {}
