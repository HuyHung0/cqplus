//! The prover module contains the implementation of the Prove function.
use ark_ec::{AffineCurve, PairingEngine};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial,
};
use ark_std::rand::Rng;
use ark_std::test_rng;
use std::cmp::min;
use std::ops::Mul;

use crate::{
    data_structures::Proof,
    derive,
    error::Error,
    kzg::Kzg,
    polynomials, srs,
    table::Table,
};
use ark_ff::{FftField, Field, One, UniformRand, Zero};

/// Compute the vector $m=(m_1,...,m_N)$ where $m_i$ is the number of time that $t_i$ appears in the vector $f$.
/// TODO: Need to check carefully the index: should start with 0 or 1?
pub fn vector_m<F: FftField>(table_t: &Table<F>, vector_f: &Vec<F>) -> Result<Vec<F>, Error> {
    // Create zero vector of length N=size table_t
    let mut m = vec![F::zero(); table_t.size];

    for fi in vector_f {
        let index = table_t.value_index_mapping.get(&fi);
        let err_str = format!("{}", fi);
        let index = index.ok_or(Error::ValueNotInTable(err_str))?;
        m[*index] = m[*index] + F::one();
    }
    Ok(m)
}

/// some TODOs:
///
/// TODO: precompute the commitment to the Lagrange basis polynomials and vanishing polynomial.
///
/// TODO: Handle the case $t_i+\beta=0$ or $f_i+\beta=0$.
///
/// TODO: Check $N_1,N_2 \geq N + \max(b_F,1)-1$.
///
/// TODO: find replacement for a secure random number generator.
pub fn prove<E: PairingEngine>(
    secret_s: E::Fr,
    big_n1: usize,
    big_n2: usize,
    table_t: &Table<E::Fr>,
    vector_f: &Vec<E::Fr>,
) -> Result<Proof<E>, Error> {
    // Create srs
    let (srs1, srs2) = srs::unsafe_setup_from_s::<E>(big_n1, big_n2, secret_s);

    // get the values N and n
    let big_n = table_t.size;
    let small_n = vector_f.len();

    // Create the group K and H of size N and n respectively
    let set_k = GeneralEvaluationDomain::<E::Fr>::new(big_n).unwrap();
    let set_h = GeneralEvaluationDomain::<E::Fr>::new(small_n).unwrap();
    let set_k: Vec<E::Fr> = set_k.elements().collect();
    let set_h: Vec<E::Fr> = set_h.elements().collect();

    // Create the random $\rho_m$
    let mut rng = test_rng();
    let random_rho_m = E::Fr::rand(&mut rng);

    // Compute the Lagrange basis polynomials and the vanishing polynomial over K and H
    let poly_lagrange_k = polynomials::poly_lagrange_basis_all(&set_k);
    // let poly_lagrange_h = polynomials::poly_lagrange_basis_all(&set_h);
    let poly_vanish_k = polynomials::poly_vanish(&set_k);
    let poly_vanish_h = polynomials::poly_vanish(&set_h);

    // Compute the commitment of Lagrange basis polynomials and the vanishing polynomial
    let commit_lagrange_k: Vec<E::G1Affine> = poly_lagrange_k
        .iter()
        .map(|poly| Kzg::<E>::commit_g1(&srs1, poly).into())
        .collect();
    // let commit_lagrange_h: Vec<E::G1Affine> = poly_lagrange_h
    //     .iter()
    //     .map(|poly| Kzg::<E>::commit_g1(&srs1, poly).into())
    //     .collect();
    let commit_vanish_k = Kzg::<E>::commit_g1(&srs1, &poly_vanish_k).into();
    // let commit_vanish_h = Kzg::<E>::commit_g1(&srs1, &poly_vanish_h).into();

    // Compute the vector $m=(m_1,\ldot,m_N)$ where $m_i$ is the number of time that $t_i$ appears in the vector $f$.
    let m = vector_m(table_t, vector_f).unwrap();

    // Compute the commitment $[m(s)]_1$
    // First way: Compute the polynomial $m(X)$ and then commit to it.
    //let commit_m = Kzg::<E>::commit_g1(&srs1, &poly_m).into();
    // Compute the polynomial $m(X)$
    //let poly_m = polynomials::poly_m(&m, &set_k, rho_m);

    // Second way: Commit to the basis lagrange polynomials and the vanishing polynomial and take the linear combination of the commitments.
    let mut commit_poly_m = E::G1Affine::zero();
    for i in 0..big_n {
        commit_poly_m = commit_poly_m + commit_lagrange_k[i].mul(m[i]).into();
    }
    commit_poly_m = commit_poly_m + commit_vanish_k.mul(random_rho_m).into();

    // Compute the polynomial $s(X)$ and the commitment $[S(s)]_1$
    // First, Random $R_s$ and $\rho_s$
    let random_r_s = E::Fr::rand(&mut rng);
    let random_rho_s = E::Fr::rand(&mut rng);
    // Second, Compute the polynomial $S(X)$ and the commitment $[S(s)]_1$
    let poly_s_s = &DensePolynomial::from_coefficients_slice(&[E::Fr::zero(), E::Fr::one()])
        .mul(random_r_s)
        + &polynomials::poly_vanish(&set_k).mul(random_rho_s);
    let commit_poly_s_s = Kzg::<E>::commit_g1(&srs1, &poly_s_s).into();

    // Random beta
    let beta = E::Fr::rand(&mut rng);

    // random rho_a
    let random_rho_a = E::Fr::rand(&mut rng);

    // Compute A_j
    let mut a = vec![E::Fr::zero(); table_t.size];
    for i in 0..table_t.size {
        let mi = m[i];
        let t_i_plus_beta_inverse = (table_t.values[i] + beta).inverse().unwrap();
        let a_j = mi * t_i_plus_beta_inverse;
        a[i] = a_j;
    }

    // Compute polynomial A(X)
    let poly_a = polynomials::poly_m(&a, &set_k, random_rho_a);

    // Compute commitment $[A(s)]_1$.
    // There are two ways to compute the commitment $[A(s)]_1$.

    // First way is to use the function Kzg::<E>::commit_g1(&srs1, &poly_a).into();
    let commit_poly_a_second_way = Kzg::<E>::commit_g1(&srs1, &poly_a).into();
    // Second way is to use the linear combination of the commitments to the Lagrange basis polynomials and the vanishing polynomial.
    // The second way is used here.
    let mut commit_poly_a = E::G1Affine::zero();
    for i in 0..big_n {
        commit_poly_a = commit_poly_a + commit_lagrange_k[i].mul(a[i]).into();
    }
    commit_poly_a = commit_poly_a + commit_vanish_k.mul(random_rho_a).into();

    if &commit_poly_a != &commit_poly_a_second_way {
        return Err(Error::CommitmentError);
    }
    else {
        dbg!("commit_poly_a is equal to commit_poly_a_second_way");
    }

    // Create random polynomial rho_b(X) of degree less than or equal 1
    let random_rho_bx_0 = E::Fr::rand(&mut rng);
    let random_rho_bx_1 = E::Fr::rand(&mut rng);
    let poly_random_rho_bx =
        DensePolynomial::from_coefficients_slice(&[random_rho_bx_0, random_rho_bx_1]);
    // let evaluate_poly_rho_bx_at_s = poly_random_rho_bx.evaluate(&secret_s);

    // Compute B_j
    let mut b = vec![E::Fr::zero(); small_n];
    for i in 0..vector_f.len() - 1 {
        let fi = vector_f[i];
        let f_i_plus_beta_inverse = (fi + beta).inverse().unwrap();
        let b_j = f_i_plus_beta_inverse;
        b[i] = b_j;
        dbg!(&i);
    }

    // Compute polynomial B(X)
    let poly_b = polynomials::poly_f(&b, &set_h, &poly_random_rho_bx);

    // Compute commitment $[B(s)]_1$
    // let mut commit_poly_b = E::G1Affine::zero();
    // for i in 0..small_n {
    //     commit_poly_b = commit_poly_b + commit_lagrange_h[i].mul(b[i]).into();
    // }
    // commit_poly_b = commit_poly_b + commit_vanish_h.mul(evaluate_poly_rho_bx_at_s).into();
    let commit_poly_b = Kzg::<E>::commit_g1(&srs1, &poly_b).into();

    // Compute polynomial $Q_B(X) = (B(X)(F(X)+\beta)-1)/\nu_H(X)$
    // First, compute poly_f
    let degree_bf = rand::thread_rng().gen_range(1..min(big_n1, big_n2) - big_n + 1); // i.e big_n1,big_n2 >= big_n+degree_bf-1
    let poly_random_rho_fx: DensePolynomial<E::Fr> = polynomials::poly_random(degree_bf);
    let poly_f = polynomials::poly_f(vector_f, &set_h, &poly_random_rho_fx);
    // Compute commitment $[F(s)]_1$
    let commit_poly_f = Kzg::<E>::commit_g1(&srs1, &poly_f).into();

    // Second, compute Q_B(X)
    let poly_f_plus_beta = poly_f + DensePolynomial::from_coefficients_slice(&[beta]);
    let numerator =
        &poly_b * &poly_f_plus_beta + DensePolynomial::from_coefficients_slice(&[-E::Fr::one()]);
    let (quotient, _reminder) = DenseOrSparsePolynomial::divide_with_q_and_r(
        &numerator.into(),
        &poly_vanish_h.clone().into(),
    )
    .unwrap();
    let poly_q_b: DensePolynomial<E::Fr> = quotient;

    // Compute commitment $[Q_B(s)]_1$
    let commit_poly_q_b = Kzg::<E>::commit_g1(&srs1, &poly_q_b).into();

    // Create random $\gamma$ and $\eta$
    let gamma = E::Fr::rand(&mut rng);
    let eta = E::Fr::rand(&mut rng);

    // Compute $B_\gamma$
    let b_gamma = poly_b.evaluate(&gamma);
    let commit_value_b_gamma =
        Kzg::<E>::commit_g1(&srs1, &DensePolynomial::from_coefficients_slice(&[b_gamma])).into();

    // Compute $D(X)=B_\gamma.(F(X)+\beta)-1-Q_B(X)\nu_H(\gamma)$
    let evaluate_vanish_h_at_gamma = poly_vanish_h.evaluate(&gamma);
    let poly_d = &poly_f_plus_beta.mul(b_gamma)
        + &DensePolynomial::from_coefficients_slice(&[-E::Fr::one()])
        + poly_q_b.mul(evaluate_vanish_h_at_gamma);

    // Compute commitment $[D(s)]_1$
    let commit_poly_d = Kzg::<E>::commit_g1(&srs1, &poly_d).into();

    // Compute P(X)=[(B(X)-B(\gamma))+\eta D(X)] / (X-\gamma)$
    let enumerator =
        &poly_b - &DensePolynomial::from_coefficients_slice(&[b_gamma]) + poly_d.mul(eta);
    let denominator = DensePolynomial::from_coefficients_slice(&[-gamma, E::Fr::one()]);
    let (quotient, _reminder) = DenseOrSparsePolynomial::divide_with_q_and_r(
        &enumerator.into(),
        &denominator.clone().into(),
    )
    .unwrap();
    let poly_p: DensePolynomial<E::Fr> = quotient;

    // commit [s-gamma]_2
    let commit_s_minus_gamma = Kzg::<E>::commit_g2(&srs2, &denominator).into();

    // Compute commitment $[P(s)]_1$
    let commit_poly_p = Kzg::<E>::commit_g1(&srs1, &poly_p).into();

    // Compute commitment $[R_C^{*}(s)]_1$
    let (ek, _, _, vartheta) = derive::derive::<E>(&srs1, &srs2, &table_t, big_n, small_n);
    let vartheta_inv = vartheta.inverse().unwrap();

    // Compute the first part of the commitment $[R_C^{*}(s)]_1$
    let mut commit_poly_r_c_part1 = E::G1Affine::zero();
    for i in 0..big_n {
        if m[i] == E::Fr::zero() {
            continue;
        } else {
            commit_poly_r_c_part1 = commit_poly_r_c_part1 + ek.commit_poly_rk[i].mul(a[i]).into();
        }
    }

    // Compute the second part of the commitment $[R_C^{*}(s)]_1$
    let mut commit_poly_r_c_part2 = E::G1Affine::zero();
    for i in 0..small_n {
        commit_poly_r_c_part2 = commit_poly_r_c_part2 + ek.commit_poly_rk[i].mul(b[i]).into();
    }

    commit_poly_r_c_part2 = commit_poly_r_c_part2.mul(-vartheta).into();

    // Compute the third part of the commitment $[R_C^{*}(s)]_1$
    let commit_poly_r_c_part3 = ek.commit_poly_u.mul(eta * random_r_s).into();

    // Finally, compute the commitment $[R_C^{*}(s)]_1$
    let commit_poly_r_c = commit_poly_r_c_part1 + commit_poly_r_c_part2 + commit_poly_r_c_part3;

    // Compute the commitment $[Q_A(s)]_1$
    // Compute the first part of the commitment $[Q_A(s)]_1$
    let mut commit_poly_q_a_part1 = E::G1Affine::zero();
    for i in 0..big_n {
        if m[i] == E::Fr::zero() {
            continue;
        } else {
            commit_poly_q_a_part1 = commit_poly_q_a_part1 + ek.commit_poly_q_j[i].mul(a[i]).into();
        }
    }

    // Compute the second part of the commitment $[Q_A(s)]_1$
    let poly_part2 = &DensePolynomial::from_coefficients_slice(&[random_rho_a])
        * &(polynomials::poly_t(&table_t, &set_k)
            + DensePolynomial::from_coefficients_slice(&[beta]))
        + DensePolynomial::from_coefficients_slice(&[random_rho_m]);
    let commit_poly_q_a_part2 = Kzg::<E>::commit_g1(&srs1, &poly_part2).into();

    // Finally, compute the commitment $[Q_A(s)]_1$
    let commit_poly_q_a = commit_poly_q_a_part1 + commit_poly_q_a_part2;

    // Compute the commitment $[Q_C(s)]_1$
    let poly = &DensePolynomial::from_coefficients_slice(&[random_rho_a])
        + &(&DensePolynomial::from_coefficients_slice(&[vartheta_inv]) * &poly_random_rho_bx);
    let commit_poly_q_c = Kzg::<E>::commit_g1(&srs1, &poly).into();

    // Compute the commitment $[Q(s)]_1$
    // Compute the third part of the commitment $[Q(s)]_1$
    let commit_poly_q_part3 = Kzg::<E>::commit_g1(
        &srs1,
        &DensePolynomial::from_coefficients_slice(&[random_rho_s]),
    )
    .into();

    // Finally, compute the commitment $[Q_C(s)]_1$
    let commit_poly_q = commit_poly_q_a
        + commit_poly_q_c.mul(eta).into()
        + commit_poly_q_part3.mul(-eta * eta).into();

    let poly_x = DensePolynomial::from_coefficients_slice(&[E::Fr::zero(), E::Fr::one()]);
    // commit to poly_x
    let commit_poly_x = Kzg::<E>::commit_g2(&srs2, &poly_x).into();


    Ok(Proof {
        commit_poly_m,
        commit_poly_s: commit_poly_s_s,
        commit_poly_a,
        commit_poly_b,
        commit_poly_q_b,
        commit_poly_p,
        commit_poly_r_c,
        commit_poly_q,
        commit_poly_d,
        commit_poly_f,
        commit_poly_x,
        commit_s_minus_gamma,
        commit_value_b_gamma,
        value_b_gamma: b_gamma,
        value_beta: beta,
        value_gamma: gamma,
        value_eta: eta,
    })
}
