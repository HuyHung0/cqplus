//! Structures of EvaluationKey, VerificationKey, and Proof
use ark_ec::PairingEngine;

/// A list of elements in G1 which are commitments of some values.
pub struct EvaluationKey<E: PairingEngine> {
    pub commit_poly_rk: Vec<E::G1Affine>,  // commitment of r_K(s)
    pub commit_poly_rh: Vec<E::G1Affine>,  // commitment of r_H(s)
    pub commit_poly_u: E::G1Affine,        // commitment of U(s)
    pub commit_poly_v: E::G1Affine,        // commitment of v_K(s)
    pub commit_s_poly_v: E::G1Affine,      // commitment of s*v_K(s)
    pub commit_poly_q_j: Vec<E::G1Affine>, // commitment of Q_j(s)
    pub commit_poly_t: E::G1Affine,        // commitment of T(s)
}

/// A list of elements in G2 which are commitments of some values.
pub struct VerificationKey<E: PairingEngine> {
    pub commit_1: E::G2Affine,        // commitment of 1 in G2
    pub commit_poly_u2: E::G2Affine,  // commitment of U(s)
    pub commit_poly_zu2: E::G2Affine, // commitment of z(s)U(s)
    pub commit_poly_vu2: E::G2Affine, // commitment of v_K(s)U(s)
    pub commit_poly_tu2: E::G2Affine, // commitment of T(s)U(s)
}

/// The proof structure is the return of Prove() function.

pub struct Proof<E: PairingEngine> {
    pub commit_poly_m: E::G1Affine,        // commitment of m(s)
    pub commit_poly_s: E::G1Affine,        // commitment of S(s)
    pub commit_poly_a: E::G1Affine,        // commitment of A(s)
    pub commit_poly_b: E::G1Affine,        // commitment of B(s)
    pub commit_poly_q_b: E::G1Affine,      // commitment of Q_B(s)
    pub commit_poly_p: E::G1Affine,        // commitment of P(s)
    pub commit_poly_r_c: E::G1Affine,      // commitment of R_C^{*}(s)
    pub commit_poly_q: E::G1Affine,        // commitment of Q(s)
    pub commit_poly_d: E::G1Affine,        // commitment of D(s)
    pub commit_poly_f: E::G1Affine,        // commitment of F(s)
    pub commit_poly_x: E::G2Affine,        // commitment of s in G2
    pub commit_s_minus_gamma: E::G2Affine, // commitment of [s - gamma]_2
    pub commit_value_b_gamma: E::G1Affine, // commitment of B_gamma = B(gamma
    pub value_b_gamma: E::Fr,              // value of B(s) at gamma
    pub value_beta: E::Fr,                 // value of random beta
    pub value_gamma: E::Fr,                // value of random gamma
    pub value_eta: E::Fr,                  // value of random eta
}
