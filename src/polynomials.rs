//! Contains some polynomials and functions using for implementation. It is similar to the derive function describe in the paper.

use ark_ff::FftField;

use crate::table::Table;

use ark_poly::{
    univariate::DensePolynomial,
    UVPolynomial,
};
// DensePolynomial: polynomial written in the form of its coefficients including zero coefficient.

/// Vanishing polynomial on a set
pub fn poly_vanish<F: FftField>(set: &[F]) -> DensePolynomial<F> {
    // polynomial = 1
    let mut polynomial = DensePolynomial::from_coefficients_slice(&[F::one()]);

    // Iterate over the set
    for (i, _) in set.iter().enumerate() {
        let x_i = set[i];

        // x_minus_xi = X-xi
        let x_minus_xi = DensePolynomial::from_coefficients_slice(&[-x_i, F::one()]);

        polynomial = &polynomial * &x_minus_xi;
    }

    // Return the vanishing polynomial
    polynomial
}

/// Vanishing polynomial multiplies by X
pub fn poly_x_vanish<F: FftField>(set: &[F]) -> DensePolynomial<F> {
    let poly = poly_vanish(set);
    &poly * &DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()])
}

/// poly 1
pub fn poly_1<F: FftField>() -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_slice(&[F::one()])
}

/// Return a vector containing all the Lagrange basis polynomials over a set. Index from 0 to #set - 1.
pub fn poly_lagrange_basis_all<F: FftField>(set: &[F]) -> Vec<DensePolynomial<F>> {
    // create a new empty vector with pre-allocated capacity to store set.len() elements
    // This vector will be used to store all the Lagrange basis polynomials
    let mut list_lagrange_polys = Vec::with_capacity(set.len());

    // index from 0 to set.len()-1
    for i in 0..set.len() {
        // l_i = 1
        let mut l_i = DensePolynomial::from_coefficients_slice(&[F::one()]);
        // DensePolynomial::from_coefficients_slice: Create a new polynomial from a list of coefficients.
        // the coefficient of $x^i$ is stored at location i in self.coeffs

        // x_i is the i-th element of the set
        let x_i = set[i];

        // j iterates through the set
        for (j, _) in set.iter().enumerate() {
            if j != i {
                // xi_minus_xj_inv = 1/(xi - xj)
                let xi_minus_xj_inv = (x_i - set[j]).inverse().unwrap();

                l_i = &l_i
                    * &DensePolynomial::from_coefficients_slice(&[
                        -set[j] * xi_minus_xj_inv,
                        xi_minus_xj_inv,
                    ]);
                // Coeff of x^0: -xj/(xi-xj)
                // Coeff of x^1: 1/(xi-xj)
                // l_i = l_i * [- xj/(xi-xj) + 1/(xi-xj) * X] = l_i * (X-xj)/(xi-xj)
            }
        }
        // Move the Lagrange basis polynomial l_i into the vector list_lagrange_polys at the end.
        list_lagrange_polys.push(l_i);
    }

    // Return the vector containing all Lagrange basis polynomials
    list_lagrange_polys
}

/// If a set is a subgroup of size n in a finite field, the vanishing polynomial on this set is poly_u = $U(X) = X^n-1$
pub fn poly_u<F: FftField>(n: usize) -> DensePolynomial<F> {
    let mut coefficients = vec![F::zero(); n];
    coefficients[0] = -F::one();
    coefficients[n - 1] = F::one();
    DensePolynomial::from_coefficients_vec(coefficients)
}

/// poly_t = $T(X)$, is the encoded polynomial of the table $t = ( t_j )_{j=1}^{N}$ over a group K with $|K|=N$.
///
/// It is defined as $T(X)=\sum_{i=1}^N t_j \lambda_j^K(X)$, where $\lambda_j^K(X)$ is the $j$-th Lagrange basis polynomial over $K$.
/// There is another way to compute $T(X)$: inverse fast fourier transform (IFFT) of the table $t$ over $K$.
pub fn poly_t<F: FftField>(table: &Table<F>, set: &[F]) -> DensePolynomial<F> {
    // Create a zero polynomial
    let mut t = DensePolynomial::from_coefficients_slice(&[F::zero()]);

    // list_lagrange_polys contains all the Lagrange basis polynomials over the domain K
    let list_lagrange_polys = poly_lagrange_basis_all(set);

    // compute T(x)
    for (i, _) in set.iter().enumerate() {
        t = &t + &(&list_lagrange_polys[i] * table.values[i]);
    }
    t
}

/// poly_f = $f(X)=\sum_{i=1}^n f_i \lambda_i^H(X) + \rho_F(X).\nu_H(X)$, where $\lambda_i^H(X)$ is the $i$-th Lagrange basis polynomial over $H$, $\rho_F(X)$ is a random polynomial.
pub fn poly_f<F: FftField>(
    vector_f: &Vec<F>,
    set: &[F],
    rho: &DensePolynomial<F>,
) -> DensePolynomial<F> {
    let mut f = DensePolynomial::from_coefficients_slice(&[F::zero()]);
    let list_lagrange_polys = poly_lagrange_basis_all(set);
    for (i, _) in set.iter().enumerate() {
        f = &f + &(&list_lagrange_polys[i] * vector_f[i]);
    }
    let result = f + rho * &poly_vanish(set);
    result
}

/// poly_random = $\rho_F(X)$, a random polynomial of predefine degree ($< b_F$)
/// TODO: make sure the coefficients are random (using some cryptographic secure random functions)
pub fn poly_random<F: FftField>(degree: usize) -> DensePolynomial<F> {
    let mut rng = ark_std::test_rng();
    let mut coefficients = Vec::with_capacity(degree + 1);
    for _ in 0..degree {
        coefficients.push(F::rand(&mut rng));
    }
    DensePolynomial::from_coefficients_vec(coefficients)
}

/// poly_m = $m(X)=\sum_{i=1}^N m_i \lambda_i^K(X) + \rho_m \vu_K(X)$,
/// where $\lambda_i^K(X)$ is the $i$-th Lagrange basis polynomial over $K$,
/// $\rho_m$ is a random element in the finite field,
/// $\vu_K(X)$ is the vanishing polynomial on the set $K$.
pub fn poly_m<F: FftField>(vector_m: &Vec<F>, set_k: &[F], rho_m: F) -> DensePolynomial<F> {
    let mut f = DensePolynomial::from_coefficients_slice(&[F::zero()]);
    let list_lagrange_polys = poly_lagrange_basis_all(set_k);
    for (i, _) in set_k.iter().enumerate() {
        f = &f + &(&list_lagrange_polys[i] * vector_m[i]);
    }
    let result = f + &DensePolynomial::from_coefficients_slice(&[rho_m]) * &poly_vanish(set_k);
    result
}