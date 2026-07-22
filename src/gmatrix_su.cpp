/*
 AGHmatrix: Relationship Matrices for Diploid and Autopolyploid Species
 Copyright (C) 2016-2026 Rodrigo Amadeu

 This file is part of AGHmatrix.

AGHmatrix is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 AGHmatrix is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

/*
 File: gmatrix_su.cpp
 Description: implement gmatrix su method

 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// gmatrix_su.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Diploid alternative allele frequency p_j = (2*n(2) + n(1)) / (2*nobs)
static inline double diploid_alt_freq(const vec& col) {
  double c2 = 0.0, c1 = 0.0;
  uword nobs = 0;
  const uword n = col.n_elem;
  for (uword i = 0; i < n; ++i) {
    const double v = col[i];
    if (!std::isnan(v)) {
      if      (v == 2.0) c2 += 1.0;
      else if (v == 1.0) c1 += 1.0;
      // v==0 (or fractional) contributes to denominator only
      nobs += 1;
    }
  }
  if (nobs == 0) return NA_REAL;
  return (2.0 * c2 + c1) / (2.0 * static_cast<double>(nobs));
}

// [[Rcpp::export]]
arma::mat Gmatrix_Su(const arma::mat& M) {
  const uword n = M.n_rows, m = M.n_cols;
  
  // 1) p_j and TwoPQ_j per marker (p is alternative-allele frequency)
  vec p(m, fill::zeros);
  for (uword j = 0; j < m; ++j) p[j] = diploid_alt_freq(M.col(j));
  vec TwoPQ = 2.0 * p % (1.0 - p);
  
  // 2) Denominator: sum_j TwoPQ_j * (1 - TwoPQ_j)    (legacy)
  const double denom = accu(TwoPQ % (1.0 - TwoPQ));
  
  //    0/2 -> 0; 1 stays 1; fractional imputed values pass through unchanged
  mat Z(n, m, fill::zeros);
  for (uword j = 0; j < m; ++j) {
    const double tpq = TwoPQ[j];
    for (uword i = 0; i < n; ++i) {
      const double v = M(i, j);
      if (!std::isnan(v)) {
        const double h = (v == 0.0 || v == 2.0) ? 0.0 : v;
        Z(i, j) = h - tpq;
      }
    }
  }
  
  // 4) G = Z Z' / denom
  return (Z * Z.t()) / denom;
}
