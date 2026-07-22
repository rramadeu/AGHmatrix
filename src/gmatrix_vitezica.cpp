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
 File: gmatrix_vitezica.cpp
 Description: implement gmatrix vitezica method

 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_Vitezica(arma::mat M, const arma::mat& FreqP) {
  const uword n = M.n_rows, m = M.n_cols;
  if (FreqP.n_rows != n || FreqP.n_cols != m)
    stop("FreqP must be n x m matrix (replicated column-wise frequencies).");
  
  rowvec F = FreqP.row(0); // vector of FreqP per marker
  vec TwoPQ = 2.0 * F.t() % (1.0 - F.t());
  const double denom = accu(square(TwoPQ));
  
  // replace NA by 0 BEFORE recoding
  na_to_zero(M);
  
  mat X(n, m, fill::zeros);
  for (uword j = 0; j < m; ++j) {
    const double p = F[j];
    const double term0 = -2.0 * p * p;
    const double term1 =  2.0 * p * (1.0 - p);
    const double term2 = -2.0 * (1.0 - p) * (1.0 - p);
    for (uword i = 0; i < n; ++i) {
      const double v = M(i,j);
      if (v == 0.0)      X(i,j) = term0;
      else if (v == 1.0) X(i,j) = term1;
      else if (v == 2.0) X(i,j) = term2;
      else               X(i,j) = 0.0; // safety for unusual codes
    }
  }
  
  return (X * X.t()) / denom;
}
