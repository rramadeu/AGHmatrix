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
 File: diploid_p0p2_TwoPQ.cpp
 Description: implement counting genotypes by column

 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// Count diploid genotypes by column (ignore NA), return p0, p2, TwoPQ
// [[Rcpp::export]]
Rcpp::List diploid_p0p2_TwoPQ_cpp(const arma::mat& M) {
  const uword n = M.n_rows, m = M.n_cols;
  arma::vec p0(m, fill::zeros), p2(m, fill::zeros);

  for (uword j = 0; j < m; ++j) {
    double c0 = 0.0, c1 = 0.0, c2 = 0.0, nobs = 0.0;
    for (uword i = 0; i < n; ++i) {
      const double v = M(i,j);
      if (!std::isnan(v)) {
        nobs += 1.0;
        if      (v == 0.0) c0 += 1.0;
        else if (v == 1.0) c1 += 1.0;
        else if (v == 2.0) c2 += 1.0;
      }
    }
    const double d = 2.0 * nobs;
    if (d > 0.0) {
      p0(j) = (2.0 * c0 + c1) / d;
      p2(j) = (2.0 * c2 + c1) / d;
    } else {
      p0(j) = 0.0;
      p2(j) = 0.0;
    }
  }
  const double TwoPQ = 2.0 * arma::dot(p0, p2);
  return Rcpp::List::create(
    Rcpp::Named("p0")    = p0,
    Rcpp::Named("p2")    = p2,
    Rcpp::Named("TwoPQ") = TwoPQ
  );
}
