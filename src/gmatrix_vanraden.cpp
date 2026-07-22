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
 File: gmatrix_vanraden.cpp
 Description: implement gmatrix vanraden method

 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_vanraden(const arma::mat& M,
                           const arma::vec& FreqP,
                           const double denom) {
  const uword n = M.n_rows, m = M.n_cols;
  if (FreqP.n_elem != m) stop("FreqP length must equal number of markers");
  if (!(denom > 0.0)) stop("Denominator (TwoPQ) must be > 0");
  
  mat Z = M;
  for (uword j = 0; j < m; ++j) {
    const double mu2p = 2.0 * FreqP[j];
    for (uword i = 0; i < n; ++i) {
      double &z = Z(i,j);
      if (!std::isnan(z)) z -= mu2p; else z = 0.0;
    }
  }
  return (Z * Z.t()) / denom;
}

// [[Rcpp::export]]
arma::mat Gmatrix_VanRaden_weighted(const arma::mat& M,
                                    const arma::vec& weights,
                                    const arma::vec& FreqP,
                                    const double denom) {
  const uword n = M.n_rows, m = M.n_cols;
  if (weights.n_elem != m) stop("weights length must equal number of markers");
  if (FreqP.n_elem   != m) stop("FreqP length must equal number of markers");
  if (!(denom > 0.0)) stop("Denominator (TwoPQ) must be > 0");
  
  mat Z = M;
  for (uword j = 0; j < m; ++j) {
    const double mu2p = 2.0 * FreqP[j];
    for (uword i = 0; i < n; ++i) {
      double &z = Z(i,j);
      if (!std::isnan(z)) z -= mu2p; else z = 0.0;
    }
  }
  
  mat W = diagmat(weights);
  return (Z * W * Z.t()) / denom;
}
