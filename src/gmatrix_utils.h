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
 File: gmatrix_utils.cpp
 Description: implement gmatrix utils

 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// gmatrix_utils.h
//
// [[Rcpp::depends(RcppArmadillo)]]
#ifndef GMATRIX_UTILS_H
#define GMATRIX_UTILS_H

#include <RcppArmadillo.h>
#include <armadillo>
using namespace arma;

// freq of 2-allele (diploid 0/1/2) with NA ignored
inline double p2_allele_freq(const colvec& x) {
  double count = 0.0; int nobs = 0;
  for (uword i = 0; i < x.n_rows; ++i) {
    double v = x[i];
    if (!std::isnan(v)) {
      if (v == 2.0) count += 2.0;
      else if (v == 1.0) count += 1.0;
      nobs++;
    }
  }
  return (nobs > 0) ? (count / (2.0 * nobs)) : 0.0;
}

// column means with NA removed
inline arma::rowvec colMeans_na_rm(const arma::mat& X) {
  arma::rowvec mu(X.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < X.n_cols; ++j) {
    double s = 0.0; int nobs = 0;
    for (arma::uword i = 0; i < X.n_rows; ++i) {
      double v = X(i,j);
      if (!std::isnan(v)) { s += v; nobs++; }
    }
    mu[j] = (nobs > 0) ? (s / nobs) : arma::datum::nan;
  }
  return mu;
}

// column sample variances (n-1 denom), NA removed; computed on centred X
inline arma::rowvec colVars_na_rm_centered(const arma::mat& Xc) {
  arma::rowvec v(Xc.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < Xc.n_cols; ++j) {
    double ss = 0.0; int nobs = 0;
    for (arma::uword i = 0; i < Xc.n_rows; ++i) {
      double z = Xc(i,j);
      if (!std::isnan(z)) { ss += z*z; nobs++; }
    }
    v[j] = (nobs > 1) ? (ss / (nobs - 1.0)) : 0.0;
  }
  return v;
}

// replace NaN by 0 in-place
inline void na_to_zero(arma::mat& X) {
  for (arma::uword i = 0; i < X.n_rows; ++i) {
    for (arma::uword j = 0; j < X.n_cols; ++j) {
      if (std::isnan(X(i, j))) {
        X(i, j) = 0.0;
      }
    }
  }
}

#endif // GMATRIX_UTILS_H
