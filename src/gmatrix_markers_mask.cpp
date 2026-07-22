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
 File: gmatrix_markers_mask.cpp
 Description: implement binary marker mask

 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_MarkersMask(const arma::mat& M) {
  const uword n = M.n_rows, m = M.n_cols;
  // Build a binary mask B where B(i,j)=1 if !is.na(M(i,j)), else 0
  mat B(n, m, fill::zeros);
  for (uword j = 0; j < m; ++j)
    for (uword i = 0; i < n; ++i)
      if (!std::isnan(M(i,j))) B(i,j) = 1.0;
  
  return B * B.t();  // tcrossprod(mask, mask)
}
