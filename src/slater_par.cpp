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
 File: amatrix_slater.cpp
 Description: implement slater par method

 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat slater_par_cpp(const arma::mat& X, const unsigned int ploidy){
  const uword n = X.n_rows, m = X.n_cols;
  const uword glev = static_cast<uword>(ploidy) + 1u;
  mat out(n, m * glev, fill::zeros);
  
  for (uword j = 0; j < m; ++j){
    const double* col = X.colptr(j);
    for (uword i = 0; i < n; ++i){
      const double v = col[i];
      if (!std::isnan(v)) {
        // accept only exact integers 0..ploidy (legacy semantics)
        const double vr = std::round(v);
        if (std::fabs(v - vr) < 1e-12) {
          const long dv = static_cast<long>(vr);
          if (dv >= 0 && dv <= static_cast<long>(ploidy)) {
            const uword cj = j * glev + static_cast<uword>(dv);
            out(i, cj) = 1.0;
          }
        }
      }
      // else: NA stays 0 (effectively treats NA as absence)
    }
  }
  return out;
}
