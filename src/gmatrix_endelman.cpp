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
 File: gmatrix_endelman.cpp
 Description: implement gmatrix based on Endelman method
 
 Written by Thiago de Paula Oliveira, 2025
 */

// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gmatrix_utils.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat Gmatrix_Endelman(const arma::mat& M, int ploidy) {
  if (ploidy != 4) stop("'Endelman' method is just implemented for ploidy=4");
  const uword n = M.n_rows, m = M.n_cols;
  
  vec freq = mean(M, 0).t() / 4.0; // colMeans / ploidy
  vec freq_comp = 1.0 - freq;
  const double SixPQ = 6.0 * dot(square(freq), square(freq_comp));
  
  mat SNP = M;
  rowvec f  = freq.t();
  rowvec f2 = square(f);
  
  mat rep1 = repmat(f2, n, 1);      // (Frequency[,1]^2) repeated by rows
  mat rep2 = repmat(f,  n, 1);      // Frequency[,1]
  
  SNP = 6.0 * rep1 - 3.0 * rep2 % M + 0.5 * M % (M - 1.0);
  
  return (SNP * SNP.t()) / SixPQ;
}
