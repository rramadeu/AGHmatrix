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
 File: dominance_matrix.cpp
 Description: implement pedigree dominance matrix

 Written by Thiago de Paula Oliveira, 2025
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix buildDominanceMatrix_cpp(NumericMatrix A, IntegerVector s, IntegerVector d) {
  int n = A.nrow();
  NumericMatrix D(n, n);
  
  for (int i = 0; i < n; ++i) {
    int si = s[i] - 1;
    int di = d[i] - 1;
    
    for (int j = i; j < n; ++j) {
      int sj = s[j] - 1;
      int dj = d[j] - 1;
      
      double u1 = (si >= 0 && sj >= 0) ? A(si, sj) : 0.0;
      double u2 = (di >= 0 && dj >= 0) ? A(di, dj) : 0.0;
      double u3 = (si >= 0 && dj >= 0) ? A(si, dj) : 0.0;
      double u4 = (sj >= 0 && di >= 0) ? A(sj, di) : 0.0;
      
      D(i, j) = D(j, i) = 0.25 * (u1 * u2 + u3 * u4);
    }
  }
  
  for (int i = 0; i < n; ++i) {
    D(i, i) = 1.0;
  }
  
  return D;
}
