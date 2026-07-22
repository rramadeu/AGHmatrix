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
 File: amatrix_kerr.cpp
 Description: implement amatrix diploid method

 Written by Thiago de Paula Oliveira, 2025
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix buildA_ploidy2_cpp(IntegerVector sire, IntegerVector dire, int n) {
  NumericMatrix A(n, n);
  
  A(0, 0) = 1.0;
  
  for (int i = 1; i < n; ++i) {
    int si = sire[i] - 1; // 1-based to 0-based
    int di = dire[i] - 1;
    
    if (sire[i] == 0 && dire[i] == 0) {
      A(i, i) = 1.0;
    } else if (sire[i] == 0) {
      A(i, i) = 1.0;
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, di);
      }
    } else if (dire[i] == 0) {
      A(i, i) = 1.0;
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, si);
      }
    } else {
      A(i, i) = 1.0 + 0.5 * A(di, si);
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * (A(j, si) + A(j, di));
      }
    }
  }
  
  return A;
}
