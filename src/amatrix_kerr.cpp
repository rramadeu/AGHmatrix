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
 Description: implement amatrix kerr method

 Written by Thiago de Paula Oliveira, 2025
 */
// Thiago de Paula Oliveira
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix buildA_kerr_cpp(IntegerVector s, IntegerVector d, double w, int v) {
  int n = s.size();
  NumericMatrix A(n, n);
  A(0, 0) = 1.0 / (2.0 * v);
  
  for (int i = 1; i < n; ++i) {
    int si = s[i] - 1;
    int di = d[i] - 1;
    
    if (s[i] == 0 && d[i] == 0) {
      A(i, i) = 1.0 / (2.0 * v);
    } else if (s[i] == 0) {
      A(i, i) = (1 + (v-1)*w + ((v-1)*(1-w)*(v*A(di, di) + 0.5 - 1))/(2.0*v - 1)) / (2.0*v);
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, di);
      }
    } else if (d[i] == 0) {
      A(i, i) = (1 + (v-1)*w + ((v-1)*(1-w)*(v*A(si, si) + 0.5 - 1))/(2.0*v - 1)) / (2.0*v);
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * A(j, si);
      }
    } else {
      A(i, i) = (1 + (v-1)*w + ((v-1)*(1-w)*(v*A(di, di) + v*A(si, si) - 1))/(2.0*v - 1)) / (2.0*v) + 0.5 * A(di, si);
      for (int j = 0; j < i; ++j) {
        A(j, i) = A(i, j) = 0.5 * (A(j, si) + A(j, di));
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      A(i, j) *= (2.0 * v);
  
  return A;
}
