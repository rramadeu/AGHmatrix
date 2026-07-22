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
 File: check_matrix.cpp
 Description: implement check_matrix method

 Written by Thiago de Paula Oliveira, 2025
 */
// Thiago de Paula Oliveira

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void check_matrix(SEXP obj) {
  if (!Rf_isMatrix(obj)) {
    stop("Input must be a matrix.");
  }
}
