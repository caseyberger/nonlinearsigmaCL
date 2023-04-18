// Casey Berger
// Created: Apr 18, 2023
// Last edited: Apr 18, 2023
#include <stdlib.h>
#include <iostream>
#include <cstdlib>

#include "lattice.h"

#ifndef ACTION_H
#define ACTION_H

void make_triangles(int i, int j, int len, int (&triangles)[8][3][2]);
double QL_triangle(int current_triangle[3][2], double *** Lattice, bool arcsin);
double Q_lattice(double *** Lattice, int len);
double A_lattice(double beta, double *** Lattice, int len);
#endif