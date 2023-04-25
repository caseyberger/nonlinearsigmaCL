// Casey Berger
// Created: Apr 18, 2023
// Last edited: Apr 18, 2023
#include <stdlib.h>
#include <iostream>
#include <cstdlib>

#include "lattice.h"

#ifndef ACTION_H
#define ACTION_H

void pick_phi(int i, int j, int (&phi)[3], double *** Lattice, bool old_lattice);
void make_triangles(int i, int j, int len, int (&triangles)[8][3][2]);
double QL_triangle(int current_triangle[3][2], double *** Lattice, bool arcsin, bool old_lattice);
double Q_lattice(double *** Lattice, int len, bool old_lattice);
double A_lattice(double beta, double *** Lattice, int len, bool old_lattice);
double S_lattice(double beta, double *** Lattice, int len, double itheta, bool old_lattice);
#endif