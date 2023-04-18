// Casey Berger
// Created: Mar 9, 2023
// Last edited: Apr 18, 2023

#include <string> //string
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos, pow
#include <iomanip> //setw

#include "lattice.h"
#include "action_suite.h"

void check_phi_magnitude(double *** Lattice, int len); 
double phi_tot(double *** Lattice, int len); 
void print_lattice(double *** Lattice, int len);
void print_value(double *** Lattice, int i, int j, int len, double value, std::string valname);
void test_triangles(int i, int j, int len);
void test_QL(double *** Lattice, int i, int j, int len);