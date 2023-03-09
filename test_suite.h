// Casey Berger
// Created: Mar 9, 2023
// Last edited: Mar 9, 2023

#include <string> //string
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos, pow
#include <iomanip> //setw

#include "test_suite.h"

void print_lattice(double Lattice[len][len][3]);
void print_value(double Lattice[len][len][3], double value[len][len]);
void test_triangles(int i, int j);
void test_QL(double QLcos, double QLsin);
