// Casey Berger
// Created: Mar 28, 2023
// Last edited: Apr 18, 2023

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>


#ifndef LATTICE_H
#define LATTICE_H

void generate_phi(double (&phi)[3]);
void lattice_init(double *** Lattice, int len);
void lattice_flush(double *** Lattice, int len);//may not even need this??
int plus_one(int i, int len);
int minus_one(int i, int len);
double dot_product(double vec1[3], double vec2[3]);
void cross_product(double vec1[3], double vec2[3], double (&cross_prod)[3]);
#endif
