#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>


#ifndef LATTICE_H
#define LATTICE_H

void lattice_init(double *** Lattice, int len);
int plus_one(int i, int len);
int minus_one(int i, int len);
double dot_product(double vec1[3], double vec2[3]);
void cross_product(double vec1[3], double vec2[3], double (&cross_prod)[3]);
#endif