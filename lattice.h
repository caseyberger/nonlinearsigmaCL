#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>


#ifndef LATTICE_H
#define LATTICE_H

void lattice_init(double *** Lattice, int len);
int plus_one(int i, int len);
int minus_one(int i, int len);

#endif