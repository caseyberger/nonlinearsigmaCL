// Casey Berger
// Created: Apr 25, 2023
// Last edited: Apr 25, 2023

#include <time.h>
#include <math.h> //exp


#include "lattice.h"

#ifndef MC_H
#define MC_H

void Metropolis_loop(double beta, double itheta, double *** Lattice, int len);

#endif