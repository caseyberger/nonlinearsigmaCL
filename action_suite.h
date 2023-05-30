// Casey Berger
// Created: Apr 18, 2023
// Last edited: Apr 18, 2023
//delete
#include <stdlib.h>
#include <iostream>
#include <cstdlib>

#include "lattice.h"

#ifndef ACTION_H
#define ACTION_H

void pick_phi(int i, int j, double (&phi)[3], double *** Lattice, bool old_lattice);
#endif