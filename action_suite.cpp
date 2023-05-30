// Casey Berger
// Created: Apr 18, 2023
// Last edited: May 17, 2023
//delete
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>

#include "action_suite.h"
#include "lattice.h"

/*
This file contains all the operations needed to calculate the 
lattice action, such as constructing the triangles, calculating
the topological charge Q_L and A_L
*/

void pick_phi(int i, int j, double (&phi)[3], double *** Lattice, bool old_lattice){
    if(old_lattice==true){
        //set phi
        phi[0] = Lattice[i][j][0];
        phi[1] = Lattice[i][j][1];
        phi[2] = Lattice[i][j][2];
    }
    else{
        //set phi
        phi[0] = Lattice[i][j][3];
        phi[1] = Lattice[i][j][4];
        phi[2] = Lattice[i][j][5];
    }
}