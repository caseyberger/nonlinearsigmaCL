// Casey Berger
// Created: Mar 28, 2023
// Last edited: May 2, 2023

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include "lattice.h"

/*
This file contains all the lattice admin operations, such as
initializing, saving, etc.
*/

void generate_phi(double (&phi)[3]){
#ifdef TEST_CONSTANT_RN
    double r = 0.5;
#else
    double r1 = ((double)rand())/((double)RAND_MAX);
    double r2 = ((double)rand())/((double)RAND_MAX);
#endif
    //generate a random polar and azimuthal angle
    double inclination =   M_PI * r1; //polar angle = inclination
    double azimuth =  2. * M_PI * r2; //azimuthal angle = azimuth
    //create unit spin vector components from angles
    phi[0] = sin(inclination) * cos(azimuth);
    phi[1] = sin(inclination) * sin(azimuth);
    phi[2] = cos(inclination);

    std:: cout << dot_product(phi, phi) << std::endl;
}

void lattice_init(double *** Lattice, int len){

//Dynamically generates a 2D square lattice with two three-component phis at each site (an old and a new)
//See NRRB via CL code for alternative random number generator if you run into issues.

#ifdef TESTING_MODE
    std:: cout << "Function: lattice_init in lattice" << std::endl;
#endif
    double phi[3];
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            generate_phi(phi);
            Lattice[i][j][0] = phi[0]; //holds current phi
            Lattice[i][j][1] = phi[1]; //holds current phi
            Lattice[i][j][2] = phi[2]; //holds current phi
            Lattice[i][j][3] = 0.0; //holds new phi
            Lattice[i][j][4] = 0.0; //holds new phi
            Lattice[i][j][5] = 0.0; //holds new phi
        }
    }
    std::cout << "Lattice initialized" << std::endl;
}

void lattice_flush(double *** Lattice, int len){
//updates old phi with new phi and sets new phi to zero
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            Lattice[i][j][0] = Lattice[i][j][3];//update
            Lattice[i][j][1] = Lattice[i][j][4];//update
            Lattice[i][j][2] = Lattice[i][j][5];//update
            Lattice[i][j][3] = 0.0; //reset
            Lattice[i][j][4] = 0.0; //reset
            Lattice[i][j][5] = 0.0; //reset
        }
    }
}

int plus_one(int i, int len){
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: plus_one in lattice" << std::endl;
#endif
    //returns site plus one, using periodic boundary conditions
    if (i==len-1){
        return 0;
    }
    else{
        return i+1;
    }
}

int minus_one(int i, int len){
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: minus_one in lattice" << std::endl;
#endif
    //returns site minus one, using periodic boundary conditions
    if (i==0){
        return len-1;
    }
    else{
        return i-1;
    }
}

double dot_product(double vec1[3], double vec2[3]){
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: dot_product in lattice" << std::endl;
#endif
    //calculates the dot product of two vectors
    double dot_prod = 0.0;
    dot_prod = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
    return dot_prod;
}

void cross_product(double vec1[3], double vec2[3],double (&cross_prod)[3]){
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: cross_product in lattice" << std::endl;
#endif
    //calculates the cross product of two vectors
    cross_prod[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    cross_prod[1] = - vec1[0]*vec2[2] + vec1[2]*vec2[0];
    cross_prod[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}
