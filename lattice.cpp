// Casey Berger
// Created: Mar 28, 2023
// Last edited: Apr 18, 2023

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

void lattice_init(double *** Lattice, int len){

//Dynamically generates a 2D square lattice with a three-component phi at each site.
//See NRRB via CL code for alternative random number generator if you run into issues.
    
#ifdef TEST_CONSTANT_RN
    double r = 1.0;
#else
    srand(time(NULL)); //seed random number
    double r = ((double)rand())/((double)RAND_MAX);
#endif
    
#ifdef TESTING_MODE
    std:: cout << "Function: lattice_init in lattice" << std::endl;
#endif

    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            //generate a random polar and azimuthal angle
            double inclination =   M_PI * r; //polar angle = inclination
            double azimuth =  2. * M_PI * r; //azimuthal angle = azimuth
            //create unit spin vector components from angles
            Lattice[i][j][0] = sin(inclination) * cos(azimuth);
            Lattice[i][j][1] = sin(inclination) * sin(azimuth);
            Lattice[i][j][2] = cos(inclination);
        }
    }
    std::cout << "Lattice initialized" << std::endl;
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
    cross_prod[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
    cross_prod[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}