//#include <stdio.h>
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
    
#ifdef TEST_CONSTANT_RNG
    double r = 1.0;
#else
    srand(time(NULL)); //seed random number
    double r = ((double)rand())/((double)RAND_MAX);
#endif
    std::cout << "Initializing fields on lattice" << std::endl;
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