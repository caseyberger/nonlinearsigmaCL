// Casey Berger
// Created: Apr 25, 2023
// Last edited: Apr 25, 2023

#include <time.h>
#include <math.h> //exp

#include "lattice.h"
#include "action_suite.h"
#include "monte_carlo.h"

/*
This file contains all the monte carlo operations, such as 
the boltzmann weight, accept-reject, etc.
*/

void Metropolis_loop(double beta, double itheta, double *** Lattice, int len){
    //to make the loop random, try shuffling the list of i and j coordinates: https://www.geeksforgeeks.org/shuffle-an-array-using-stl-in-c/
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: Metropolis loop in monte_carlo" << std::endl;
    std:: cout << "Note: Currently, looping over the lattice in sequential order. Update to make this a random order." << std::endl;
#endif
    double S_old, S_new, delta_S;
    double phi_new[3];
    bool old_lattice = true;
    bool new_lattice = false;
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
#ifdef TEST_CONSTANT_RN
            double r = 1.0;
#else
            srand(time(NULL)); //seed random number
            double r = ((double)rand())/((double)RAND_MAX);
#endif
            S_old = S_lattice(beta, Lattice, len, itheta, old_lattice);
            generate_phi(phi_new);
            Lattice[i][j][3] = phi_new[0];
            Lattice[i][j][4] = phi_new[1];
            Lattice[i][j][5] = phi_new[2];
            S_new = S_lattice(beta, Lattice, len, itheta, new_lattice);
            delta_S = S_new - S_old;
            if(delta_S < 0 || exp(delta_S) > r){
                Lattice[i][j][0] = Lattice[i][j][3];
                Lattice[i][j][1] = Lattice[i][j][4];
                Lattice[i][j][2] = Lattice[i][j][5];
            }//accept-reject step
        }//loop over j
    }//loop over i
}


